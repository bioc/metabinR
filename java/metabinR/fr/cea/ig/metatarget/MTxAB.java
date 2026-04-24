/*
 *
 * MetaTarget MTxAB
 *
 * Copyright (C) 2022 Anestis Gkanogiannis <anestis@gkanogiannis.com>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 *
 */
package fr.cea.ig.metatarget;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.io.IOUtils;

import fr.cea.ig.metatarget.datastructures.ClusterPoisson;
import fr.cea.ig.metatarget.datastructures.ClusterVectorAB;
import fr.cea.ig.metatarget.datastructures.Dictionary;
import fr.cea.ig.metatarget.datastructures.FastaManager;
import fr.cea.ig.metatarget.datastructures.Sequence;
import fr.cea.ig.metatarget.datastructures.SequenceProcessor;
import fr.cea.ig.metatarget.datastructures.VectorUtils;
import fr.cea.ig.metatarget.utils.Utils;
import gnu.trove.iterator.TIntLongIterator;
import gnu.trove.map.hash.TIntLongHashMap;

public class MTxAB implements MetaBin {
    private static String version;

    private static Dictionary dictionary = null;
    private static FastaManager frm = null;
    private StringBuilder sb = null;

    // List of output writers for each AB cluster
    private List<BufferedOutputStream> bos;
    // Number of sequences per AB cluster
    private List<AtomicInteger> spc;

    /**
     * Typed entry point for the abundance-based binning pipeline. Called from R
     * via rJava.
     *
     * @param inputFastaFiles         input FASTA/FASTQ paths (gzip allowed)
     * @param excludeMin              exclude k-mers with global count &lt; this
     * @param excludeMax              exclude k-mers with global count &gt;= this
     *                                (0 = disabled)
     * @param kAB                     k-mer length for abundance binning
     * @param numOfClustersAB         number of AB clusters
     * @param outputPrefix            output file prefix (path); required
     * @param keepQualities           keep FASTQ qualities in output
     * @param dryRun                  do not write any output files
     * @param compressOut             gzip output files
     * @param numOfThreads            number of worker threads
     * @return tab-separated assignments (header + one line per read)
     */
    public String run(String[] inputFastaFiles, int excludeMin, int excludeMax,
                      int kAB, int numOfClustersAB, String outputPrefix,
                      boolean keepQualities, boolean dryRun, boolean compressOut,
                      int numOfThreads) throws Exception {
        version = new Date(Utils.classBuildTimeMillis(getClass())).toString();
        System.out.println("version MTxAB =" + version);

        if (inputFastaFiles == null || inputFastaFiles.length == 0) {
            throw new IllegalArgumentException("No input FASTA/FASTQ files provided.");
        }
        if (outputPrefix == null || outputPrefix.isEmpty()) {
            throw new IllegalArgumentException("outputPrefix must be a non-empty string.");
        }
        List<String> inputFastaFileNames = Arrays.asList(inputFastaFiles);

        // return String of assignments
        sb = new StringBuilder();

        // Process sequences for the 1st time and count kmers for Abundance Based Binning
        processSequencesAB_count(numOfThreads, kAB, inputFastaFileNames, excludeMin, excludeMax);
        removeMinKmers(excludeMin);
        System.out.println(Utils.RAMInfo(Runtime.getRuntime()));

        // histogram of kmer counts
        TIntLongHashMap countsHistoAB = dictionary.getCountsHisto();
        if (!dryRun)
            save_histogram(countsHistoAB, outputPrefix);

        // Create Poisson models for each Abundance Based Binning cluster
        ClusterPoisson[] clusterPoissonsAB = VectorUtils.createABClusterPoissonsEMsync(
                numOfClustersAB, excludeMin, excludeMax, dictionary);

        // Create vector structure for each Abundance Based Binning cluster
        ClusterVectorAB[] clusterVectorsAB = dictionary.createABClusterVectors(
                clusterPoissonsAB, excludeMin, excludeMax);

        // Dictionary is not needed anymore
        dictionary.clear();
        dictionary = null;

        // Create output writers for each Abundance Based Binning cluster
        prepareOutputFiles(numOfClustersAB, outputPrefix, keepQualities, compressOut, dryRun);

        // Process sequences for the 2nd time, assign a cluster to each one and write it to the output
        processSequencesAB_assign(numOfThreads, kAB, inputFastaFileNames, clusterVectorsAB, keepQualities);
        for (ClusterVectorAB cv : clusterVectorsAB) {
            cv.clear();
            cv = null;
        }
        clusterVectorsAB = null;

        // Close output writers
        closeOutputFiles(numOfClustersAB);

        System.out.println(Utils.RAMInfo(Runtime.getRuntime()));

        return sb.toString();
    }

    private void processSequencesAB_count(int numOfThreads, int kAB, List<String> inputFastaFileNames, int excludeMin, int excludeMax) {
        try {
            int cpus = Runtime.getRuntime().availableProcessors();
            int usingThreads = (cpus < numOfThreads ? cpus : numOfThreads);
            System.out.println("cpus=" + cpus);
            System.out.println("using=" + usingThreads);

            dictionary = new Dictionary(1024 * 1024, usingThreads, excludeMin, excludeMax);

            CountDownLatch startSignal = new CountDownLatch(1);
            CountDownLatch doneSignal = new CountDownLatch(usingThreads + 1);

            System.out.println(Utils.time() + " START of AB Counting");

            ExecutorService pool = Executors.newFixedThreadPool(usingThreads + 1);

            Map<Integer, SequenceProcessor> sequenceProcessors = new HashMap<Integer, SequenceProcessor>();
            if (frm != null) {
                frm.clear();
                frm = null;
            }
            frm = new FastaManager(false, inputFastaFileNames, startSignal, doneSignal);
            pool.execute(frm);

            SequenceProcessor.resetCounters();
            // Starting threads
            for (int i = 0; i < usingThreads; i++) {
                SequenceProcessor sp = new SequenceProcessor(this, dictionary, frm,
                        SequenceProcessor.MODE.AB_KMERCOUNT, kAB, startSignal, doneSignal, null);
                sequenceProcessors.put(sp.getId(), sp);
                pool.execute(sp);
            }

            doneSignal.await();
            pool.shutdown();

            System.out.println(Utils.time() + " END of AB Counting");
            System.out.println(Utils.time() + " Loaded sequences: " + SequenceProcessor.getSequenceCount().get());
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void removeMinKmers(int excludeMin) {
        try {
            System.out.println(Utils.time() + " Total kmers(before remove):\t" + dictionary.getKmerCodes().size() + "\n");
            // remove kmers<minCount
            if (excludeMin > 1) {
                System.out.println(Utils.time() + " Removing kmer with global count < " + excludeMin);
                dictionary.removeAll(excludeMin);
            }
        } catch (Exception e) {
            e.printStackTrace(System.err);
        }
    }

    private void save_histogram(TIntLongHashMap countsHistoAB, String outputClustersFileNamePrefix) {
        try {
            File f = new File(outputClustersFileNamePrefix + "__AB.histogram.tsv").getCanonicalFile();
            f.getParentFile().mkdirs();
            BufferedOutputStream bo = new BufferedOutputStream(new FileOutputStream(f));
            IOUtils.write("counts\tfrequency\n", bo, StandardCharsets.UTF_8);
            for (TIntLongIterator it = countsHistoAB.iterator(); it.hasNext();) {
                it.advance();
                IOUtils.write(it.key() + "\t" + it.value() + "\n", bo, StandardCharsets.UTF_8);
            }
            bo.flush();
            bo.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void prepareOutputFiles(int numOfClustersAB, String outputClustersFileNamePrefix,
            boolean keepQualities, boolean compressOut, boolean dryRun) {
        try {
            bos = Arrays.asList(new BufferedOutputStream[numOfClustersAB + 1]);
            spc = new ArrayList<AtomicInteger>(numOfClustersAB);
            for (int i = 0; i < numOfClustersAB; i++) {
                if (!dryRun) {
                    File f;
                    if (compressOut)
                        f = new File(outputClustersFileNamePrefix + "__AB." + (i + 1) + (frm.isFastq && keepQualities ? ".fastq.gz" : ".fasta.gz")).getCanonicalFile();
                    else
                        f = new File(outputClustersFileNamePrefix + "__AB." + (i + 1) + (frm.isFastq && keepQualities ? ".fastq" : ".fasta")).getCanonicalFile();
                    f.getParentFile().mkdirs();
                    BufferedOutputStream bo = null;
                    if (compressOut)
                        bo = new BufferedOutputStream(new GZIPOutputStream(new FileOutputStream(f)));
                    else
                        bo = new BufferedOutputStream(new FileOutputStream(f));
                    bos.set(i, bo);
                }
                spc.add(new AtomicInteger(0));
            }
            // Create output writers for assignments
            List<Integer> numbers = Stream.iterate(1, n -> n + 1).limit(numOfClustersAB).collect(Collectors.toList());
            String line = "read_id\tAB\tAB." +
                    numbers.stream().map(String::valueOf).collect(Collectors.joining("\tAB.")) +
                    "\n";
            sb.append(line);
            if (!dryRun) {
                File f = new File(outputClustersFileNamePrefix + "__AB.assignments.tsv").getCanonicalFile();
                f.getParentFile().mkdirs();
                bos.set(numOfClustersAB, new BufferedOutputStream(new FileOutputStream(f)));
                IOUtils.write(line, bos.get(numOfClustersAB), StandardCharsets.UTF_8);
            }
        } catch (Exception e) {
            e.printStackTrace(System.err);
        }
    }

    private void processSequencesAB_assign(int numOfThreads, int kAB, List<String> inputFastaFileNames, ClusterVectorAB[] ABClusterVectors,
            boolean keepQualities) {
        try {
            int cpus = Runtime.getRuntime().availableProcessors();
            int usingThreads = (cpus < numOfThreads ? cpus : numOfThreads);
            System.out.println("cpus=" + cpus);
            System.out.println("using=" + usingThreads);

            CountDownLatch startSignal = new CountDownLatch(1);
            CountDownLatch doneSignal = new CountDownLatch(usingThreads + 1);

            System.out.println(Utils.time() + " START of AB Binning");

            ExecutorService pool = Executors.newFixedThreadPool(usingThreads + 1);

            Map<Integer, SequenceProcessor> sequenceProcessors = new HashMap<Integer, SequenceProcessor>();
            if (frm != null) {
                frm.clear();
                frm = null;
            }
            frm = new FastaManager(false, inputFastaFileNames, startSignal, doneSignal);
            pool.execute(frm);

            SequenceProcessor.resetCounters();
            // Starting threads
            for (int i = 0; i < usingThreads; i++) {
                SequenceProcessor sp = new SequenceProcessor(this, dictionary, frm,
                        SequenceProcessor.MODE.AB_BINNING, kAB, startSignal, doneSignal, new Object[] { keepQualities });
                sp.setABClusterVectors(ABClusterVectors);
                sequenceProcessors.put(sp.getId(), sp);
                pool.execute(sp);
            }

            doneSignal.await();
            pool.shutdown();

            System.out.println(Utils.time() + " END of AB Binning");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void closeOutputFiles(int numOfClustersAB) {
        try {
            System.out.println("\tClustered reads:");
            BufferedOutputStream bo;
            for (int i = 0; i < numOfClustersAB; i++) {
                bo = bos.get(i);
                if (bo != null) {
                    bo.flush();
                    bo.close();
                }
                System.out.println("\t\tAB Cluster " + (i + 1) + ": " + spc.get(i).get());
            }
            bo = bos.get(numOfClustersAB);
            if (bo != null) {
                bo.flush();
                bo.close();
            }
        } catch (Exception e) {
            e.printStackTrace(System.err);
        }
    }

    @Override
    public synchronized boolean saveSeqToCluster(Sequence sequence, boolean keepQualities) {
        try {
            int clusterId = sequence.getAssignedCluster();
            if (clusterId <= 0) {
                return false;
            }

            spc.get(clusterId - 1).incrementAndGet();

            String line = sequence.getShortName() + "\t" + (clusterId) + "\t" +
                    Arrays.stream(sequence.getDistancesToClusters()).
                            mapToObj(String::valueOf).collect(Collectors.joining("\t")) +
                    "\n";
            sb.append(line);
            if (bos.get(bos.size() - 1) != null)
                IOUtils.write(line, bos.get(bos.size() - 1), StandardCharsets.UTF_8);

            BufferedOutputStream bo = bos.get(clusterId - 1);
            if (bo == null) return false;

            if (sequence.getHeader() != null) {
                IOUtils.write(sequence.getHeader(), bo);
                IOUtils.write("\n", bo, StandardCharsets.UTF_8);
            }

            if (sequence.getSeq() != null) {
                IOUtils.write(sequence.getSeq(), bo);
                IOUtils.write("\n", bo, StandardCharsets.UTF_8);
            }

            if (keepQualities) {
                if (sequence.getQual() != null) {
                    IOUtils.write("+\n", bo, StandardCharsets.UTF_8);
                    IOUtils.write(sequence.getQual(), bo);
                    IOUtils.write("\n", bo, StandardCharsets.UTF_8);
                }
            }
        } catch (Exception e) {
            e.printStackTrace(System.err);
            return false;
        }

        return true;
    }
}
