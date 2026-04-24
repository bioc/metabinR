/*
 *
 * MetaTarget MTxCB
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

import fr.cea.ig.metatarget.datastructures.FastaManager;
import fr.cea.ig.metatarget.datastructures.Sequence;
import fr.cea.ig.metatarget.datastructures.SequenceProcessor;
import fr.cea.ig.metatarget.kmeans.ClusterVectorCB;
import fr.cea.ig.metatarget.kmeans.ConcurrentKMeans;
import fr.cea.ig.metatarget.utils.Utils;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;

public class MTxCB implements MetaBin {
    private static String version;

    private static FastaManager frm = null;
    private StringBuilder sb = null;
    private TIntIntHashMap sequenceAssignmentsCB;
    private TIntObjectHashMap<double[]> sequenceDistancesCB;

    // List of output writers for each CB cluster
    private List<BufferedOutputStream> bos;
    // Number of sequences per CB cluster
    private List<AtomicInteger> spc;

    /**
     * Typed entry point for the composition-based binning pipeline. Called from R
     * via rJava.
     *
     * @param inputFastaFiles input FASTA/FASTQ paths (gzip allowed)
     * @param kCB             k-mer length for composition binning
     * @param numOfClustersCB number of CB clusters
     * @param outputPrefix    output file prefix (path); required
     * @param keepQualities   keep FASTQ qualities in output
     * @param dryRun          do not write any output files
     * @param compressOut     gzip output files
     * @param numOfThreads    number of worker threads
     * @return tab-separated assignments (header + one line per read)
     */
    public String run(String[] inputFastaFiles, int kCB, int numOfClustersCB,
                      String outputPrefix, boolean keepQualities, boolean dryRun,
                      boolean compressOut, int numOfThreads) throws Exception {
        version = new Date(Utils.classBuildTimeMillis(getClass())).toString();
        System.out.println("version MTxCB =" + version);

        if (inputFastaFiles == null || inputFastaFiles.length == 0) {
            throw new IllegalArgumentException("No input FASTA/FASTQ files provided.");
        }
        if (outputPrefix == null || outputPrefix.isEmpty()) {
            throw new IllegalArgumentException("outputPrefix must be a non-empty string.");
        }
        List<String> inputFastaFileNames = Arrays.asList(inputFastaFiles);

        sequenceAssignmentsCB = new TIntIntHashMap();
        sequenceDistancesCB = new TIntObjectHashMap<double[]>();

        // return String of assignments
        sb = new StringBuilder();

        // Process sequences for the 1st time and create vector for each sequence
        processSequencesCB_count(numOfThreads, kCB, inputFastaFileNames);
        System.out.println(Utils.RAMInfo(Runtime.getRuntime()));

        // Create vector structure for each Composition Based Binning cluster
        ClusterVectorCB[] clusterVectorsCB = createCBClusterVectors(numOfThreads, numOfClustersCB, kCB, frm.getStaticSequences());
        parseAssignmentsCB(clusterVectorsCB);

        // Create output writers for each Composition Based Binning cluster
        prepareOutputFiles(numOfClustersCB, outputPrefix, keepQualities, compressOut, dryRun);

        // Process sequences for the 2nd time and write them to the output with the corresponding writer
        processSequencesCB_assign(numOfThreads, inputFastaFileNames, keepQualities);

        // Close output writers
        closeOutputFiles(numOfClustersCB);

        System.out.println(Utils.RAMInfo(Runtime.getRuntime()));

        return sb.toString();
    }

    private void processSequencesCB_count(int numOfThreads, int kCB, List<String> inputFastaFileNames) {
        try {
            int cpus = Runtime.getRuntime().availableProcessors();
            int usingThreads = (cpus < numOfThreads ? cpus : numOfThreads);
            System.out.println("cpus=" + cpus);
            System.out.println("using=" + usingThreads);

            CountDownLatch startSignal = new CountDownLatch(1);
            CountDownLatch doneSignal = new CountDownLatch(usingThreads + 1);

            System.out.println(Utils.time() + " START of CB Counting");

            ExecutorService pool = Executors.newFixedThreadPool(usingThreads + 1);

            Map<Integer, SequenceProcessor> sequenceProcessors = new HashMap<Integer, SequenceProcessor>();
            if (frm != null) {
                frm.clear();
                frm = null;
            }
            frm = new FastaManager(false, inputFastaFileNames, startSignal, doneSignal);
            pool.execute(frm);

            SequenceProcessor.resetCounters();
            for (int i = 0; i < usingThreads; i++) {
                SequenceProcessor sp = new SequenceProcessor(this, null, frm,
                        SequenceProcessor.MODE.CB_SEQUENCEVECTORBUILD, kCB, startSignal, doneSignal, null);
                sequenceProcessors.put(sp.getId(), sp);
                pool.execute(sp);
            }

            doneSignal.await();
            pool.shutdown();

            System.out.println(Utils.time() + " END of CB Counting");
            System.out.println(Utils.time() + " Loaded sequences: " + SequenceProcessor.getSequenceCount().get());
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private ClusterVectorCB[] createCBClusterVectors(int numOfThreads, int numOfClustersCB, int kCB,
            List<Sequence> sequences) {
        try {
            int cpus = Runtime.getRuntime().availableProcessors();
            int usingThreads = (cpus < numOfThreads ? cpus : numOfThreads);
            System.out.println("cpus=" + cpus);
            System.out.println("using=" + usingThreads);
            CountDownLatch startSignal = new CountDownLatch(1);
            CountDownLatch doneSignal = new CountDownLatch(1);
            ExecutorService pool = Executors.newFixedThreadPool(1);

            System.out.println(Utils.time() + " START of Creating CB Clusters\tSize=" + sequences.size());
            ConcurrentKMeans kMeans = new ConcurrentKMeans(numOfClustersCB, kCB, sequences, 25, System.currentTimeMillis(), usingThreads, startSignal, doneSignal);
            pool.execute(kMeans);

            startSignal.countDown();
            doneSignal.await();
            pool.shutdown();

            System.out.println(Utils.time() + " END of Creating CB Clusters.");

            return kMeans.getClusters();
        } catch (Exception e) {
            e.printStackTrace(System.err);
            return null;
        }
    }

    private void parseAssignmentsCB(ClusterVectorCB[] CBClusters) {
        try {
            for (int i = 0; i < CBClusters.length; i++) {
                for (int seq_index : CBClusters[i].getMemberIndexes()) {
                    sequenceAssignmentsCB.put(seq_index + 1, (i + 1));
                }
            }
            for (Sequence s : frm.getStaticSequences()) {
                sequenceDistancesCB.put(s.getSequenceId(), s.getDistancesToClusters());
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void prepareOutputFiles(int numOfClustersCB, String outputClustersFileNamePrefix,
            boolean keepQualities, boolean compressOut, boolean dryRun) {
        try {
            bos = Arrays.asList(new BufferedOutputStream[numOfClustersCB + 1]);
            spc = new ArrayList<AtomicInteger>();
            for (int i = 0; i < numOfClustersCB; i++) {
                if (!dryRun) {
                    File f;
                    if (compressOut)
                        f = new File(outputClustersFileNamePrefix + "__CB." + (i + 1) + (frm.isFastq && keepQualities ? ".fastq.gz" : ".fasta.gz")).getCanonicalFile();
                    else
                        f = new File(outputClustersFileNamePrefix + "__CB." + (i + 1) + (frm.isFastq && keepQualities ? ".fastq" : ".fasta")).getCanonicalFile();
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
            List<Integer> numbers = Stream.iterate(1, n -> n + 1).limit(numOfClustersCB).collect(Collectors.toList());
            String line = "read_id\tCB\tCB." +
                    numbers.stream().map(String::valueOf).collect(Collectors.joining("\tCB.")) +
                    "\n";
            sb.append(line);
            if (!dryRun) {
                File f = new File(outputClustersFileNamePrefix + "__CB.assignments.tsv").getCanonicalFile();
                f.getParentFile().mkdirs();
                bos.set(numOfClustersCB, new BufferedOutputStream(new FileOutputStream(f)));
                IOUtils.write(line, bos.get(numOfClustersCB), StandardCharsets.UTF_8);
            }
        } catch (Exception e) {
            e.printStackTrace(System.err);
        }
    }

    private void processSequencesCB_assign(int numOfThreads, List<String> inputFastaFileNames, boolean keepQualities) {
        try {
            int cpus = Runtime.getRuntime().availableProcessors();
            int usingThreads = (cpus < numOfThreads ? cpus : numOfThreads);
            System.out.println("cpus=" + cpus);
            System.out.println("using=" + usingThreads);

            CountDownLatch startSignal = new CountDownLatch(1);
            CountDownLatch doneSignal = new CountDownLatch(usingThreads + 1);

            System.out.println(Utils.time() + " START of CB Binning");

            ExecutorService pool = Executors.newFixedThreadPool(usingThreads + 1);

            Map<Integer, SequenceProcessor> sequenceProcessors = new HashMap<Integer, SequenceProcessor>();
            if (frm != null) {
                frm.clear();
                frm = null;
            }
            frm = new FastaManager(false, inputFastaFileNames, startSignal, doneSignal);
            pool.execute(frm);

            SequenceProcessor.resetCounters();
            for (int i = 0; i < usingThreads; i++) {
                SequenceProcessor sp = new SequenceProcessor(this, null, frm,
                        SequenceProcessor.MODE.CB_BINNING, -1, startSignal, doneSignal, new Object[] { keepQualities });
                sequenceProcessors.put(sp.getId(), sp);
                pool.execute(sp);
            }

            doneSignal.await();
            pool.shutdown();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void closeOutputFiles(int numOfClustersCB) {
        try {
            System.out.println("\tClustered reads:");
            BufferedOutputStream bo;
            for (int i = 0; i < numOfClustersCB; i++) {
                bo = bos.get(i);
                if (bo != null) {
                    bo.flush();
                    bo.close();
                }
                System.out.println("\t\tCB Cluster " + (i + 1) + ": " + spc.get(i).get());
            }
            bo = bos.get(numOfClustersCB);
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
            int clusterId = sequenceAssignmentsCB.get(sequence.getSequenceId());
            if (clusterId <= 0) {
                return false;
            }

            spc.get(clusterId - 1).incrementAndGet();

            String line = sequence.getShortName() + "\t" + (clusterId) + "\t" +
                    Arrays.stream(sequenceDistancesCB.get(sequence.getSequenceId())).
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
