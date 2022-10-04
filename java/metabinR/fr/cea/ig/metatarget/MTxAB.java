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

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
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

public class MTxAB implements MetaBin{
	private static String version;
	
	private static Dictionary dictionary = null;
	private static FastaManager frm = null;
	private StringBuilder sb = null;
	
	//List of output writers for each AB cluster
	private List<BufferedOutputStream> bos;
	//Number of sequences per AB cluster
	private List<AtomicInteger> spc;
	
	public String go(String[] args) throws Exception {
		version = new Date(Utils.classBuildTimeMillis(getClass())).toString();
		System.out.println("version MTxAB ="+version);
		
		int numOfThreads = 1;
		int excludeMin = 1;
		int excludeMax = 0;
    	int kAB = 10;
    	int numOfClustersAB = 3;
    	List<String> inputFastaFileNames = null;
    	String outputClustersFileNamePrefix = null;
    	boolean keepQualities = false;
    	boolean compressOut = false;
    	boolean dryRun = false;
    	
	    CommandLineParser parser = new BasicParser();
	    Options options = createOptions();
	    HelpFormatter formatter = new HelpFormatter();
	    CommandLine cmd = null;
	    try {
	        cmd = parser.parse( options, args );
	        
	    	if(args.length<1){
	    		throw new Exception();
	    	}
	        
	        if(cmd.hasOption("t")){
	        	numOfThreads = Integer.parseInt(cmd.getOptionValue("t"));
	        }
	        if(cmd.hasOption("i")){
	        	inputFastaFileNames = Arrays.asList(cmd.getOptionValues("i"));
	        }
	        if(cmd.hasOption("kAB")){
	        	kAB = Integer.parseInt(cmd.getOptionValue("kAB"));
	        }
	        if(cmd.hasOption("nAB")){
	    		numOfClustersAB = Integer.parseInt(cmd.getOptionValue("nAB"));
	    	}	
	        if(cmd.hasOption("oAB")){
	        	outputClustersFileNamePrefix = cmd.getOptionValue("oAB");
	        }
	    	if(cmd.hasOption("eMin")){
	    		excludeMin = Integer.parseInt(cmd.getOptionValue("eMin"));
	    	}
	    	if(cmd.hasOption("eMax")){
	    		excludeMax = Integer.parseInt(cmd.getOptionValue("eMax"));
	    	}
	    	if(cmd.hasOption("q")){
	    		keepQualities = true;
	    	}
	    	if(cmd.hasOption("z")){
	    		compressOut = true;
	    	}
	    	if(cmd.hasOption("d")){
	    		dryRun = true;
	    	}
	    }
	    catch(Exception exp ) {
	    	formatter.printHelp( "MTxAB", options, true);
	    	System.exit(0);
	    }
	    
	    //return String of assignments
	    sb = new StringBuilder();
	    
	    //Process sequences for the 1st time and count kmers for Abundance Based Binning 
	    processSequencesAB_count(numOfThreads, kAB, inputFastaFileNames, excludeMin, excludeMax);
	    removeMinKmers(excludeMin);
	    System.out.println(Utils.RAMInfo(Runtime.getRuntime()));
	    
	    //histogram of kmer counts
		TIntLongHashMap countsHistoAB = dictionary.getCountsHisto();
		if(!dryRun)
			save_histogram(countsHistoAB, outputClustersFileNamePrefix);
		
		//Create Poisson models for each Abundance Based Binning cluster
	    ClusterPoisson[] clusterPoissonsAB = VectorUtils.createABClusterPoissonsEMsync(numOfClustersAB, excludeMin, excludeMax, dictionary);
	    
	    //Create vector structure for each Abundance Based Binning cluster
	    ClusterVectorAB[] clusterVectorsAB = dictionary.createABClusterVectors(clusterPoissonsAB, excludeMin, excludeMax);
	    
	    //Dictionary is not needed anymore
	    dictionary.clear();
	    dictionary = null;
	    
	    //Create output writers for each Abundance Based Binning cluster
	    prepareOutputFiles(numOfClustersAB, outputClustersFileNamePrefix, keepQualities, compressOut, dryRun);
	    
		//Process sequences for the 2nd time, assign a cluster to each one and write it to the output with the corresponding writer
		processSequencesAB_assign(numOfThreads, kAB, inputFastaFileNames, clusterVectorsAB, keepQualities);
		for(ClusterVectorAB cv : clusterVectorsAB){
	    	cv.clear();
	    	cv = null;
	    }
	    clusterVectorsAB = null;
		
	    //Close output writers
	    closeOutputFiles(numOfClustersAB);
	    
	    System.out.println(Utils.RAMInfo(Runtime.getRuntime()));
	    
	    return sb.toString();
	}

	@SuppressWarnings("static-access")
	private Options createOptions(){
		Options options = new Options();
		Option t = OptionBuilder.withArgName("numOfThreads").withLongOpt("numOfThreads").hasArg().withDescription("Number of threads to use.").isRequired(false).create("t");
		options.addOption(t);
		Option k = OptionBuilder.withArgName("kMerSizeAB").withLongOpt("kMerSizeAB").hasArg().withDescription("k-mer length for AB.").isRequired(false).create("kAB");
		options.addOption(k);
		Option input = OptionBuilder.withArgName("input").withLongOpt("input").hasArgs(Integer.MAX_VALUE).withDescription("Input Fasta/q files paths.").isRequired(true).create("i");
		options.addOption(input);
		Option excludeMinOption = OptionBuilder.withArgName("eMin").withLongOpt("eMin").hasArg().withDescription("ExcludeMin kmers.").isRequired(false).create("eMin");
		options.addOption(excludeMinOption);
		Option excludeMaxOption = OptionBuilder.withArgName("eMax").withLongOpt("eMax").hasArg().withDescription("ExcludeMax kmers.").isRequired(false).create("eMax");
		options.addOption(excludeMaxOption);
		Option output_AB = OptionBuilder.withArgName("outputAB").withLongOpt("outputAB").hasArg().withDescription("Output Abundance Based Clusters files prefix.").isRequired(true).create("oAB");
		options.addOption(output_AB);
		Option numOfClusters = OptionBuilder.withArgName("numOfClustersAB").withLongOpt("numOfClustersAB").hasArg().withDescription("Number of Clusters for AB.").isRequired(false).create("nAB");
		options.addOption(numOfClusters);
		Option keepQualities = OptionBuilder.withArgName("q").withLongOpt("quality").hasArg(false).withDescription("Keep qualities.").isRequired(false).create("q");
		options.addOption(keepQualities);
		Option compressOut = OptionBuilder.withArgName("z").withLongOpt("gzip").hasArg(false).withDescription("Compress output.").isRequired(false).create("z");
		options.addOption(compressOut);
		Option dryRun = OptionBuilder.withArgName("d").withLongOpt("dry").hasArg(false).withDescription("Dry run.").isRequired(false).create("d");
		options.addOption(dryRun);
		
		return options;
	}
	
	private void processSequencesAB_count(int numOfThreads, int kAB, List<String> inputFastaFileNames, int excludeMin, int excludeMax){
		try{
			int cpus = Runtime.getRuntime().availableProcessors();
		    int usingThreads = (cpus < numOfThreads ? cpus : numOfThreads);
		    System.out.println("cpus="+cpus);
		    System.out.println("using="+usingThreads);

		    dictionary = new Dictionary(1024*1024, usingThreads, excludeMin, excludeMax);
		    
		    CountDownLatch startSignal = new CountDownLatch(1);
		    CountDownLatch doneSignal = new CountDownLatch(usingThreads + 1);
		    
		    System.out.println(Utils.time()+" START of AB Counting");
		    
			ExecutorService pool = Executors.newFixedThreadPool(usingThreads + 1);
			
			Map<Integer, SequenceProcessor> sequenceProcessors = new HashMap<Integer, SequenceProcessor>();
			if(frm!=null){
				frm.clear();
				frm = null;
			}
			frm = new FastaManager(false, inputFastaFileNames, startSignal, doneSignal);
			pool.execute(frm);
			
			SequenceProcessor.resetCounters();
			//Starting threads
			for(int i=0; i<usingThreads; i++){
				SequenceProcessor sp = new SequenceProcessor(this, dictionary, frm, 
						SequenceProcessor.MODE.AB_KMERCOUNT, kAB, startSignal, doneSignal, null);
				sequenceProcessors.put(sp.getId(), sp);
				pool.execute(sp);
			}
			
			doneSignal.await();
			pool.shutdown();
			
			System.out.println(Utils.time()+" END of AB Counting");
			System.out.println(Utils.time()+" Loaded sequences: "+SequenceProcessor.getSequenceCount().get());
		}
		catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	private void removeMinKmers(int excludeMin) {
		try {
			System.out.println(Utils.time()+" Total kmers(before remove):\t"+dictionary.getKmerCodes().size()+"\n");
			//remove kmers<minCount
			if(excludeMin>1){
				System.out.println(Utils.time()+" Removing kmer with global count < "+excludeMin);
				dictionary.removeAll(excludeMin);
			}
		} 
		catch (Exception e) {
			e.printStackTrace(System.err);
		}
	}
	
	private void save_histogram(TIntLongHashMap countsHistoAB, String outputClustersFileNamePrefix) {
		try {
			File f = new File(outputClustersFileNamePrefix+"__AB.histogram.tsv").getCanonicalFile();
			f.getParentFile().mkdirs();
			BufferedOutputStream bo = new BufferedOutputStream(new FileOutputStream(f));
			IOUtils.write("counts\tfrequency\n", bo, StandardCharsets.UTF_8);
			for ( TIntLongIterator it = countsHistoAB.iterator(); it.hasNext(); ) {
				it.advance();
				IOUtils.write(it.key()+"\t"+it.value()+"\n", bo, StandardCharsets.UTF_8);
			}
			bo.flush(); bo.close();
		}
		catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	private void prepareOutputFiles(int numOfClustersAB, String outputClustersFileNamePrefix,
			boolean keepQualities, boolean compressOut, boolean dryRun) {
		try {
			bos = Arrays.asList(new BufferedOutputStream[numOfClustersAB+1]);
			spc = new ArrayList<AtomicInteger>(numOfClustersAB);
			for(int i=0; i<numOfClustersAB; i++){
				if(!dryRun) {
					File f;
					if(compressOut) f = new File(outputClustersFileNamePrefix+"__AB."+(i+1)+(frm.isFastq&&keepQualities?".fastq.gz":".fasta.gz")).getCanonicalFile();
					else f = new File(outputClustersFileNamePrefix+"__AB."+(i+1)+(frm.isFastq&&keepQualities?".fastq":".fasta")).getCanonicalFile();
					f.getParentFile().mkdirs();
					BufferedOutputStream bo = null;
					if(compressOut) bo = new BufferedOutputStream(new GZIPOutputStream(new FileOutputStream(f)));
					else bo = new BufferedOutputStream(new FileOutputStream(f));
					bos.set(i,bo);
				}
				spc.add(new AtomicInteger(0));
			}
			//Create output writers for assignments
			List<Integer> numbers = Stream.iterate(1, n -> n + 1).limit(numOfClustersAB).collect(Collectors.toList());
			String line = "read_id\tAB\tAB."+
					numbers.stream().map(String::valueOf).collect(Collectors.joining("\tAB."))+
						"\n";
			sb.append(line);
			if(!dryRun) {
				File f = new File(outputClustersFileNamePrefix+"__AB.assignments.tsv").getCanonicalFile();
				f.getParentFile().mkdirs();
				bos.set(numOfClustersAB, new BufferedOutputStream(new FileOutputStream(f)));
				IOUtils.write(line, bos.get(numOfClustersAB), StandardCharsets.UTF_8);
			}
		} 
		catch (Exception e) {
			e.printStackTrace(System.err);
		}
	}
	
	private void processSequencesAB_assign(int numOfThreads, int kAB, List<String> inputFastaFileNames, ClusterVectorAB[] ABClusterVectors, 
			boolean keepQualities){
		try{
			int cpus = Runtime.getRuntime().availableProcessors();
		    int usingThreads = (cpus < numOfThreads ? cpus : numOfThreads);
		    System.out.println("cpus="+cpus);
		    System.out.println("using="+usingThreads);
		    
		    CountDownLatch startSignal = new CountDownLatch(1);
		    CountDownLatch doneSignal = new CountDownLatch(usingThreads+1);
			
			System.out.println(Utils.time()+" START of AB Binning");
			
			ExecutorService pool = Executors.newFixedThreadPool(usingThreads+1);
			
			Map<Integer, SequenceProcessor> sequenceProcessors = new HashMap<Integer, SequenceProcessor>();
			if(frm!=null){
				frm.clear();
				frm = null;
			}
			frm = new FastaManager(false, inputFastaFileNames, startSignal, doneSignal);
			pool.execute(frm);
			
			SequenceProcessor.resetCounters();
			//Starting threads
			for(int i=0; i<usingThreads; i++){
				SequenceProcessor sp = new SequenceProcessor(this, dictionary, frm, 
						SequenceProcessor.MODE.AB_BINNING, kAB, startSignal, doneSignal, new Object[]{keepQualities});
				sp.setABClusterVectors(ABClusterVectors);
				sequenceProcessors.put(sp.getId(), sp);
				pool.execute(sp);
			}
			
			doneSignal.await();
			pool.shutdown();
			
			System.out.println(Utils.time()+" END of AB Binning");
		}
		catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	private void closeOutputFiles(int numOfClustersAB) {
		try {
		    System.out.println("\tClustered reads:");
		    BufferedOutputStream bo;
			for(int i=0; i<numOfClustersAB; i++){
				bo = bos.get(i);
				if(bo!=null) {
					bo.flush();bo.close();
				}
				System.out.println("\t\tAB Cluster "+(i+1)+": "+spc.get(i).get());
			}
			bo = bos.get(numOfClustersAB);
			if(bo!=null) {
				bo.flush();bo.close();
			}
		} 
		catch (Exception e) {
			e.printStackTrace(System.err);
		}
	}
	
	@Override
	public synchronized boolean saveSeqToCluster(Sequence sequence, boolean keepQualities) {
		try {
			int clusterId = sequence.getAssignedCluster();
			if(clusterId <= 0){
				return false;
			}
			
			spc.get(clusterId-1).incrementAndGet();
			
			String line = sequence.getShortName()+"\t"+(clusterId)+"\t"+
					Arrays.stream(sequence.getDistancesToClusters()).
						mapToObj(String::valueOf).collect(Collectors.joining("\t"))+
							"\n";
			sb.append(line);
			if(bos.get(bos.size()-1)!=null)
				IOUtils.write(line, bos.get(bos.size()-1), StandardCharsets.UTF_8);
			
			BufferedOutputStream bo = bos.get(clusterId-1);
			if(bo==null) return false;
			
			if(sequence.getHeader()!=null){
				IOUtils.write(sequence.getHeader(), bo);
				IOUtils.write("\n", bo, StandardCharsets.UTF_8);
			}
			
			if(sequence.getSeq()!=null){
				IOUtils.write(sequence.getSeq(), bo);
				IOUtils.write("\n", bo, StandardCharsets.UTF_8);
			}
			
			if(keepQualities){
				if(sequence.getQual()!=null){
					IOUtils.write("+\n", bo, StandardCharsets.UTF_8);
					IOUtils.write(sequence.getQual(), bo);
					IOUtils.write("\n", bo, StandardCharsets.UTF_8);
				}
			}
		} 
		catch (Exception e) {
			e.printStackTrace(System.err);
			return false;
		}
		
		return true;
	}
	
	/*
	private void saveABClusters(ClusterVectorTrove[] ABClusterVectors, ArrayList<Integer> readAssignmentsAB, List<String> inputFastaFileNames, String outputClustersFileNamePrefix){
		try{
			//parse again the reads
			CountDownLatch startSignal = new CountDownLatch(1);
			CountDownLatch doneSignal = new CountDownLatch(1);
			ExecutorService pool = Executors.newFixedThreadPool(1);
			if(frm!=null){
				frm.clear();
				frm = null;
			}
			frm = new FastaManager(keepQualities, inputFastaFileNames, startSignal, doneSignal, true);
			pool.execute(frm);	
			doneSignal.await();
			pool.shutdown();
			
			System.out.println(Utils.time()+" START of Saving AB Clusters");
				
			ArrayList<BufferedOutputStream> bos = new ArrayList<BufferedOutputStream>();
			ArrayList<AtomicInteger> readsPerCluster = new ArrayList<AtomicInteger>();	
			for(int i=0; i<ABClusterVectors.length; i++){
				File f = new File(outputClustersFileNamePrefix+"__"+String.format("%03d", (i+1))+(frm.isFastq?".fastq":".fasta")).getCanonicalFile();
				f.getParentFile().mkdirs();
				BufferedOutputStream bo = new BufferedOutputStream(new FileOutputStream(f));
				bos.add(bo);
				readsPerCluster.add(new AtomicInteger(0));
			}
			
			List<Integer> sequenceIds = frm.getSequenceIds();
			for (int i=0; i<sequenceIds.size(); i++) {
				Sequence sequence = frm.getStaticSequences().get(sequenceIds.get(i));
				int cluster = readAssignmentsAB.get(i);
				if(cluster < 0){
					continue;
				}
				BufferedOutputStream bo = bos.get(cluster);
				readsPerCluster.get(cluster).incrementAndGet();
				byte[] header = sequence.getHeader();
				if(header!=null){
					IOUtils.write(header, bo);
					IOUtils.write("\n", bo, StandardCharsets.UTF_8);
				}
				byte[] seq = sequence.getSeq();
				if(seq!=null){
					IOUtils.write(seq, bo);
					IOUtils.write("\n", bo, StandardCharsets.UTF_8);
				}
				if(keepQualities){
					byte[] qual = sequence.getQual();
					if(qual!=null){
						IOUtils.write("+\n", bo, StandardCharsets.UTF_8);
						IOUtils.write(qual, bo);
						IOUtils.write("\n", bo, StandardCharsets.UTF_8);
					}
				}
				sequence.clearFull();
			}
			
			System.out.println(Utils.time()+" END of Saving AB Clusters");
			
			System.out.println("\tClassified reads:");
			for(int i=0; i<ABClusterVectors.length; i++){
				BufferedOutputStream bo = bos.get(i);
				bo.flush();bo.close();
				System.out.println("\t\tAB Cluster "+String.format("%03d", (i+1))+": "+readsPerCluster.get(i).get());
			}
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}
	*/
	
	public static void main(String[] args) {
		MTxAB MT_AB = new MTxAB();
		try {
			System.out.println(MT_AB.go(args));
		} 
		catch (Exception e) {
			e.printStackTrace(System.err);
		}
	}
	
}
