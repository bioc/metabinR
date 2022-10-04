/*
 *
 * MetaTarget MTxABxCB
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
import fr.cea.ig.metatarget.kmeans.ClusterVectorCB;
import fr.cea.ig.metatarget.kmeans.ConcurrentKMeans;
import fr.cea.ig.metatarget.utils.Utils;
import gnu.trove.iterator.TIntLongIterator;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntLongHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;

public class MTxABxCB implements MetaBin{
	private static String version;
	
	private static Dictionary dictionary;
	private FastaManager frm;
	private StringBuilder sb = null;
	private TIntIntHashMap sequenceAssignmentsAB;
	private TIntIntHashMap sequenceAssignmentsCB;
	private TIntObjectHashMap<double[]> sequenceDistancesCB;
	private AtomicInteger currCBid;
	private StringBuilder clustersDesign = null;
	private boolean finalStep = false;
	
	//List of output writers for each output cluster
	private List<BufferedOutputStream> bos;
	//Number of sequences per output cluster
	private List<AtomicInteger> spc;
	
	public String go(String[] args) throws Exception {
		version = new Date(Utils.classBuildTimeMillis(getClass())).toString();
		System.out.println("version MTxABxCB ="+version);
		
		int numOfThreads = 1;
		int excludeMin = 1;
		int excludeMax = 0;
    	int kAB = 10;
    	int kCB = 4;
    	int numOfClustersAB = 3;
    	int readLength = 400;
    	long genomeSize = (long)3.0E06;
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
	        if(cmd.hasOption("kCB")){
	        	kCB = Integer.parseInt(cmd.getOptionValue("kCB"));
	        }
	        if(cmd.hasOption("nAB")){
	    		numOfClustersAB = Integer.parseInt(cmd.getOptionValue("nAB"));
	    	}
	        if(cmd.hasOption("r")){
	    		readLength = Integer.parseInt(cmd.getOptionValue("r"));
	    	}
	        if(cmd.hasOption("g")){
	    		genomeSize = Integer.parseInt(cmd.getOptionValue("g"));
	    	}
	        if(cmd.hasOption("oC")){
	        	outputClustersFileNamePrefix = cmd.getOptionValue("oC");
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
	    	formatter.printHelp( "MTxABxCB", options, true);
	    	System.exit(0);
	    }
	    
	    sequenceAssignmentsAB = new TIntIntHashMap();
	    sequenceAssignmentsCB = new TIntIntHashMap();
	    sequenceDistancesCB = new TIntObjectHashMap<double[]>();
	    currCBid = new AtomicInteger(0);
	    clustersDesign = new StringBuilder();
	    clustersDesign.append("ABxCB\tAB\tsize\n");
	    
	    //return String of assignments
	    sb = new StringBuilder();
	    
	    //Abundance based
	    //Process sequences for the 1st time and count kmers for Abundance Based Binning 
	    processSequencesAB_count(numOfThreads, kAB, inputFastaFileNames, excludeMin, excludeMax);
	    removeMinKmers(excludeMin);
	    System.out.println(Utils.RAMInfo(Runtime.getRuntime()));
	    
	    //histogram of kmer counts
	  	TIntLongHashMap countsHistoAB = dictionary.getCountsHisto();
	  	readLength = guessSeqLength(kAB);
	  	if(!dryRun)
	  		save_histogram(countsHistoAB, outputClustersFileNamePrefix);
	  	
	  	//Create Poisson models for each Abundance Based Binning cluster
		ClusterPoisson[] clusterPoissonsAB = VectorUtils.createABClusterPoissonsEMsync(numOfClustersAB, excludeMin, excludeMax, dictionary);
		//Filter and merge AB poisson distributions
		System.out.println(Utils.time()+"\tFilter before="+clusterPoissonsAB.length);
	    ClusterPoisson[] filtered_clusterPoissonsAB = VectorUtils.filterClusterPoissons(clusterPoissonsAB);
	    numOfClustersAB = filtered_clusterPoissonsAB.length;
	    System.out.println(Utils.time()+"\tFilter after="+filtered_clusterPoissonsAB.length);
	    
	    //Create vector structure for each Abundance Based Binning cluster (after filter and merge)
	    ClusterVectorAB[] clusterVectorsAB = dictionary.createABClusterVectors(filtered_clusterPoissonsAB, excludeMin, excludeMax);
	    
	    //Dictionary is not needed anymore
	    dictionary.clear();
	    dictionary = null;
	    
	    //Process sequences for the 2nd time, assign an AB cluster to each one and put it on a list for all AB assignments
	    processSequencesAB_assign(numOfThreads, kAB, inputFastaFileNames, clusterVectorsAB, keepQualities);	
	    for(ClusterVectorAB cv : clusterVectorsAB){
	    	cv.clear();
	    	cv = null;
	    }
	    clusterVectorsAB = null;
	    
	    //Compositional based
	    //Process sequences for the 3rd time and create vector for each sequence
	    processSequencesCB_count(numOfThreads, kCB, inputFastaFileNames);
	    System.out.println(Utils.RAMInfo(Runtime.getRuntime()));
	    
	    //We have now again the sequences
	    //Create a list for each AB cluster and add the corresponding sequences
	    //Create a list for each AB cluster with the corresponding sequence id mapping (index on the list to sequenceId)
	    Object[] listsAB = getSequenceListsAB(numOfClustersAB);
	    @SuppressWarnings("unchecked")
		List<List<Sequence>> sequenceListsAB = (List<List<Sequence>>) listsAB[0];
	    @SuppressWarnings("unchecked")
		List<List<Integer>> sequenceIdsMappingsAB = (List<List<Integer>>) listsAB[1];
	    
	    //Run CB for each AB and assign each sequence to a CB bin
	    CBBinning(numOfThreads, kCB, sequenceListsAB, sequenceIdsMappingsAB, filtered_clusterPoissonsAB, readLength, genomeSize, outputClustersFileNamePrefix);
	    System.out.println(Utils.RAMInfo(Runtime.getRuntime()));
	    
	    //Create output writers for each Composition Based Binning cluster
	    int numOfClustersCB = currCBid.get();
	    prepareOutputFiles(numOfClustersCB, outputClustersFileNamePrefix, keepQualities, compressOut, dryRun);
	    
	    //Process sequences for the 4th time and write them to the output with the corresponding writer
	  	finalStep = true;
	    processSequencesCB_assign(numOfThreads, inputFastaFileNames, keepQualities);
	  	    
	    //Close output writers
	    closeOutputFiles(numOfClustersCB);
	    
	    System.out.println(Utils.RAMInfo(Runtime.getRuntime()));
	    
	    return sb.toString();
	}
	
	@SuppressWarnings("static-access")
	private Options createOptions(){
		Options options = new Options();
		Option t = OptionBuilder.withArgName("numOfThreads").withLongOpt("numOfThreads").hasArg().withDescription("Number of threads to use.").isRequired(false).create("t");
		options.addOption(t);
		Option kAB = OptionBuilder.withArgName("kMerSizeAB").withLongOpt("kMerSizeAB").hasArg().withDescription("k-mer length for AB.").isRequired(false).create("kAB");
		options.addOption(kAB);
		Option kCB = OptionBuilder.withArgName("kMerSizeCB").withLongOpt("kMerSizeCB").hasArg().withDescription("k-mer length for CB.").isRequired(false).create("kCB");
		options.addOption(kCB);
		Option input = OptionBuilder.withArgName("input").withLongOpt("input").hasArgs(Integer.MAX_VALUE).withDescription("Input Fasta/q files paths.").isRequired(true).create("i");
		options.addOption(input);
		Option excludeMinOption = OptionBuilder.withArgName("eMin").withLongOpt("eMin").hasArg().withDescription("ExcludeMin kmers.").isRequired(false).create("eMin");
		options.addOption(excludeMinOption);
		Option excludeMaxOption = OptionBuilder.withArgName("eMax").withLongOpt("eMax").hasArg().withDescription("ExcludeMax kmers.").isRequired(false).create("eMax");
		options.addOption(excludeMaxOption);
		Option output_C = OptionBuilder.withArgName("outputC").withLongOpt("outputC").hasArg().withDescription("Output Clusters files prefix.").isRequired(true).create("oC");
		options.addOption(output_C);
		Option numOfClustersAB = OptionBuilder.withArgName("numOfClustersAB").withLongOpt("numOfClustersAB").hasArg().withDescription("Number of Clusters for AB.").isRequired(false).create("nAB");
		options.addOption(numOfClustersAB);
		Option readLength = OptionBuilder.withArgName("readLength").withLongOpt("readLength").hasArg().withDescription("Read length.").isRequired(false).create("r");
		options.addOption(readLength);
		Option genomeSize = OptionBuilder.withArgName("genomeSize").withLongOpt("genomeSize").hasArg().withDescription("Average genome size.").isRequired(false).create("g");
		options.addOption(genomeSize);
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
	
	private int guessSeqLength(int kAB) {
		try {
			return (int)(dictionary.getTotalKmers()/(2.0*frm.getStaticSequences().size())) + (kAB-1) ;
		} 
		catch (Exception e) {
			e.printStackTrace(System.err);
			return 150;
		}
	}
	
	private void removeMinKmers(int excludeMin) {
		try {
			System.out.println(Utils.time()+" Distinct kmers(before remove):\t"+dictionary.getKmerCodes().size()+"\n");
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
			File f = new File(outputClustersFileNamePrefix+"__ABxCB.histogram.tsv").getCanonicalFile();
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
	
	private void processSequencesCB_count(int numOfThreads, int kCB, List<String> inputFastaFileNames){
		try{
			int cpus = Runtime.getRuntime().availableProcessors();
		    int usingThreads = (cpus < numOfThreads ? cpus : numOfThreads);
		    System.out.println("cpus="+cpus);
		    System.out.println("using="+usingThreads);
		    
		    CountDownLatch startSignal = new CountDownLatch(1);
		    CountDownLatch doneSignal = new CountDownLatch(usingThreads+1);
		    
		    System.out.println(Utils.time()+" START of CB Counting");
		    
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
				SequenceProcessor sp = new SequenceProcessor(this, null, frm, 
						SequenceProcessor.MODE.CB_SEQUENCEVECTORBUILD, kCB, startSignal, doneSignal, null);
				sequenceProcessors.put(sp.getId(), sp);
				pool.execute(sp);
			}
			
			doneSignal.await();
			pool.shutdown();
			
			System.out.println(Utils.time()+" END of CB Counting");
			System.out.println(Utils.time()+" Loaded sequences: "+SequenceProcessor.getSequenceCount().get());
		}
		catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	private Object[] getSequenceListsAB(int numOfClustersAB){
		try{
			List<List<Sequence>> sequenceLists = new ArrayList<List<Sequence>>();
			for(int i=0; i<numOfClustersAB; i++){
				sequenceLists.add(new ArrayList<Sequence>());
			}
			
			List<List<Integer>> sequenceIdsMappings = new ArrayList<List<Integer>>();
			for(int i=0; i<numOfClustersAB; i++){
				sequenceIdsMappings.add(new ArrayList<Integer>());
			}
			
			for(Sequence s : frm.getStaticSequences()) {
				s.setAssignedCluster(sequenceAssignmentsAB.get(s.getSequenceId()));
				if(s.getAssignedCluster() <= 0) {
					continue;
				}
				sequenceLists.get(s.getAssignedCluster()-1).add(s);
				sequenceIdsMappings.get(s.getAssignedCluster()-1).add(s.getSequenceId());
			}
			
			return new Object[] {sequenceLists, sequenceIdsMappings};
		}
		catch(Exception e){
			e.printStackTrace();
			return null;
		}
	}
	
	private void CBBinning(int numOfThreads, int kCB, 
			List<List<Sequence>> sequenceListsAB, 
			List<List<Integer>> sequenceIdsMappingsAB, 
			ClusterPoisson[] clusterPoissonsAB, int readLength, long genomeSize,
			String outputClustersFileNamePrefix){
		try{
			for(int i=0; i<sequenceListsAB.size(); i++){
				System.out.println(Utils.time()+" AB Cluster="+(i+1)+"\tSize="+sequenceListsAB.get(i).size());
			}
			System.out.println("");
			
			for(int i=0; i<sequenceListsAB.size(); i++){
				int cpus = Runtime.getRuntime().availableProcessors();
			    int usingThreads = (cpus < numOfThreads ? cpus : numOfThreads);
			    System.out.println("cpus="+cpus);
			    System.out.println("using="+usingThreads);
				CountDownLatch startSignal = new CountDownLatch(1);
				CountDownLatch doneSignal = new CountDownLatch(1);
			    ExecutorService pool = Executors.newFixedThreadPool(1);
				//CountDownLatch doneSignal = new CountDownLatch(usingThreads);
			    //ExecutorService pool = Executors.newFixedThreadPool(usingThreads);
			    
			    List<Sequence> sequences = sequenceListsAB.get(i);
			    List<Integer> sequencesIdsMapping = sequenceIdsMappingsAB.get(i);
				ClusterPoisson clusterPoisson = clusterPoissonsAB[i];
				int numOfClustersCB = guessNumOfClustersCB(i, clusterPoisson.getGenomeAbundance(), clusterPoisson.getGenomeLength(), sequences.size(), readLength, genomeSize);
				if(sequences.size() < 10){
			    	System.out.println(Utils.time()+" AB Cluster="+(i+1)+"\tSize="+sequences.size()+" quiting CB phase.");
			    	continue;
			    }
				System.out.println(Utils.time()+" START of Creating CB Clusters for AB Cluster="+(i+1)+"\tSize="+sequences.size());
				ConcurrentKMeans kMeans =  new ConcurrentKMeans(numOfClustersCB, kCB, sequences, 25, System.currentTimeMillis(), usingThreads, startSignal, doneSignal);
				pool.execute(kMeans);
				
				startSignal.countDown();
				doneSignal.await();
				pool.shutdown();
				
				System.out.println(Utils.time()+" END of Creating CB Clusters for AB Cluster="+(i+1));
			
				ClusterVectorCB[] clusterVectorsCB = kMeans.getClusters();
				//prosoxh edw, ta index tou kathe CB cluster einai diaforetika
				parseAssignmentsCB(clusterVectorsCB, sequencesIdsMapping, (i+1));
			}
		}
		catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	private int guessNumOfClustersCB(int ABid, double abundance, double EMLength, long size, long readLength, long genomeSize){
		double newLength = (size * readLength) / (abundance);
		
		int estimatedEMSpecies = (int) (EMLength/genomeSize);
		if(estimatedEMSpecies<1){
			estimatedEMSpecies = 1;
		}
		int estimatednewSpecies = (int) (newLength/genomeSize);
		if(estimatednewSpecies<1){
			estimatednewSpecies = 1;
		}
		
		//estimatedEMSpecies += 1;
		//estimatednewSpecies += 1;
		
		System.out.println(String.format("%5s","ABid")+"\t"+String.format("%8s","size")+"\t"+String.format("%12.12s","abundance")+"\t"+String.format("%10.10s","EMLength")+"\t"+String.format("%10s","newLength")+"\t"+String.format("%10s","EMspecies")+"\t"+String.format("%10s","newspecies"));
		System.out.println(String.format("%5s","-----")+"\t"+String.format("%8s","----")+"\t"+String.format("%12.12s","---------")+"\t"+String.format("%10.10s","--------")+"\t"+String.format("%10s","---------")+"\t"+String.format("%10s","---------")+"\t"+String.format("%10s","----------"));
		System.out.println(String.format("%5d",(ABid+1))+"\t"+String.format("%,8d",size)+"\t"+String.format("%12.5f",abundance)+"\t"+String.format("%,10d",(long)EMLength)+"\t"+String.format("%,10d",(long)newLength)+"\t"+String.format("%5d",estimatedEMSpecies)+"\t\t"+String.format("%5d",estimatednewSpecies));
	
		//return estimatednewSpecies;
		return estimatedEMSpecies;
	}
	
	private void parseAssignmentsCB(ClusterVectorCB[] CBClusters, List<Integer> sequencesIdsMapping, int currABid) {
		try {
			//Parse CB assignments
		    for(int i=0; i<CBClusters.length; i++){
		    	currCBid.incrementAndGet();
				for(int seq_index : CBClusters[i].getMemberIndexes()) {
					sequenceAssignmentsCB.put(sequencesIdsMapping.get(seq_index), currCBid.get());
				}
				clustersDesign.append(currCBid.get()+"\t"+currABid+"\t"+CBClusters[i].getMemberIndexes().length+"\n");
		    }
		    //Parse CB distances
		    for(Sequence s : frm.getStaticSequences()) {
		    	//Prosoxh edw, oi distances einai gia ta CB sto current AB
		    	sequenceDistancesCB.put(s.getSequenceId(), s.getDistancesToClusters());
		    }
		} 
		catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	private void prepareOutputFiles(int numOfClustersCB, String outputClustersFileNamePrefix, 
			boolean keepQualities, boolean compressOut, boolean dryRun) {
		try {
			bos = Arrays.asList(new BufferedOutputStream[numOfClustersCB+2]);
			spc = new ArrayList<AtomicInteger>();
			for(int i=0; i<numOfClustersCB; i++){
				if(!dryRun) {
					File f;
					if(compressOut) f = new File(outputClustersFileNamePrefix+"__ABxCB."+(i+1)+(frm.isFastq&&keepQualities?".fastq.gz":".fasta.gz")).getCanonicalFile();
					else f = new File(outputClustersFileNamePrefix+"__ABxCB."+(i+1)+(frm.isFastq&&keepQualities?".fastq":".fasta")).getCanonicalFile();
					f.getParentFile().mkdirs();
					BufferedOutputStream bo = null;
					if(compressOut) bo = new BufferedOutputStream(new GZIPOutputStream(new FileOutputStream(f)));
					else bo = new BufferedOutputStream(new FileOutputStream(f));
					bos.set(i,bo);
				}
				spc.add(new AtomicInteger(0));
			}
			//Create output writers for assignments
			List<Integer> numbers = Stream.iterate(1, n -> n + 1).limit(numOfClustersCB).collect(Collectors.toList());
			String line = "read_id\tABxCB\tABxCB."+
					numbers.stream().map(String::valueOf).collect(Collectors.joining("\tABxCB."))+
						"\n";
			sb.append(line);
			if(!dryRun) {
				File f = new File(outputClustersFileNamePrefix+"__ABxCB.assignments.tsv").getCanonicalFile();
				f.getParentFile().mkdirs();
				bos.set(numOfClustersCB, new BufferedOutputStream(new FileOutputStream(f)));
				IOUtils.write(line, bos.get(numOfClustersCB), StandardCharsets.UTF_8);
				
				f = new File(outputClustersFileNamePrefix+"__ABxCB.clustersDesign.tsv").getCanonicalFile();
				f.getParentFile().mkdirs();
				bos.set(numOfClustersCB+1, new BufferedOutputStream(new FileOutputStream(f)));
			}
		} 
		catch (Exception e) {
			e.printStackTrace(System.err);
		}
	}
	
	private void processSequencesCB_assign(int numOfThreads, List<String> inputFastaFileNames, boolean keepQualities){
		try{
			int cpus = Runtime.getRuntime().availableProcessors();
		    int usingThreads = (cpus < numOfThreads ? cpus : numOfThreads);
		    System.out.println("cpus="+cpus);
		    System.out.println("using="+usingThreads);
		    
		    CountDownLatch startSignal = new CountDownLatch(1);
		    CountDownLatch doneSignal = new CountDownLatch(usingThreads+1);
			
			System.out.println(Utils.time()+" START of CB Binning");
			
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
				SequenceProcessor sp = new SequenceProcessor(this, null, frm, 
						SequenceProcessor.MODE.CB_BINNING, -1, startSignal, doneSignal, new Object[]{keepQualities});
				sequenceProcessors.put(sp.getId(), sp);
				pool.execute(sp);
			}
			
			doneSignal.await();
			pool.shutdown();
			
			
			
			/*
			System.out.println(Utils.time()+" START of Saving CB Clusters.");
			for(int j=0; j<CBClusters.length; j++){
				ClusterVectorCB CBCluster = CBClusters[j];
				System.out.println("\tCB Cluster="+(j+1)+"\tSize="+CBCluster.getMemberIndexes().length);
				File f = new File(outputClustersFileNamePrefix+"_CBcluster_"+String.format("%03d", (j+1))+(frm.isFastq?".fastq":".fasta")).getCanonicalFile();
				f.getParentFile().mkdirs();
				BufferedOutputStream bo = new BufferedOutputStream(new FileOutputStream(f));
					
				for(int n : CBCluster.getMemberIndexes()){
					Sequence read = frm.getReads().get(n);
					byte[] header = read.getHeader();
					if(header!=null){
						IOUtils.write(header, bo);
						IOUtils.write("\n", bo);
					}
					byte[] seq = read.getSeq();
					if(seq!=null){
						IOUtils.write(seq, bo);
						IOUtils.write("\n", bo);
					}
					if(keepQualities){
						byte[] qual = read.getQual();
						if(qual!=null){
							IOUtils.write("+\n", bo);
							IOUtils.write(qual, bo);
							IOUtils.write("\n", bo);
						}
					}
					read.clearFull();
				}
					
				bo.flush();
				bo.close();
			}
			System.out.println(Utils.time()+" END of Saving CB Clusters.");
			*/
			
			/*
			System.out.println("\tClassified reads:");
			for(int j=0; j<CBClusters.length; j++){
				System.out.println("\t\tCluster "+String.format("%03d", j)+": "+CBClusters[j].getMemberIndexes().length);
			}
			*/
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}
	
	private void closeOutputFiles(int numOfClustersCB) {
		try {
			System.out.println("\tClustered reads:");
		    BufferedOutputStream bo;
			for(int i=0; i<numOfClustersCB; i++){
				bo = bos.get(i);
				if(bo!=null) {
					bo.flush();bo.close();
				}
				System.out.println("\t\tABxCB Cluster "+(i+1)+": "+spc.get(i).get());
			}
			bo = bos.get(numOfClustersCB);
			if(bo!=null) {
				bo.flush();bo.close();
			}
			bo = bos.get(numOfClustersCB+1);
			if(bo!=null) {
				IOUtils.write(clustersDesign.toString(), bo, StandardCharsets.UTF_8);
				bo.flush();bo.close();
			}
		} 
		catch (Exception e) {
			e.printStackTrace(System.err);
		}
	}
	
	public synchronized boolean saveSeqToCluster(Sequence sequence, boolean keepQualities) {
		try {
			//we are in the AB assignment step
			if(!finalStep) {
				sequenceAssignmentsAB.put(sequence.getSequenceId(),
						sequence.getAssignedCluster());
			}
			//we are in the final step (CB assignment is done and ready to write output fasta)
			else {
				int clusterId = sequenceAssignmentsCB.get(sequence.getSequenceId());
				if(clusterId <= 0){
					return false;
				}
				
				spc.get(clusterId-1).incrementAndGet();
				
				String[] distancesLine = new String[currCBid.get()];
				Arrays.fill(distancesLine, "NA");
				double[] distances = sequenceDistancesCB.get(sequence.getSequenceId());
				int indxMin = Utils.indexOfSmallest(distances);
				for(int i=0; i<distances.length; i++)
					distancesLine[clusterId-1-indxMin+i] = String.valueOf(distances[i]);
				
				String line = sequence.getShortName()+"\t"+(clusterId)+"\t"+
						Arrays.stream(distancesLine).
							collect(Collectors.joining("\t"))+
								"\n";
				sb.append(line);
				if(bos.get(bos.size()-2)!=null)
					IOUtils.write(line, bos.get(bos.size()-2), StandardCharsets.UTF_8);
				
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
		}
		catch (Exception e) {
			e.printStackTrace(System.err);
			return false;
		}
		
		return true;
	}
	
	public static void main(String[] args) {
		MTxABxCB MT_AB_CB = new MTxABxCB();
		try {
			System.out.println(MT_AB_CB.go(args));
		} 
		catch (Exception e) {
			e.printStackTrace(System.err);
		}
	}
	
}
