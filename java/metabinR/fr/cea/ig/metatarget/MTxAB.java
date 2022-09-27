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

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.io.IOUtils;

import fr.cea.ig.metatarget.datastructures.ClusterPoisson;
import fr.cea.ig.metatarget.datastructures.ClusterVectorTrove;
import fr.cea.ig.metatarget.datastructures.Dictionary;
import fr.cea.ig.metatarget.datastructures.FastaManager;
import fr.cea.ig.metatarget.datastructures.Sequence;
import fr.cea.ig.metatarget.datastructures.SequenceProcessor;
import fr.cea.ig.metatarget.datastructures.VectorUtils;
import fr.cea.ig.metatarget.utils.Utils;
import gnu.trove.map.hash.TIntLongHashMap;

public class MTxAB {
	private static String version;
	
	private static int excludeMin = 1;
	private static int excludeMax = 0;
	
	private static Dictionary dictionary = null;
	private static FastaManager frm = null;

	private static boolean keepQualities = false;
	
	public void go(String[] args) {
		version = new Date(Utils.classBuildTimeMillis(getClass())).toString();
		System.out.println("version MTxAB ="+version);
		
		int numOfThreads = 1;
    	int kAB = 20;
    	int numOfClustersAB = 5;
    	List<String> inputFastaFileNames = null;
    	String outputClustersFileNamePrefix = null;
    	
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
	        if(cmd.hasOption("kAB")){
	        	kAB = Integer.parseInt(cmd.getOptionValue("kAB"));
	        }
	        if(cmd.hasOption("nAB")){
	    		numOfClustersAB = Integer.parseInt(cmd.getOptionValue("nAB"));
	    	}	
	        if(cmd.hasOption("i")){
	        	inputFastaFileNames = Arrays.asList(cmd.getOptionValues("i"));
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
	    		keepQualities = Boolean.parseBoolean(cmd.getOptionValue("q"));
	    	}
	    }
	    catch(Exception exp ) {
	    	formatter.printHelp( "MTxAB", options, true);
	    	System.exit(0);
	    }
	    
	    //Abundance based
	    processSequencesAB(numOfThreads, kAB, inputFastaFileNames);
	    @SuppressWarnings("unused")
		TIntLongHashMap countsHistoAB = dictionary.getCountsHisto();
	    ClusterPoisson[] clusterPoissons = VectorUtils.createABClusterPoissonsEMsync(numOfClustersAB, excludeMin, excludeMax, dictionary);
	    ClusterVectorTrove[] ABClusterVectors = dictionary.createABClusterVectors(clusterPoissons, excludeMin, excludeMax);
	    dictionary.clear();
	    dictionary = null;
	    
	    ArrayList<Integer> readAssignmentsAB = ABBinning(numOfThreads, kAB, inputFastaFileNames, ABClusterVectors);
	  
	    //Save
	    saveABClusters(ABClusterVectors, readAssignmentsAB, inputFastaFileNames, outputClustersFileNamePrefix);
	    
	    System.out.println(Utils.RAMInfo(Runtime.getRuntime()));
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
		Option keepQualities = OptionBuilder.withArgName("q").withLongOpt("quality").hasArg().withDescription("Keep qualities.").isRequired(false).create("q");
		options.addOption(keepQualities);
		
		return options;
	}
	
	private void processSequencesAB(int numOfThreads, int kAB, List<String> inputFastaFileNames){
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
			frm = new FastaManager(false, inputFastaFileNames, startSignal, doneSignal, false);
			pool.execute(frm);
			
			SequenceProcessor.resetCounters();
			//Starting threads
			for(int i=0; i<usingThreads; i++){
				SequenceProcessor sp = new SequenceProcessor(dictionary, frm, SequenceProcessor.MODE.KMERCOUNT, kAB, startSignal, doneSignal, null);
				sequenceProcessors.put(sp.getId(), sp);
				pool.execute(sp);
			}
			
			doneSignal.await();
			pool.shutdown();
			
			System.out.println(Utils.time()+" Total kmers(before remove):\t"+dictionary.getKmerCodes().size()+"\n");
			//remove kmers<minCount
			if(excludeMin>1){
				System.out.println(Utils.time()+" Removing kmer with global count < "+excludeMin);
				dictionary.removeAll(excludeMin);
			}
			
			System.out.println(Utils.time()+" END of AB Counting");
			System.out.println(Utils.time()+" Loaded sequences: "+SequenceProcessor.getSequenceCount().get());
		}
		catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	private ArrayList<Integer> ABBinning(int numOfThreads, int kAB, List<String> inputFastaFileNames, ClusterVectorTrove[] ABClusterVectors){
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
			frm = new FastaManager(false, inputFastaFileNames, startSignal, doneSignal, false);
			pool.execute(frm);
			
			SequenceProcessor.resetCounters();
			//Starting threads
			for(int i=0; i<usingThreads; i++){
				SequenceProcessor sp = new SequenceProcessor(dictionary, frm, SequenceProcessor.MODE.ABBINNING, kAB, startSignal, doneSignal, null);
				sp.setABClusterVectors(ABClusterVectors);
				sequenceProcessors.put(sp.getId(), sp);
				pool.execute(sp);
			}
			
			doneSignal.await();
			pool.shutdown();
			
			ArrayList<Integer> readAssignmentsAB = new ArrayList<Integer>();
			
			List<Integer> sequenceIds = frm.getSequenceIds();
			for (int i=0; i<sequenceIds.size(); i++) {	
				readAssignmentsAB.add(frm.getStaticSequences().get(sequenceIds.get(i)).getAbundanceCluster());
			}
			
			System.out.println(Utils.time()+" END of AB Binning");
			
			return readAssignmentsAB;
		}
		catch (Exception e) {
			e.printStackTrace();
			return null;
		}
	}
	
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
	
	public static void main(String[] args) {
		MTxAB MT_AB = new MTxAB();
		MT_AB.go(args);
	}
}
