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

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.io.IOUtils;

import fr.cea.ig.metatarget.datastructures.FastaManager;
import fr.cea.ig.metatarget.datastructures.Sequence;
import fr.cea.ig.metatarget.datastructures.SequenceProcessor;
import fr.cea.ig.metatarget.kmeans.ClusterVectorCB;
import fr.cea.ig.metatarget.kmeans.ConcurrentKMeans;
import fr.cea.ig.metatarget.utils.Utils;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;

public class MTxCB implements MetaBin{
	private static String version;
	
	private static FastaManager frm = null;
	private StringBuilder sb = null;
	private TIntIntHashMap sequenceAssignmentsCB;
	private TIntObjectHashMap<double[]> sequenceDistancesCB;
	
	//List of output writers for each CB cluster
	private List<BufferedOutputStream> bos;
	//Number of sequences per CB cluster
	private List<AtomicInteger> spc;
	
	public String go(String[] args) throws Exception {
		version = new Date(Utils.classBuildTimeMillis(getClass())).toString();
		System.out.println("version MTxCB ="+version);
		
		int numOfThreads = 1;
    	int kCB = 4;
    	int numOfClustersCB = 5;
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
	        if(cmd.hasOption("kCB")){
	        	kCB = Integer.parseInt(cmd.getOptionValue("kCB"));
	        }
	        if(cmd.hasOption("nCB")){
	    		numOfClustersCB = Integer.parseInt(cmd.getOptionValue("nCB"));
	    	}
	        if(cmd.hasOption("i")){
	        	inputFastaFileNames = Arrays.asList(cmd.getOptionValues("i"));
	        }     
	        if(cmd.hasOption("oCB")){
	        	outputClustersFileNamePrefix = cmd.getOptionValue("oCB");
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
	    	formatter.printHelp( "MTxCB", options, true);
	    	System.exit(0);
	    }
	    
	    sequenceAssignmentsCB = new TIntIntHashMap();
	    sequenceDistancesCB = new TIntObjectHashMap<double[]>();
	    
	    //return String of assignments
	    sb = new StringBuilder();
	    
	    //Process sequences for the 1st time and create vector for each sequence
	    processSequencesCB_count(numOfThreads, kCB, inputFastaFileNames);
	    System.out.println(Utils.RAMInfo(Runtime.getRuntime()));
	    
	    //Create vector structure for each Composition Based Binning cluster
	    ClusterVectorCB[] clusterVectorsCB = createCBClusterVectors(numOfThreads, numOfClustersCB, kCB, frm.getStaticSequences());
	    parseAssignmentsCB(clusterVectorsCB);
	    
	    //Create output writers for each Composition Based Binning cluster
	    prepareOutputFiles(numOfClustersCB, outputClustersFileNamePrefix, keepQualities, compressOut, dryRun);
	    
		//Process sequences for the 2nd time and write them to the output with the corresponding writer
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
		Option kCB = OptionBuilder.withArgName("kMerSizeCB").withLongOpt("kMerSizeCB").hasArg().withDescription("k-mer length for CB.").isRequired(false).create("kCB");
		options.addOption(kCB);
		Option input = OptionBuilder.withArgName("input").withLongOpt("input").hasArgs(Integer.MAX_VALUE).withDescription("Input Fasta/q files paths.").isRequired(true).create("i");
		options.addOption(input);
		Option output_CB = OptionBuilder.withArgName("outputCB").withLongOpt("outputCB").hasArg().withDescription("Output Composition Based Clusters files prefix.").isRequired(true).create("oCB");
		options.addOption(output_CB);
		Option numOfClustersCB = OptionBuilder.withArgName("numOfClustersCB").withLongOpt("numOfClustersCB").hasArg().withDescription("Number of Clusters for CB.").isRequired(false).create("nCB");
		options.addOption(numOfClustersCB);
		Option keepQualities = OptionBuilder.withArgName("q").withLongOpt("quality").hasArg().withDescription("Keep qualities.").isRequired(false).create("q");
		options.addOption(keepQualities);
		Option compressOut = OptionBuilder.withArgName("z").withLongOpt("gzip").hasArg(false).withDescription("Compress output.").isRequired(false).create("z");
		options.addOption(compressOut);
		Option dryRun = OptionBuilder.withArgName("d").withLongOpt("dry").hasArg(false).withDescription("Dry run.").isRequired(false).create("d");
		options.addOption(dryRun);
		
		
		return options;
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
	
	private ClusterVectorCB[] createCBClusterVectors(int numOfThreads, int numOfClustersCB, int kCB, 
			List<Sequence> sequences){
		try{
			int cpus = Runtime.getRuntime().availableProcessors();
		    int usingThreads = (cpus < numOfThreads ? cpus : numOfThreads);
		    System.out.println("cpus="+cpus);
		    System.out.println("using="+usingThreads);
			CountDownLatch startSignal = new CountDownLatch(1);
			CountDownLatch doneSignal = new CountDownLatch(1);
			ExecutorService pool = Executors.newFixedThreadPool(1);
			//CountDownLatch doneSignal = new CountDownLatch(usingThreads);
			//ExecutorService pool = Executors.newFixedThreadPool(usingThreads);
			
			System.out.println(Utils.time()+" START of Creating CB Clusters\tSize="+sequences.size());
			ConcurrentKMeans kMeans = new ConcurrentKMeans(numOfClustersCB, kCB, sequences, 25, System.currentTimeMillis(), usingThreads, startSignal, doneSignal); //, dictionary);
			pool.execute(kMeans);	
			
			startSignal.countDown();
			doneSignal.await();
			pool.shutdown();
			
			System.out.println(Utils.time()+" END of Creating CB Clusters.");
			
			return kMeans.getClusters();
		}
		catch (Exception e) {
			e.printStackTrace(System.err);
			return null;
		}
	}
	
	private void parseAssignmentsCB(ClusterVectorCB[] CBClusters) {
		try {
			//Parse CB assignments
		    for(int i=0; i<CBClusters.length; i++){
				for(int seq_index : CBClusters[i].getMemberIndexes()) {
					sequenceAssignmentsCB.put(seq_index+1, (i+1));
				}
		    }
		    //Parse CB distances
		    for(Sequence s : frm.getStaticSequences()) {
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
			bos = Arrays.asList(new BufferedOutputStream[numOfClustersCB+1]);
			spc = new ArrayList<AtomicInteger>();
			for(int i=0; i<numOfClustersCB; i++){
				if(!dryRun) {
					File f;
					if(compressOut) f = new File(outputClustersFileNamePrefix+"__CB."+(i+1)+(frm.isFastq&&keepQualities?".fastq.gz":".fasta.gz")).getCanonicalFile();
					else f = new File(outputClustersFileNamePrefix+"__CB."+(i+1)+(frm.isFastq&&keepQualities?".fastq":".fasta")).getCanonicalFile();
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
			String line = "read_id\tCB\tCB."+
					numbers.stream().map(String::valueOf).collect(Collectors.joining("\tCB."))+
						"\n";
			sb.append(line);
			if(!dryRun) {
				File f = new File(outputClustersFileNamePrefix+"__CB.assignments.tsv").getCanonicalFile();
				f.getParentFile().mkdirs();
				bos.set(numOfClustersCB, new BufferedOutputStream(new FileOutputStream(f)));
				IOUtils.write(line, bos.get(numOfClustersCB), StandardCharsets.UTF_8);
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
				System.out.println("\t\tCB Cluster "+(i+1)+": "+spc.get(i).get());
			}
			bo = bos.get(numOfClustersCB);
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
			int clusterId = sequenceAssignmentsCB.get(sequence.getSequenceId());
			if(clusterId <= 0){
				return false;
			}
			
			spc.get(clusterId-1).incrementAndGet();
			
			String line = sequence.getShortName()+"\t"+(clusterId)+"\t"+
					Arrays.stream(sequenceDistancesCB.get(sequence.getSequenceId())).
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
	
	public static void main(String[] args) {
		MTxCB MT_CB = new MTxCB();
		try {
			System.out.println(MT_CB.go(args));
		} 
		catch (Exception e) {
			e.printStackTrace(System.err);
		}
	}
	
}
