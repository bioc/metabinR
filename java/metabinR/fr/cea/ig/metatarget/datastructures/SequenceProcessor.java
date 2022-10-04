/*
 *
 * MetaTarget SequenceProcessor
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
package fr.cea.ig.metatarget.datastructures;

import java.util.concurrent.CountDownLatch;
import java.util.concurrent.atomic.AtomicInteger;

import fr.cea.ig.metatarget.MetaBin;
import fr.cea.ig.metatarget.utils.Coder;
import fr.cea.ig.metatarget.utils.Utils;

public class SequenceProcessor implements Runnable {

	private static AtomicInteger sequenceCount = new AtomicInteger(0);
	
	private static AtomicInteger taskCount = new AtomicInteger(0);
	
	private final int id = taskCount.getAndIncrement();
	private MetaBin caller;
	private Dictionary dictionary = null;
	private FastaManager frm = null;
	private MODE mode = null;
	private int k;
	private CountDownLatch startSignal = null;
	private CountDownLatch doneSignal = null;
	private Object[] parametersExtra = null;
	
	private ClusterVectorAB[] ABClusterVectors = null;

	public static enum MODE {
		AB_KMERCOUNT,
		AB_BINNING,
		CB_SEQUENCEVECTORBUILD,
		CB_BINNING,
		KMERCOUNT_SEQUENCEVECTORBUILD
	}
	
	public SequenceProcessor(MetaBin caller, Dictionary dictionary, FastaManager frm, MODE mode, int k, CountDownLatch startSignal, CountDownLatch doneSignal, Object[] parametersExtra) {
		this.caller = caller;
		this.dictionary = dictionary;
		this.frm = frm;
		this.mode = mode;
		this.k = k;
		this.startSignal = startSignal;
		this.doneSignal = doneSignal;
		this.parametersExtra = parametersExtra;
	}
	
	public static void resetCounters(){
		sequenceCount = new AtomicInteger(0);
		taskCount = new AtomicInteger(0);
	}
	
	public static AtomicInteger getSequenceCount() {
		return sequenceCount;
	}

	public MODE getMode() {
		return mode;
	}

	public void setABClusterVectors(ClusterVectorAB[] ABClusterVectors) {
		this.ABClusterVectors = ABClusterVectors;
	}

	public void setStartSignal(CountDownLatch startSignal) {
		this.startSignal = startSignal;
	}

	public void setDoneSignal(CountDownLatch doneSignal) {
		this.doneSignal = doneSignal;
	}

	public void setMode(MODE mode) {
		this.mode = mode;
	}

	public int getId() {
		return id;
	}

	@Override
	public void run() {
		try{
			startSignal.await();
			boolean done = false;
			System.out.println(Utils.time()+" SequenceProcessor: "+id+"\t"+mode.toString()+" START");
			while(!done){
				Sequence sequence = frm.getNextSequence();
				if(sequence==null){
					if(!frm.hasMore()){
						System.out.println(Utils.time()+" SequenceProcessor: "+id+"\t"+mode.toString()+" EXIT");
						done = true;
						break;
					}
					continue;
				}
//				SequenceD2 sequenceD2 = null;
//				if(sequence!=null){
//					sequenceD2 = new SequenceD2(sequence);
//				}
//				if(sequenceD2==null){
//					if(!frm.hasMore()){
//						done = true;
//						break;
//					}
//					continue;
//				}

				switch(mode){
					case AB_KMERCOUNT :
						processSequence_AB_KMERCOUNT(sequence);
						break;
					case AB_BINNING :
						processSequence_AB_BINNING(sequence);
						break;
					case CB_SEQUENCEVECTORBUILD :
						processSequence_CB_SEQUENCEVECTORBUILD(sequence);
						break;
					case CB_BINNING :
						processSequence_CB_BINNING(sequence);
						break;
					case KMERCOUNT_SEQUENCEVECTORBUILD :
						System.out.println("UNSUPPORTED");
						break;
					default:
						System.err.println(Utils.time()+" Mode+\"" + mode.toString() + "\" not known.");
						break;
				}
			}
			doneSignal.countDown();
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}
	
	private void processSequence_AB_KMERCOUNT(Sequence sequence){
		try{
			sequenceCount.incrementAndGet();
			if(sequenceCount.get() % 100000 == 0){
				System.out.println(Utils.time()+" SequenceProcessor: "+id+"\t"+mode.toString()+" "+sequenceCount.get());
			}
			
			//long kmerCode;
			
			for(int i=0; i<=sequence.getLength()-k; i++){				
				//kmerCode = Coder.encodeToLong(sequence, i, i+k, false);
				dictionary.insert(Coder.encodeToLong(sequence, i, i+k, false));
				
				//reverse
				//kmerCode = Coder.encodeToLong(sequence, i, i+k, true);
				dictionary.insert(Coder.encodeToLong(sequence, i, i+k, true));
			}
			
			sequence.clearFull();
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}
	
	private void processSequence_AB_BINNING(Sequence sequence){
		try{
			sequenceCount.incrementAndGet();
			if(sequenceCount.get() % 100000 == 0){
				System.out.println(Utils.time()+" SequenceProcessor: "+id+"\t"+mode.toString()+" "+sequenceCount.get());
			}
			
			//long kmerCode;
			
			for(int i=0; i+k<=sequence.getLength(); i++){				
				//kmerCode = Coder.encodeToLong(sequence, i, i+k, false);
				sequence.insertKmerCount(Coder.encodeToLong(sequence, i, i+k, false), 1);
				
				//reverse
				//kmerCode = Coder.encodeToLong(sequence, i, i+k, true);
				sequence.insertKmerCount(Coder.encodeToLong(sequence, i, i+k, true), 1);			
			}
			
			//Double maxSim = Double.NEGATIVE_INFINITY;
			double maxSim = 0.0;
			int abundanceCluster = -1;
			double similarity;
			double[] distances = new double[ABClusterVectors.length];
			for(int c=0; c<ABClusterVectors.length; c++){
				similarity = VectorUtils.cosineSequence2Cluster(sequence, ABClusterVectors[c]);
				if( similarity > maxSim){
					maxSim = similarity;
					abundanceCluster = (c+1);
				}
				distances[c] = similarity;
			}
			sequence.setAssignedCluster(abundanceCluster);
			sequence.setDistancesToClusters(distances);
			
			if(parametersExtra!=null && ((Boolean)parametersExtra[0])==true){
				caller.saveSeqToCluster(sequence, true);
			}
			else {
				caller.saveSeqToCluster(sequence, false);
			}
			
			sequence.clearFull();
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}
	
	private void processSequence_CB_SEQUENCEVECTORBUILD(Sequence sequence){
		try{
			sequenceCount.incrementAndGet();
			if(sequenceCount.get() % 100000 == 0){
				System.out.println(Utils.time()+" SequenceProcessor: "+id+"\t"+mode.toString()+" "+sequenceCount.get());
			}
			
			//long kmerCode;
			
			for(int i=0; i+k<=sequence.getLength(); i++){				
				//kmerCode = Coder.encodeToLong(sequence, i, i+k, false);
				sequence.insertKmerCount(Coder.encodeToLong(sequence, i, i+k, false), 1);
				
				//reverse
				//kmerCode = Coder.encodeToLong(sequence, i, i+k, true);
				sequence.insertKmerCount(Coder.encodeToLong(sequence, i, i+k, true), 1);
			}
			
			//(boolean)parametersExtra[0] true, clear header and sequence
			//(boolean)parametersExtra[0] false, keep header and sequence
			//if(parametersExtra!=null && ((Boolean)parametersExtra[0])==true){
				sequence.clearHeadSeq();
			//}
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}
	
	private void processSequence_CB_BINNING(Sequence sequence){
		try{
			sequenceCount.incrementAndGet();
			if(sequenceCount.get() % 100000 == 0){
				System.out.println(Utils.time()+" SequenceProcessor: "+id+"\t"+mode.toString()+" "+sequenceCount.get());
			}
			
			if(parametersExtra!=null && ((Boolean)parametersExtra[0])==true){
				caller.saveSeqToCluster(sequence, true);
			}
			else {
				caller.saveSeqToCluster(sequence, false);
			}
			
			sequence.clearFull();
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}
	
	/*
	private void processSequence_KMERCOUNT_SEQUENCEVECTORBUILD(Sequence sequence){
		try{
			sequenceCount.incrementAndGet(); 
			if(sequenceCount.get() % 1000000 == 0){
				System.out.println(Utils.time()+" SequenceProcessor: "+id+"\t"+mode.toString()+" "+(sequence.getSequenceId()));
			}
			
			sequence.initKmerCounts();
			
			long kmerCode;
			
			for(int i=0; i+k<=sequence.getLength(); i++){				
				kmerCode = Coder.encodeToLong(sequence, i, i+k, false);
				if(kmerCode>=0L){
					dictionary.insert(kmerCode);
					sequence.insertKmer(kmerCode, 1);
				}
				
				//reverse
				kmerCode = Coder.encodeToLong(sequence, i, i+k, true);
				if(kmerCode>=0L){
					dictionary.insert(kmerCode);
					sequence.insertKmer(kmerCode, 1);
				}
			}
			
			//(boolean)parametersSequenceVectorClear[0] true, clear header and sequence
			//(boolean)parametersSequenceVectorClear[0] false, keep header and sequence
			if(parametersSequenceVectorClear!=null && ((Boolean)parametersSequenceVectorClear[0])==true){
				sequence.clearHeadSeq();
			}
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}
	*/
	
}
