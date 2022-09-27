/*
 *
 * MetaTarget SequenceProcessorD2
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
package fr.cea.ig.metatarget.datastructures.d2;

import java.util.concurrent.CountDownLatch;
import java.util.concurrent.atomic.AtomicInteger;

import fr.cea.ig.metatarget.datastructures.Dictionary;
import fr.cea.ig.metatarget.datastructures.DictionaryFilter;
import fr.cea.ig.metatarget.datastructures.FastaManager;
import fr.cea.ig.metatarget.datastructures.Sequence;
import fr.cea.ig.metatarget.utils.Coder;
import fr.cea.ig.metatarget.utils.Utils;

public class SequenceProcessorD2 implements Runnable {

	private static AtomicInteger sequenceCount = new AtomicInteger(0);
	
	private static AtomicInteger taskCount = new AtomicInteger(0);
	
	private final int id = taskCount.getAndIncrement();
	private Dictionary dictionary = null;
	private FastaManager frm = null;
	//For filtering with long kFilter
	private DictionaryFilter dictionaryFilter = null;
	private MODE mode = null;
	private int k;
	//For filtering with long kFilter
	private int kFilter = 0;
	private CountDownLatch startSignal = null;
	private CountDownLatch doneSignal = null;
	
	private SequenceD2Centroid rpSampleVector = null;
	
	public static enum MODE {
		KMER_COUNTING_SEQUENCE_D2,
		KMER_COUNTING_SAMPLE_D2
	}
	
	public SequenceD2Centroid getrpSampleVector(){
		return rpSampleVector;
	}
	
	public static void resetCounters(){
		sequenceCount = new AtomicInteger(0);
		taskCount = new AtomicInteger(0);
	}
	
	public SequenceProcessorD2(Dictionary dictionary, FastaManager frm, DictionaryFilter dictionaryFilter, MODE mode, int k, int kFilter, CountDownLatch startSignal, CountDownLatch doneSignal) {
		this.dictionary = dictionary;
		this.frm = frm;
		this.dictionaryFilter = dictionaryFilter;
		this.mode = mode;
		this.k = k;
		this.kFilter = kFilter;
		this.startSignal = startSignal;
		this.doneSignal = doneSignal;
	}
	
	public static AtomicInteger getSequenceCount() {
		return sequenceCount;
	}

	public MODE getMode() {
		return mode;
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
			if(rpSampleVector==null){
				rpSampleVector = new SequenceD2Centroid(new Sequence(id, null, null, true));
			}
			
			startSignal.await();
			boolean done = false;
			System.out.println(Utils.time()+" SequenceProcessorD2: "+id+"\t"+mode.toString()+" START");
			while(!done){
				Sequence sequence = frm.getNextSequence();
				SequenceD2 sequenceD2 = null;
				if(sequence!=null){
					sequenceD2 = new SequenceD2(sequence);
				}
				if(sequenceD2==null) {
					if(!frm.hasMore()){
						done = true;
						break;
					}
					continue;
				}
				switch(mode){
					case KMER_COUNTING_SEQUENCE_D2 :
						processSequence_KMER_COUNTING_SEQUENCE_D2(sequenceD2);
						break;
					case KMER_COUNTING_SAMPLE_D2 :
						processSequence_KMER_COUNTING_SAMPLE_D2(sequenceD2);
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
	
	
	private void processSequence_KMER_COUNTING_SEQUENCE_D2(SequenceD2 sequence){
		try{
			sequenceCount.incrementAndGet(); 
			if(sequenceCount.get() % 1000000 == 0){
				System.out.println(Utils.time()+" SequenceProcessorD2: "+id+"\t"+mode.toString()+" "+sequenceCount.get());
			}
			
			long kmerCode;
			long kmerCodeFilter;	
			
			for(int i=0; i+k<=sequence.getLength(); i++){	
				int oldAs = sequence.as;
				int oldTs = sequence.ts;
				int oldCs = sequence.cs;
				int oldGs = sequence.gs;
				kmerCode = Coder.encodeToLong(sequence, i, i+k, false, true);
				if(kmerCode>=0L){
					if(dictionary!=null){
						dictionary.insert(kmerCode);
					}
					sequence.insertKmerCount(kmerCode, 1);
					sequence.insertKmerProb(kmerCode, (short)(sequence.as-oldAs), (short)(sequence.ts-oldTs), (short)(sequence.cs-oldCs), (short)(sequence.gs-oldGs));
				}
				
				//reverse
				oldAs = sequence.as;
				oldTs = sequence.ts;
				oldCs = sequence.cs;
				oldGs = sequence.gs;
				kmerCode = Coder.encodeToLong(sequence, i, i+k, true, true);
				if(kmerCode>=0L){
					if(dictionary!=null){
						dictionary.insert(kmerCode);
					}
					sequence.insertKmerCount(kmerCode, 1);
					sequence.insertKmerProb(kmerCode, (short)(sequence.as-oldAs), (short)(sequence.ts-oldTs), (short)(sequence.cs-oldCs), (short)(sequence.gs-oldGs));
				}
			}
		
			//calculate per sequence probs
			sequence.calculateProbs(k);
			rpSampleVector.addWith((SequenceD2Interface)sequence, true);
			sequence.clear();
			

			/////////////////////
			//long kmer fitlering
			if(kFilter>0 && dictionaryFilter!=null){
				for(int i=0; i+kFilter<=sequence.getLength(); i++){	
					kmerCodeFilter = Coder.encodeToLong(sequence, i, i+kFilter, false, false);
					if(kmerCodeFilter>=0L){
						dictionaryFilter.insert(kmerCodeFilter, sequence.getSequenceId());
					}
					
					//reverse
					kmerCodeFilter = Coder.encodeToLong(sequence, i, i+kFilter, true, false);
					if(kmerCodeFilter>=0L){
						dictionaryFilter.insert(kmerCodeFilter, sequence.getSequenceId());
					}
				}
			}
			
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}
	
	private void processSequence_KMER_COUNTING_SAMPLE_D2(SequenceD2 sequence){
		try{
			sequenceCount.incrementAndGet();
			if(sequenceCount.get() % 1000000 == 0){
				System.out.println(Utils.time()+" SequenceProcessorCMG: "+id+"\t"+mode.toString()+" "+sequenceCount.get());
			}
			
			long kmerCode;
			
			for(int i=0; i+k<=sequence.getLength(); i++){
				int oldAs = sequence.as;
				int oldTs = sequence.ts;
				int oldCs = sequence.cs;
				int oldGs = sequence.gs;
				kmerCode = Coder.encodeToLong(sequence, i, i+k, false, true);
				if(kmerCode>=0L){
					sequence.insertKmerCount(kmerCode, 1);
					sequence.insertKmerProb(kmerCode, (short)(sequence.as-oldAs), (short)(sequence.ts-oldTs), (short)(sequence.cs-oldCs), (short)(sequence.gs-oldGs));
				}
				
				//reverse
				oldAs = sequence.as;
				oldTs = sequence.ts;
				oldCs = sequence.cs;
				oldGs = sequence.gs;
				kmerCode = Coder.encodeToLong(sequence, i, i+k, true, true);
				if(kmerCode>=0L){
					sequence.insertKmerCount(kmerCode, 1);
					sequence.insertKmerProb(kmerCode, (short)(sequence.as-oldAs), (short)(sequence.ts-oldTs), (short)(sequence.cs-oldCs), (short)(sequence.gs-oldGs));
				}
			}
			
			rpSampleVector.addWith((SequenceD2Interface)sequence, false);
			sequence.clear();
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}
	
}
