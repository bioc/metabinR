/*
 *
 * MetaTarget DictionaryFilter
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

import fr.cea.ig.metatarget.utils.Utils;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;

import java.util.List;
import java.util.Set;
import java.util.concurrent.ConcurrentMap;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.MapMaker;

public class DictionaryFilter {

	private  ConcurrentMap<Long, TIntList> kmer2sequenceMap;
	private  HashMultimap<Integer, Integer> sequence2sequenceMap;
	
	public DictionaryFilter(int initialCapacity, int concurrencyLevel){
		kmer2sequenceMap = new MapMaker().concurrencyLevel(concurrencyLevel).initialCapacity(initialCapacity).makeMap();
	}
	
	/*
	public HashMultimap<Integer, Integer> getR2R(){
		return sequence2sequenceMap;
	}
	*/
	
	public List<Integer> getR2RListForId(int sequenceId){
		return Utils.asSortedList(sequence2sequenceMap.get(sequenceId));
	}
	
	public Set<Integer> getR2RSetForId(int sequenceId){
		return sequence2sequenceMap.get(sequenceId);
	}
	
	public synchronized void insert(long kmerId, final int sequenceId){
		if(kmerId<0L){
			return;
		}
		
		//Collection<Integer> old = kmer2sequenceMap.putIfAbsent(kmerId, new ConcurrentSkipListSet<Integer>(ImmutableList.of(sequenceId)));
		//Collection<Integer> old = kmer2sequenceMap.putIfAbsent(kmerId, new ConcurrentLinkedDeque<Integer>(ImmutableList.of(sequenceId)));
		TIntList old = kmer2sequenceMap.putIfAbsent(kmerId, new TIntArrayList(){{add(sequenceId);}});
		if(old != null) {
			if(!old.contains(sequenceId)){
				old.add(sequenceId);
			}
		}
		
		return;
	}
	
	public void populateSequence2SequenceRelation(int usingThsequences){
		if(sequence2sequenceMap!=null){
			sequence2sequenceMap.clear();
			sequence2sequenceMap = null;
		}
		sequence2sequenceMap = HashMultimap.create();
			
		/*
		//PARALLEL
		System.out.println("PARALLEL");
		ParallelForTrove.blockingFor(usingThsequences, kmer2sequenceMap.values(), 
	    		 // The operation to perform with each item
	    		 new ParallelForTrove.Operation<TIntList>() {
	    		    public void perform(TIntList list) {
	    		    	list.sort();
	    				for(int i=0; i<list.size(); i++){
	    					for(int j=i; j<list.size(); j++){
	    						sequence2sequenceMap.get(list.get(i)).add(list.get(j));
	    					}
	    				}
	    		    };
	    		});
		*/
		
		
		//SERIAL
		for(TIntList list : kmer2sequenceMap.values()){
			//System.out.println("sortinglist "+(++count)+" size="+list.size());
			//list.sort();
			for(int i=0; i<list.size(); i++){
				for(int j=i+1; j<list.size(); j++){
					sequence2sequenceMap.get(list.get(i)).add(list.get(j));
				}
			}
		}
	}
	
	public boolean areRelated(Sequence sequence1, Sequence sequence2){
		return sequence2sequenceMap.get(sequence1.getSequenceId()).contains(sequence2.getSequenceId()) || sequence2sequenceMap.get(sequence2.getSequenceId()).contains(sequence1.getSequenceId());
	}
	
	public boolean areRelated(int sequence1Id, int sequence2Id){
		return sequence2sequenceMap.get(sequence1Id).contains(sequence2Id) || sequence2sequenceMap.get(sequence2Id).contains(sequence1Id);
	}
	
	public String toString(int kFilter) {
		StringBuilder sb = new StringBuilder();
		int numOflists = kmer2sequenceMap.values().size();
		sb.append("numOflists:\t\t"+numOflists+"\n");
		
		return sb.toString();
	}
	
}
