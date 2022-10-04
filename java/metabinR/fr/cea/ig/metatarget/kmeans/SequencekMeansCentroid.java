/*
 *
 * MetaTarget SequencekMeansCentroid
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
package fr.cea.ig.metatarget.kmeans;

import java.util.Map;
import java.util.TreeMap;
import java.util.Map.Entry;

import fr.cea.ig.metatarget.datastructures.Sequence;
import fr.cea.ig.metatarget.utils.Utils;
import gnu.trove.iterator.TLongIntIterator;
import gnu.trove.map.hash.TLongDoubleHashMap;

public class SequencekMeansCentroid {
	
	private TLongDoubleHashMap kmerValues = null;
	
	private int[] ranks = null;
	
	public SequencekMeansCentroid(Sequence read) {
		initKmerValues();
		addWith(read);
	}

	public void initRanks(Map<Long, Integer> spaceRanks){
		int spaceSize = spaceRanks.size(); 
		ranks = new int[spaceSize];
	
		//put present kmers in a hashmap
		TreeMap<Long, Double> unsorted = new TreeMap<Long, Double>();
		for(long kmerCode : getKmerCodes()){
			unsorted.put(kmerCode, getValueForKmerCode(kmerCode));
		}
		//System.out.println("unsorted1 "+unsorted.size());
		//System.out.println(unsorted);
		
		//insert zeros for non present kmers
		for(long kmerCode : spaceRanks.keySet()){
			if(!unsorted.containsKey(kmerCode)){
				unsorted.put(kmerCode, 0.0);
			}
		}
		//System.out.println("unsorted2 "+unsorted.size());
		//System.out.println(unsorted);
		
		//sort the map
		Map<Long, Double> sorted =  Utils.sortByValue(unsorted);	
		//System.out.println("sorted "+sorted.size());
		//System.out.println(sorted);
		
		//write ranks
		int current = 0;
		for(Entry<Long, Double> entry : sorted.entrySet()){
			ranks[spaceRanks.get(entry.getKey())] = current++;
		}
		//System.out.println("ranks "+ranks.length);
		//System.out.println(Arrays.toString(ranks));
		//System.out.println("space "+spaceRanks.size());
		//System.out.println(spaceRanks.toString());
		
		unsorted.clear();
		sorted.clear();
		unsorted = null;
		sorted = null;
		//System.exit(0);
	}
	
	public int[] getRanks(){
		return ranks;
	}
	
	public void initKmerValues(){
		kmerValues = new TLongDoubleHashMap();
	}
	
	//public Collection<Long> getKmerCodes(){
	public long[] getKmerCodes(){
		//return kmerValues.keySet();
		return kmerValues.keys();
	}
	
	public double getValueForKmerCode(long kmerCode){
		return kmerValues.get(kmerCode);
	}
	
	public void insertKmerValue(long kmerCode, double value){
		if(kmerCode<0L){
			return;
		}
		kmerValues.adjustOrPutValue(kmerCode, value, value);
	}
	
	public void clear(){
		if(kmerValues!=null){
			kmerValues.clear();
			kmerValues = null;
		}
	}
	
	public void addWith(Sequence other){
		for ( TLongIntIterator it = other.iteratorCounts(); it.hasNext(); ) {
			it.advance();
			insertKmerValue(it.key(), (double)it.value());
		}
	}
	
	public void divideWith(double div){
		for(long kmerCode : getKmerCodes()){
			kmerValues.put(kmerCode, kmerValues.get(kmerCode)/div);
		}
	}
	
	public static double distanceSpearman(SequencekMeansCentroid center, Sequence read, Map<Long, Integer> spaceRanks) {
    	double distance = 0.0;
    	int spaceSize = spaceRanks.size();
    	for(int i=0; i<spaceSize; i++){
    		distance += (double)Math.abs(center.getRanks()[i] - read.getRanks()[i]);
		}
        return (double)distance;// /(double)(spaceSize*spaceSize);
    }
	
    public static double distanceEuclid(SequencekMeansCentroid center, Sequence read, Map<Long, Integer> spaceRanks) {
    	double sumSquared = 0.0;
    	for(long kmerCode : spaceRanks.keySet()){
    		double v = (double)read.getCountForKmerCode(kmerCode) - center.getValueForKmerCode(kmerCode);
    		sumSquared += v*v;	
		}
        return Math.sqrt(sumSquared);
    }
	
}
