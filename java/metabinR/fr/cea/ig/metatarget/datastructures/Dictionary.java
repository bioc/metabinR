/*
 *
 * MetaTarget Dictionary
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
import gnu.trove.map.hash.TIntDoubleHashMap;
import gnu.trove.map.hash.TIntLongHashMap;

import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.atomic.AtomicInteger;

//import com.google.common.collect.MapMaker;

public class Dictionary implements DictionaryInterface{
	@SuppressWarnings("unused")
	private static int _excludeMin = 1;
	@SuppressWarnings("unused")
	private static int _excludeMax = 0;
	
	private ConcurrentMap<Long, AtomicInteger> kmerCountMap;
	
	private TIntLongHashMap countsHisto;
	
	private long minKmerCode;
	private long maxKmerCode;
	
	private long uniqueKmers = 0;
	private long distinctKmers = 0;
	private long totalKmers = 0;
	private int maxKmerCount = 0;
	
	public Dictionary(int initialCapacity, int concurrencyLevel, int excludeMin, int excludeMax){
		_excludeMin = excludeMin;
		_excludeMax = excludeMax;
		
		System.out.println(Utils.time()+" Using Dictionary ORIGINAL");
		kmerCountMap = new ConcurrentHashMap<Long, AtomicInteger>(initialCapacity, 0.75f, concurrencyLevel);
		//kmerCountMap = new MapMaker().concurrencyLevel(concurrencyLevel).initialCapacity(initialCapacity).makeMap();
		//kmerCountMap = Utils.getDB().getTreeMap("dictionary");
	
		minKmerCode = Long.MAX_VALUE;
		maxKmerCode = Long.MIN_VALUE;
	}
	
	public long getUniqueKmers() {
		return uniqueKmers;
	}

	public long getDistinctKmers() {
		return distinctKmers;
	}

	public long getTotalKmers() {
		return totalKmers;
	}

	public void clear(){
		if(countsHisto!=null){
			countsHisto.clear();
		}
		countsHisto = null;
		
		if(kmerCountMap!=null){
			kmerCountMap.clear();
		}
		kmerCountMap = null;
	}
	
	public Map<Long, AtomicInteger> getHashMap(){
		return kmerCountMap;
	}
	
	public Set<Long> getKmerCodes(){
		return kmerCountMap.keySet();
	}
	
	/*
	public synchronized Map<Long, Integer> getKmerRanks(){
		if(kmerRanks==null){
			kmerRanks = new TreeMap<Long, Integer>();
			int rank = 0;
			for(long kmer : getKmerCodes()){
				kmerRanks.put(kmer, rank++);
			}
		}
		return kmerRanks;
	}
	*/
	
	public void insert(long kmerCode){
		if(kmerCode<0L){
			return;
		}
		
		AtomicInteger oldCount = kmerCountMap.putIfAbsent(kmerCode, new AtomicInteger(1));
		if(oldCount != null) {
			oldCount.incrementAndGet();
		}
	}
	
	public void removeAll(int minCount){
		for(long kmer : getKmerCodes()){
			if(getGlobalCountFor(kmer) < minCount){
				kmerCountMap.remove(kmer);
			}
		}
	}
	
	public int getMaxCount() {
		if(maxKmerCount == 0){
			createCountsHisto();
		}
		return maxKmerCount;
	}

	public TIntLongHashMap getCountsHisto() {
		if(countsHisto == null){
			createCountsHisto();
		}
		return countsHisto;
	}

	public int getGlobalCountFor(long kmer){
		AtomicInteger gc = kmerCountMap.get(kmer);
		if(gc != null){
			return gc.get();
		}
		else{
			return 0;
		}
	}
	
	/*
	public void printKeys(){
		for(Entry<Long, AtomicInteger> entry : kmerCountMap.entrySet()){
			System.out.println(entry.getKey()+"\t"+entry.getValue().intValue());
		}
	}
	*/
	
	private TIntLongHashMap createCountsHisto(){
		StringBuilder sb = new StringBuilder();

		countsHisto = new TIntLongHashMap();
		uniqueKmers = 0L;
		distinctKmers = 0L;
		totalKmers = 0L;
		maxKmerCount = 0;
		for(Entry<Long, AtomicInteger> entry : kmerCountMap.entrySet()){
			distinctKmers++;
			int count = entry.getValue().get();
			countsHisto.adjustOrPutValue(count, 1L, 1L);
			
			if(count == 1){
				uniqueKmers++;
			}
			totalKmers += (long)count;
			if(count > maxKmerCount){
				maxKmerCount = count;
			}
			
			long kmerCode = entry.getKey();
			if(kmerCode<minKmerCode){
				minKmerCode = kmerCode;
			}
			if(kmerCode>maxKmerCode){
				maxKmerCode = kmerCode;
			}
		}
		sb.append("Unique:\t\t"+uniqueKmers+"\n");
		sb.append("Distinct:\t"+distinctKmers+"\n");
		sb.append("Total:\t\t"+totalKmers+"\n");
		sb.append("MaxCount:\t"+maxKmerCount+"\n");
		System.out.println(Utils.time());
		System.out.println(sb.toString());
		
		return countsHisto;
	}
	
	public ClusterVectorAB[] createABClusterVectors(ClusterPoisson[] clusterPoissons, int excludeMin, int excludeMax){
		int numberOfClusters = clusterPoissons.length;
		System.out.println(Utils.time()+" START of Creating AB Cluster Vectors");
		ClusterVectorAB[] clusterVectors = new ClusterVectorAB[numberOfClusters];
		try {
			CountDownLatch doneSignal = new CountDownLatch(numberOfClusters);
			for (int i = 0; i < numberOfClusters; i++){
				System.out.println("\tCluster "+(i+1)+"\tAbundance="+clusterPoissons[i].getGenomeAbundance()+"\tLength="+clusterPoissons[i].getGenomeLength()+"\tLowLimit="+clusterPoissons[i].getLowLimit()+"\tHighLimit="+clusterPoissons[i].getHighLimit());
				Thread t = new Thread(new CreateABClusterVectorThread(i, clusterPoissons, clusterVectors, doneSignal, excludeMin, excludeMax, this));
				t.start();
			}
			doneSignal.await();
		} 
		catch (InterruptedException e) {
			e.printStackTrace();
		}
		System.out.println(Utils.time()+" END of Creating AB Cluster Vectors");
		
		return clusterVectors;
	}
	
	private class CreateABClusterVectorThread implements Runnable {

		private int clusterId;
		private ClusterPoisson[] clusterPoissons;
		private ClusterVectorAB[] clusterVectors;
		private CountDownLatch doneSignal;
		private int excludeMin;
		private int excludeMax;
		@SuppressWarnings("unused")
		private Dictionary dictionary;
		
		private CreateABClusterVectorThread(int clusterId, ClusterPoisson[] clusterPoissons, ClusterVectorAB[] clusterVectors, CountDownLatch doneSignal, int excludeMin, int excludeMax, Dictionary dictionary){
			this.clusterId = clusterId;
			this.clusterPoissons = clusterPoissons;
			this.clusterVectors = clusterVectors;
			this.doneSignal = doneSignal;
			this.excludeMin = excludeMin;
			this.excludeMax = excludeMax;
			this.dictionary = dictionary;
		}
		
		@Override
		public void run() {
			//calculate probabilities
			TIntDoubleHashMap probabilities = new TIntDoubleHashMap();
			TIntLongHashMap countsHisto = getCountsHisto();
			for(int count : countsHisto.keys()){
				probabilities.put(count, clusterPoissons[clusterId].getProbability(count, clusterPoissons, clusterId));
			}
			
			long expectedSize = distinctKmers;
			if(expectedSize<=0) expectedSize = Integer.MAX_VALUE;
			//System.out.println(Utils.time()+"\tExcpected Vector Size = "+expectedSize);
			//System.out.println(Utils.time()+"\tsplitKmerCode = "+(maxKmerCode+minKmerCode)/2L);
			
			ClusterVectorAB cv = new ClusterVectorAB(expectedSize, (maxKmerCode+minKmerCode)/2L);
			//ClusterVectorAbstract cv = new ClusterVectorKoloboke(expectedSize, (maxKmerCode+minKmerCode)/2L);
			//ClusterVectorAbstract cv = new ClusterVectorHPPC(expectedSize);
			//ClusterVectorAbstract cv = new ClusterVectorGS(expectedSize);
			
			clusterVectors[clusterId] = cv;
			
			int count;
			//double probability;
			//int highLimit = clusterPoissons[clusterId].getHighLimit();
			//int lowLimit = clusterPoissons[clusterId].getLowLimit();
			for(Entry<Long, AtomicInteger> e : kmerCountMap.entrySet()){
				count = e.getValue().get();
				if(count >= excludeMin && (excludeMax==0 || count <= excludeMax)){
					//if(count>=lowLimit && count<=highLimit){
						cv.insertKMer(e.getKey(), probabilities.get(count));
					//}
				}
			}
			
			double norm = cv.getNorm(); //dictionary);
			System.out.println(Utils.time()+"\tCluster "+(clusterId+1)+"\tsize in kmers="+cv.size()+"\tnorm="+norm);
			doneSignal.countDown();
		}
		
	}
	
	
	/*
	public String toString(boolean verbose, int k) {
		long unique = 0;
		long distinct = 0;
		long total = 0;
		long max_count = 0;
		
		StringBuilder sb = new StringBuilder();
		if(verbose){
			sb.append("DictionaryCMG [\n");
		}
		//LongMapIterator<AtomicInteger> iter = hashmap.longMapIterator();
		//while(iter.moveToNext()){
		for(Entry<Long, AtomicInteger> e : kmerCountMap.entrySet()){
			distinct++;
			AtomicInteger count = e.getValue();
			//AtomicInteger count = iter.value();
			if(count.get() == 1)
				unique++;
			total += (long)count.get();
			if((long)count.get() > max_count)
				max_count = (long)count.get();
			if(verbose){
				//sb.append(e.getKey().toString(k, coder)+"=[count="+count+"],\n");
			}
		}
		if(verbose)
			sb.append("]\n\n");
		sb.append("Unique:\t\t"+unique+"\n");
		sb.append("Distinct:\t"+distinct+"\n");
		sb.append("Total:\t\t"+total+"\n");
		sb.append("Max_count:\t"+max_count+"\n");

		return sb.toString();
	}
	*/
	
}
