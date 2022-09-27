/*
 *
 * MetaTarget VectorUtils
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

import java.util.HashSet;

import fr.cea.ig.metatarget.utils.Utils;
import gnu.trove.iterator.TLongIntIterator;

public class VectorUtils {

	public static double cosineSequence2Cluster(Sequence sequence, ClusterVectorTrove clusterVector) {
		if(clusterVector == null || clusterVector.isEmpty()){
			return Double.NEGATIVE_INFINITY;
		}
		double similarity = 0.0;
		double clusterNorm = clusterVector.getNorm(); //this);
		double sequenceNorm = 0.0;
		for ( TLongIntIterator it = sequence.iteratorCounts(); it.hasNext(); ) {
			it.advance();
			long kmerCode = it.key();
			int sequenceKmerCount = it.value();
			//int kmerGlobalCount = dictionary.getGlobalCountFor(kmerCode);
			double sequenceKmerWeight = (double)sequenceKmerCount;
			//double sequenceKmerWeight = (double)sequenceKmerCount / Math.sqrt((double)kmerGlobalCount);
			sequenceNorm += sequenceKmerWeight * sequenceKmerWeight;
			double clusterKmerWeight = clusterVector.getCountForKmerCode(kmerCode);
			if(clusterKmerWeight != 0.0){
				similarity += clusterKmerWeight * sequenceKmerWeight;
			}
		}
		sequenceNorm = Math.sqrt(sequenceNorm);
		
		if(clusterNorm==0.0) clusterNorm=1.0;
		if(sequenceNorm==0.0) sequenceNorm=1.0;
		return similarity / (clusterNorm * sequenceNorm);
	}
	
	public static ClusterPoisson[] createABClusterPoissonsEMsync(int numOfClusters, int excludeMin, int excludeMax, Dictionary dictionary){
		ClusterPoisson[] clusterPoissons = new ClusterPoisson[numOfClusters];

		//EM
		EMsync em = new EMsync(numOfClusters, excludeMin, excludeMax, dictionary);
		System.out.println(Utils.time()+" START of EMsync");
		em.performEM();
		System.out.println(Utils.time()+" END of EMsync");
		
		double[] clusterAbundances = em.getClusterAbundances();
		double[] clusterSizes = em.getClusterSizes();
		for(int i=0; i<numOfClusters; i++){
			double abundance = clusterAbundances[i];
			double size = clusterSizes[i];
			ClusterPoisson cp = new ClusterPoisson(abundance, size);
			clusterPoissons[i] = cp;
			
			double sigma = Math.sqrt(abundance);
			int highLimit = (int)(abundance + 5.0*sigma);
			int lowLimit = (int)(abundance - 5.0*sigma);
			
			if(lowLimit>highLimit){
				highLimit = lowLimit;
			}
			
			if(lowLimit<excludeMin){
				lowLimit=excludeMin;
			}
			
			if(highLimit>excludeMax && excludeMax!=0){
				highLimit=excludeMax;
			}
			
			cp.setHighLimit(highLimit);
			cp.setLowLimit(lowLimit);
		}

		return clusterPoissons;
	}
	
	public static ClusterPoisson[] createABClusterPoissonsEMasync(int numOfClusters, int excludeMin, int excludeMax, Dictionary dictionary){
		ClusterPoisson[] clusterPoissons = new ClusterPoisson[numOfClusters];

		//EM
		EMasync em = new EMasync(numOfClusters, excludeMin, excludeMax, dictionary);
		System.out.println(Utils.time()+" START of EMasync");
		em.performEM();
		System.out.println(Utils.time()+" END of EMasync");
		
		double[] clusterAbundances = em.getClusterAbundances();
		double[] clusterSizes = em.getClusterSizes();
		for(int i=0; i<numOfClusters; i++){
			double abundance = clusterAbundances[i];
			double size = clusterSizes[i];
			ClusterPoisson cp = new ClusterPoisson(abundance, size);
			clusterPoissons[i] = cp;
			
			double sigma = Math.sqrt(abundance);
			int highLimit = (int)(abundance + 3.0*sigma);
			int lowLimit = (int)(abundance - 3.0*sigma);
			
			if(lowLimit>highLimit){
				highLimit = lowLimit;
			}
			
			if(lowLimit<excludeMin){
				lowLimit=excludeMin;
			}
			
			if(highLimit>excludeMax && excludeMax!=0){
				highLimit=excludeMax;
			}
			
			cp.setHighLimit(highLimit);
			cp.setLowLimit(lowLimit);
		}

		return clusterPoissons;
	}
	
	public static ClusterPoisson[] filterClusterPoissons(ClusterPoisson[] oldClusterPoissons) {
		HashSet<ClusterPoisson> newClusterPoissons = new HashSet<ClusterPoisson>();
		
		for(ClusterPoisson oldCp : oldClusterPoissons){
			if(!newClusterPoissons.isEmpty()){
				boolean mrgd = false;
				for(ClusterPoisson newCp : newClusterPoissons){
					if(Math.abs((newCp.getGenomeAbundance()-oldCp.getGenomeAbundance())/newCp.getGenomeAbundance())<0.01){
						newClusterPoissons.remove(newCp);
						ClusterPoisson merged = merge(newCp,oldCp); 
						newClusterPoissons.add(merged);
						mrgd = true;
						break;
					}
				}
				if(!mrgd){
					newClusterPoissons.add(oldCp);
				}
			}
			else{
				newClusterPoissons.add(oldCp);
			}
		}
		
		ClusterPoisson[] ret = new ClusterPoisson[newClusterPoissons.size()];
		int i=0;
		for(ClusterPoisson cp : newClusterPoissons){
			ret[i++] = cp;
		}
		return ret;
	}
	
	private static ClusterPoisson merge(ClusterPoisson a, ClusterPoisson b){
		return new ClusterPoisson(a.getGenomeAbundance()>b.getGenomeAbundance()?a.getGenomeAbundance():b.getGenomeAbundance(), a.getGenomeLength()>b.getGenomeLength()?a.getGenomeLength():b.getGenomeLength());
	}
}
