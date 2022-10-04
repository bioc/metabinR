/*
 *
 * MetaTarget EMasync
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

import fr.cea.ig.metatarget.utils.Utils;
import gnu.trove.iterator.TIntLongIterator;
import gnu.trove.map.hash.TIntLongHashMap;

public class EMasync {
	private static final int maxRuns = 25;
	private static final double CONVERGE_THRESHOLD = 0.00001;
	
	private final int numberOfClusters;
	private final TIntLongHashMap countsHisto;
	private final int excludeMin;
	private final int excludeMax;
	private int kmerHighestCount;
	
	@SuppressWarnings("unused")
	private double kmerThreshold = 0.5;
	@SuppressWarnings("unused")
	private double readThreshold = 0.5;

	private double[] clusterSizes;
	private double[] clusterAbundances;
	private double[] currAbundances;
	private double[] currClassSize;
	private double[] prevClassSize;
	private double[] prevClassAbundance;	
	
	private double[][] countPoissonClass;
	//private boolean[] usability;
	
	public double[] getClusterAbundances() {
		return clusterAbundances;
	}
	
	public double[] getClusterSizes() {
		return clusterSizes;
	}

	public EMasync(int numberOfClusters, int excludeMin, int excludeMax, Dictionary dictionary){
		this.numberOfClusters = numberOfClusters;
		this.countsHisto = dictionary.getCountsHisto();
		this.excludeMin = excludeMin;
		this.excludeMax = excludeMax;
		this.kmerHighestCount = dictionary.getMaxCount();
		if(excludeMax!=0 && kmerHighestCount>excludeMax){
			kmerHighestCount = excludeMax;
		}
		
		int i, j;
		kmerThreshold = 1.0 / (double)numberOfClusters;

		clusterSizes = new double[numberOfClusters];
		clusterAbundances = new double[numberOfClusters];
		currAbundances = new double[numberOfClusters];
		currClassSize = new double[numberOfClusters];
		prevClassSize = new double[numberOfClusters];
		prevClassAbundance = new double[numberOfClusters];
		
		j = 10 * (numberOfClusters - 1);
		for (i = 0; i < numberOfClusters; i++){
			// Set abundance value with exact values
			clusterAbundances[i] = (double)j;
			//clusterAbundances[i] = (double)1;
			prevClassAbundance[i] = clusterAbundances[i];
			j = j - 10;
			if (j == 0){
				j = 1;
			}

			clusterSizes[i] = (double)1000000;
			//clusterSizes[i] = (double)1;
			prevClassSize[i] = clusterSizes[i];
		}
		
		// Initialize the probability array of every counts
		countPoissonClass = new double[kmerHighestCount + 1][numberOfClusters];
		//usability = new boolean[kmerHighestCount + 1];
	}
	
	private double safeDivideWith(double div){
		return (div==0.0 ? 1.0:div);
	}
	
	public void performEM(){	
		try {
			CountDownLatch doneSignal = new CountDownLatch(numberOfClusters);
			for (int i = 0; i < numberOfClusters; i++){
				Thread t = new Thread(new EMThread(i, doneSignal));
				t.start();
			}
			doneSignal.await();
		} 
		catch (InterruptedException e) {
			e.printStackTrace();
		}
	}
	
	private class EMThread implements Runnable {
		
		private int clusterId;
		private CountDownLatch doneSignal;
		
		private EMThread(int clusterId, CountDownLatch doneSignal){
			this.clusterId = clusterId;
			this.doneSignal = doneSignal;
		}
		
		@Override
		public void run() {
			int j, k, cnt;
			int zeroflag, convflag;
			double fTemp;
			k = 0;
			while (true){
				k++;
				System.out.println("Cluster="+clusterId+"\tRun="+k);
				
				currAbundances[clusterId] = clusterAbundances[clusterId];
				currClassSize[clusterId] = clusterSizes[clusterId];
				// Expectation
				// For each count except count = 0, calculate the probability that the given profile belongs to each class
				for (cnt = 1; cnt <= kmerHighestCount; cnt++){
					if(cnt >= excludeMin && (excludeMax==0 || cnt <= excludeMax)){
						countPoissonClass[cnt][clusterId] = 0.0;
						for (j = 0; j < numberOfClusters; j++){
							double lnSum = ClusterPoisson.lnPoissonProbabilitySum(currAbundances[j], currAbundances[clusterId], cnt, currClassSize[j]);
							double tempSum = Math.exp(lnSum);
							if(Double.isInfinite(tempSum)){
								countPoissonClass[cnt][clusterId] = Double.POSITIVE_INFINITY;
								break;
							}
							
							/*
							double temp = currClassSize[j] 
									 * Math.pow(currAbundances[j] / safeDivideWith(currAbundances[clusterId]), (double)cnt) 
									 * Math.pow(Math.E, currAbundances[clusterId] - currAbundances[j]);
							*/
							
							countPoissonClass[cnt][clusterId] += tempSum;
						}
						countPoissonClass[cnt][clusterId] = currClassSize[clusterId] / countPoissonClass[cnt][clusterId];
					}
				}
				
				
				// Maximization
				// classSize and classAbundance
				prevClassSize[clusterId] = clusterSizes[clusterId];
				clusterSizes[clusterId] = 0.0;
				prevClassAbundance[clusterId] = clusterAbundances[clusterId];
				clusterAbundances[clusterId] = 0.0;
				for ( TIntLongIterator it = countsHisto.iterator(); it.hasNext(); ) {
				    it.advance();
				    int count = it.key();
				    long countAbundance = it.value();
					if(count >= excludeMin && (excludeMax==0 || count <= excludeMax)){
						clusterSizes[clusterId] += (double)countAbundance * countPoissonClass[count][clusterId];
						clusterAbundances[clusterId] += (double)countAbundance * (countPoissonClass[count][clusterId] * (double)count);
					}
				}
				clusterAbundances[clusterId] /= safeDivideWith(clusterSizes[clusterId]);
				
				
				// Check for convergence
				zeroflag = 0;
				convflag = 0;
				if (clusterAbundances[clusterId] == 0.0 || clusterSizes[clusterId] == 0.0){
					zeroflag = 1;
				}
				else{
					fTemp = (clusterAbundances[clusterId] - prevClassAbundance[clusterId]) / safeDivideWith(prevClassAbundance[clusterId]);
					if (fTemp < 0.0){
						fTemp = fTemp * -1.0;
					}
					if (fTemp < CONVERGE_THRESHOLD){
						convflag++;
					}
					fTemp = (clusterSizes[clusterId] - prevClassSize[clusterId]) / safeDivideWith(prevClassSize[clusterId]);
					if (fTemp < 0.0){
						fTemp = fTemp * -1.0;
					}
					if (fTemp < CONVERGE_THRESHOLD){
						convflag++;
					}
				}

				if (zeroflag == 1 || convflag == 2){
					break;
				}
				
				if (k >= maxRuns){
					break;
				}
			}
			
			System.out.println(Utils.time()+"\tCluster "+clusterId+" abundance="+clusterAbundances[clusterId]);
			System.out.println(Utils.time()+"\tCluster "+clusterId+" length="+clusterSizes[clusterId]);
			System.out.println(Utils.time()+"\tCluster "+clusterId+" Runs="+k);
			
			doneSignal.countDown();
		}

	}
	
}
