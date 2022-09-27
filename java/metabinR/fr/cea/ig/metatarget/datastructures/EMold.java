package fr.cea.ig.metatarget.datastructures;

import java.util.Map.Entry;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.atomic.AtomicInteger;

import fr.cea.ig.metatarget.utils.Utils;

public class EMold {
	private static final int maxRuns = 25;
	private static final double CONVERGE_THRESHOLD = 0.00001;
	
	private int numberOfClusters;
	private Dictionary dictionary;
	private int excludeMin;
	private int excludeMax;
	
	@SuppressWarnings("unused")
	private double kmerThreshold = 0.5;
	@SuppressWarnings("unused")
	private double readThreshold = 0.5;
	private int kmerHighestCount;

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

	public EMold(int numberOfClusters, int excludeMin, int excludeMax, Dictionary dictionary){
		this.numberOfClusters = numberOfClusters;
		this.dictionary = dictionary;
		this.excludeMin = excludeMin;
		this.excludeMax = excludeMax;
		kmerHighestCount = this.dictionary.getMaxCount();
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
	
	/*
	private void checkNaN(double value, String msg){
		if(Double.isNaN(value)){
			System.out.println("NaN "+msg);
			System.exit(0);
		}
	}
	*/
	
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
							//checkNaN(clusterSizes[j], "A");
							//checkNaN(Math.pow(clusterAbundances[j] / safeDivideWith(clusterAbundances[clusterId]), (double)cnt), "B");
							//checkNaN(Math.pow(Math.E, clusterAbundances[clusterId] - clusterAbundances[j]), "C");
							
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
							//checkNaN(countPoissonClass[cnt][clusterId], "00");
						}
						countPoissonClass[cnt][clusterId] = currClassSize[clusterId] / countPoissonClass[cnt][clusterId];
						//checkNaN(countPoissonClass[cnt][clusterId], "0");
					}
				}
				
				// Maximization
				// classSize and classAbundance
				//checkNaN(clusterAbundances[clusterId], "1");
				prevClassSize[clusterId] = clusterSizes[clusterId];
				clusterSizes[clusterId] = 0.0;
				prevClassAbundance[clusterId] = clusterAbundances[clusterId];
				clusterAbundances[clusterId] = 0.0;
				for(Entry<Long, AtomicInteger> entry : dictionary.getHashMap().entrySet()){
					int count = entry.getValue().get();
					if(count >= excludeMin && (excludeMax==0 || count <= excludeMax)){
						clusterSizes[clusterId] += countPoissonClass[count][clusterId];
						clusterAbundances[clusterId] += (countPoissonClass[count][clusterId] * (double)count);
						//checkNaN(clusterAbundances[clusterId], "2");
					}
				}
				//checkNaN(clusterAbundances[clusterId], "3");
				clusterAbundances[clusterId] /= safeDivideWith(clusterSizes[clusterId]);
				//checkNaN(clusterAbundances[clusterId], "4");
				
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
