package fr.cea.ig.metatarget.datastructures;

import java.util.concurrent.CountDownLatch;

import fr.cea.ig.metatarget.utils.Utils;
import gnu.trove.iterator.TIntLongIterator;
import gnu.trove.map.hash.TIntLongHashMap;

public class EMsync {
	private static final int maxRuns = 25;
	private static final double CONVERGE_THRESHOLD = 0.00001;
	
	private int numberOfClusters;
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

	public EMsync(int numberOfClusters, int excludeMin, int excludeMax, Dictionary dictionary){
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

			//clusterSizes[i] = (double)100000;
			clusterSizes[i] = (double)1;
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
			CountDownLatch doneSignal = new CountDownLatch(1);
			Thread t = new Thread(new EMThread(doneSignal));
			t.start();
			doneSignal.await();
		} 
		catch (InterruptedException e) {
			e.printStackTrace();
		}
	}
	
	private class EMThread implements Runnable {
		
		private CountDownLatch doneSignal;
		
		private EMThread(CountDownLatch doneSignal){
			this.doneSignal = doneSignal;
		}
		
		@Override
		public void run() {
			int i, j, k, cnt;
			int zeroflag, convflag;
			double fTemp;
			k = 0;
			while (true){
				k++;
				System.out.println("\tRun="+k);
				
				for (i = 0; i < numberOfClusters; i++){
					currAbundances[i] = clusterAbundances[i];
					currClassSize[i] = clusterSizes[i];
				}
				// Expectation
				// For each count except count = 0, calculate the probability that the given profile belongs to each class
				for (cnt = 1; cnt <= kmerHighestCount; cnt++){
					if(cnt >= excludeMin && (excludeMax==0 || cnt <= excludeMax)){
						for (i = 0; i < numberOfClusters; i++){
							countPoissonClass[cnt][i] = 0.0;
							for (j = 0; j < numberOfClusters; j++){
								double lnSum = ClusterPoisson.lnPoissonProbabilitySum(currAbundances[j], currAbundances[i], cnt, currClassSize[j]);
								double tempSum = Math.exp(lnSum);
								if(Double.isInfinite(tempSum)){
									countPoissonClass[cnt][i] = Double.POSITIVE_INFINITY;
									break;
								}
								countPoissonClass[cnt][i] += tempSum;
							}
							countPoissonClass[cnt][i] = currClassSize[i] / countPoissonClass[cnt][i];
						}
					}
				}
				
				
				// Maximization
				// classSize and classAbundance
				for (i = 0; i < numberOfClusters; i++){
					prevClassSize[i] = clusterSizes[i];
					clusterSizes[i] = 0.0;
					prevClassAbundance[i] = clusterAbundances[i];
					clusterAbundances[i] = 0.0;
				}
				for ( TIntLongIterator it = countsHisto.iterator(); it.hasNext(); ) {
				    it.advance();
				    int count = it.key();
				    long countAbundance = it.value();
					if(count >= excludeMin && (excludeMax==0 || count <= excludeMax)){
						for (i = 0; i < numberOfClusters; i++){
							clusterSizes[i] += (double)countAbundance * countPoissonClass[count][i];
							clusterAbundances[i] += (double)countAbundance * (countPoissonClass[count][i] * (double)count);
						}
					}
				}
				for (i = 0; i < numberOfClusters; i++){
					clusterAbundances[i] /= safeDivideWith(clusterSizes[i]);	
				}
				
				
				// Check for convergence
				zeroflag = 0;
				convflag = 0;
				for (i = 0; i < numberOfClusters; i++){
					if (clusterAbundances[i] == 0.0 || clusterSizes[i] == 0.0){
						zeroflag = 1;
					}
					else{
						fTemp = (clusterAbundances[i] - prevClassAbundance[i]) / safeDivideWith(prevClassAbundance[i]);
						if (fTemp < 0.0){
							fTemp = fTemp * -1.0;
						}
						if (fTemp < CONVERGE_THRESHOLD){
							convflag++;
						}
						fTemp = (clusterSizes[i] - prevClassSize[i]) / safeDivideWith(prevClassSize[i]);
						if (fTemp < 0.0){
							fTemp = fTemp * -1.0;
						}
						if (fTemp < CONVERGE_THRESHOLD){
							convflag++;
						}
					}
				}

				if (zeroflag == 1 || convflag == 2*numberOfClusters){
					break;
				}
				
				if (k >= maxRuns){
					break;
				}
			}
			
			System.out.println(Utils.time()+"\tRuns="+k);
			
			doneSignal.countDown();
		}

	}
	
}
