/*
 *
 * MetaTarget Measures
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
package fr.cea.ig.metatarget.utils;

public class Measures {
	
	// class       index c [1,numOfClasses]
	// bin/cluster index k [1,numOfBins]
	// class       index 0 not used
	// bin/cluster index 0 is the null bin/cluster
	// N is totalNumOfReads
	// N_k is reads in cluster k
	// N_c is reads in class   c that are clustered
	public static double ARI(int numOfClasses, int numOfBins, long totalNumOfReads, long[] trueNumReadsPerClass, long[][] countsBinPerClass){
		double t1 = 0.0;
		double t2 = 0.0;
		double t3 = 0.0;
		double t4 = 0.0;
		for(int c=1; c<=numOfClasses; c++){
			t1 += ((double)trueNumReadsPerClass[c] * ((double)trueNumReadsPerClass[c]) - 1.0) / 2.0;
		}
		for(int k=1; k<=numOfBins; k++){
			int N_k = 0;
			for(int c=1; c<=numOfClasses; c++){
				N_k += countsBinPerClass[k][c];
			}
			t2 += ((double)N_k * ((double)N_k) - 1.0) / 2.0;
		}
		t3 = (2.0 * t1 * t2) / ((double)totalNumOfReads * ((double)totalNumOfReads - 1.0));
		for(int c=1; c<=numOfClasses; c++){
			for(int k=1; k<=numOfBins; k++){
				t4 += ((double)countsBinPerClass[k][c] * ((double)countsBinPerClass[k][c]) - 1.0) / 2.0;
			}
		}
		
		double ARI = (t4-t3) / (0.5*(t1+t2)-t3);
		return ARI;
	}
	
	// class       index c [1,numOfClasses]
	// bin/cluster index k [1,numOfBins]
	// class       index c=0 not used
	// bin/cluster index k=0 is the null bin/cluster
	// N is totalNumOfReads
	// N_k is reads in cluster k
	// N_c is reads in class   c that are clustered
	public static double Homogeneity(int numOfClasses, int numOfBins, long totalNumOfReads, long[] trueNumReadsPerClass, long[][] countsBinPerClassPre){
		// homogeneity
		// H(Classes|Clusters)(C|K) conditional entropy of the classes given the cluster assignments
		
		int totalNumOfReadsClustered = 0;
		
		//for the null cluster
		long[][] countsBinPerClass = new long[numOfBins+1][numOfClasses+1];
		for(int k=1; k<=numOfBins; k++){
			for(int c=1; c<=numOfClasses; c++){
				countsBinPerClass[k][c]  = countsBinPerClassPre[k][c];
				countsBinPerClass[0][c] += countsBinPerClassPre[k][c];
				totalNumOfReadsClustered += countsBinPerClassPre[k][c];
			}
		}
		for(int c=1; c<=numOfClasses; c++){
			countsBinPerClass[0][c] = (trueNumReadsPerClass[c] - countsBinPerClass[0][c]);
			//System.out.println(countsBinPerClass[0][c]);
		}
		
		//H(C|K)
		double homogeneity = 0.0;
		double Hck = 0.0;
		for(int k=1; k<=numOfBins; k++){
			int N_k = 0;
			for(int c=1; c<=numOfClasses; c++){
				N_k += countsBinPerClass[k][c];
			}
			for(int c=1; c<=numOfClasses; c++){
				if(countsBinPerClass[k][c] != 0 && N_k != 0){
					//H(C|K) = Sum_K(Sum_C(a_ck/N * log(a_ck/N_k)))
					Hck -= ((double)countsBinPerClass[k][c] / (double)totalNumOfReadsClustered) * Math.log((double)countsBinPerClass[k][c] / (double)N_k);
				}
			}
		}
		//Hck = -1.0 * Hck;
		
		//H(C)
		double Hc = 0.0;
		for(int c=1; c<=numOfClasses; c++){
			int N_c = 0;
			for(int k=1; k<=numOfBins; k++){
				N_c += countsBinPerClass[k][c];
			}
			if(N_c != 0){
				//H(C) = Sum_C(N_c/N * log(N_c/N))
				Hc -= ((double)N_c / (double)totalNumOfReadsClustered) * Math.log((double)N_c / (double)totalNumOfReadsClustered);
			}
		}
		//Hc = -1.0 * Hc;
		
		//System.out.println("H "+Hck+" "+Hc);
		if(Hck != 0.0){
			homogeneity = 1.0 - (Hck / Hc); 
		}
		else{
			homogeneity = 1.0;
		}
		
		return homogeneity;
	}
	
	
	// class       index c [1,numOfClasses]
	// bin/cluster index k [1,numOfBins]
	// class       index 0 not used
	// bin/cluster index 0 is the null bin/cluster
	// N is totalNumOfReads
	// N_k is reads in cluster k
	// N_c is reads in class   c that are clustered
	public static double Completeness(int numOfClasses, int numOfBins, long totalNumOfReads, long[] trueNumReadsPerClass, long[][] countsBinPerClassPre){
		// completeness
		// H(Clusters|Classes)(K|C) conditional entropy of the cluster assignments given the classes 
		
		@SuppressWarnings("unused")
		int totalNumOfReadsClustered = 0;
		
		//for the null cluster
		long[][] countsBinPerClass = new long[numOfBins+1][numOfClasses+1];
		for(int k=1; k<=numOfBins; k++){
			for(int c=1; c<=numOfClasses; c++){
				countsBinPerClass[k][c]  = countsBinPerClassPre[k][c];
				countsBinPerClass[0][c] += countsBinPerClassPre[k][c];
				totalNumOfReadsClustered += countsBinPerClassPre[k][c];
			}
		}
		for(int c=1; c<=numOfClasses; c++){
			countsBinPerClass[0][c] = (trueNumReadsPerClass[c] - countsBinPerClass[0][c]);
			//System.out.println(countsBinPerClass[0][c]);
		}
		
		//H(K|C)
		double completeness = 0.0;
		double Hkc = 0.0;
		for(int c=1; c<=numOfClasses; c++){
			int N_c = 0;
			for(int k=0; k<=numOfBins; k++){
				N_c += countsBinPerClass[k][c];
			}
			for(int k=1; k<=numOfBins; k++){
				if(countsBinPerClass[k][c] != 0 && N_c != 0){
					//H(K|C) = Sum_C(Sum_K(a_ck/N * log(a_ck/N_c)))
					Hkc -= ((double)countsBinPerClass[k][c] / (double)totalNumOfReads) * Math.log((double)countsBinPerClass[k][c] / (double)N_c);
				}
			}
		}
		//Hkc = -1.0 * Hkc;
		
		//H(K)
		double Hk = 0.0;
		for(int k=1; k<=numOfBins; k++){
			int N_k = 0;
			for(int c=1; c<=numOfClasses; c++){
				N_k += countsBinPerClass[k][c];
			}
			if(N_k != 0){
				//H(K) = Sum_K(N_k/N * log(N_k/N))
				Hk -= ((double)N_k / (double)totalNumOfReads) * Math.log((double)N_k / (double)totalNumOfReads);
			}
		}
		//Hk = -1.0 * Hk;
		
		//System.out.println("C "+Hkc+" "+Hk);
		if(Hkc != 0.0){
			completeness = 1.0 - (Hkc / Hk); 
		}
		else{
			completeness = 1.0;
		}
		
		return completeness;
	}

}
