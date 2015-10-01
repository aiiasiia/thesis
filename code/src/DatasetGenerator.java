/** 
  * DatasetGenerator.java
  * 
  * Randomly plants (l,d) motif once in each of
  *  t (def 20) randomly generated DNA sequences, each
  *  n (def 600) bases long.
  * 
  * @author Juliet Q. Nabos
  * @version 2.1 3/02/2014
 */

import java.util.*;
import java.io.*;
//import java.util.Scanner;
// sent as MF_DataSetGenerator_v2_1.java
public class DatasetGenerator {
	static int t = 20; //number of sequences
	static int n = 600; //length of sequence
	static int l = 15; //length of motif
	static int d = 4; //number of mutations
	static int set;
	static Random rn = new Random();
	static char[][] seqS;
	static char[] motif;
	static int[] mutationSite;
	static char[][] motifAlignment;
	static int[] alignment;
	static char base = 0;


	public static void main (String[] args) throws Exception {
		if( args.length == 4 ) {
			t = Integer.parseInt( args[0] );
			n = Integer.parseInt( args[1] );
			l = Integer.parseInt( args[2] );
			d = Integer.parseInt( args[3] );
		}
		else if( args.length == 3) {
			t = 20;
			n = 600;
			l = Integer.parseInt( args[0] );
			d = Integer.parseInt( args[1] );
			set = Integer.parseInt( args[2] );
		}
		else {
			System.out.println("Usage:");
			System.out.println("\t    java DatasetGenerator <t> <n> <l> <d>");
			System.out.println("\t or java DatasetGenerator <l> <d> <set>,\t t=20, n=600");
			return;
		}

		seqS = new char[t][n];
		motif = new char[l];
		mutationSite = new int[d];
		motifAlignment = new char[t][l];
		alignment = new int[t];

		String filename = l + "," + d + "," + set;
		PrintStream outToFile = new PrintStream( new File (filename));	
		outToFile.println(t);
		outToFile.println(n);
		outToFile.println(l);
		outToFile.println(d);
		//generate data set:
		generateSeq(t, n);	
		generateMotif(l);
		//mutate motif 
		char[] gappedMotif = new char[l];
		for (int i = 0; i < t; i++) {
			generateMutationSite(l,d);	
			//copy motif[] to gappedMotif[] for mutation
			for (int ii = 0; ii < l; ii++) {
				gappedMotif[ii] = motif[ii];
			}
			//change mutation sites
			for (int j = 0; j < d; j++) {
				double randomBase = rn.nextDouble();
				selectBase(randomBase);
				if (base == gappedMotif[mutationSite[j]]) {
					j--;
					continue;
				}
				else {
					gappedMotif[mutationSite[j]] = base;
				}
			}
			for (int j = 0; j < l; j++) {
				motifAlignment[i][j] = gappedMotif[j];
			}
		}
		
		//plant gapped motifs
		for (int i = 0; i < t; i++) {
			int k = 0;
			int plantingSite = rn.nextInt(n-l+1);
			alignment[i] = plantingSite;
			for (int j = plantingSite; j < (plantingSite + l) ; j++) {
				seqS[i][j] = motifAlignment[i][k];
				k++;
			}
		}
		
		//print alignment
		for (int i = 0; i < t; i++) {
			outToFile.print(alignment[i] + " ");
		}
		outToFile.println();
		
		//print motif
		for (int i = 0; i < l; ++i) {
			outToFile.print(motif[i]);
		}
		outToFile.println();
		
		//write seqS to file
		for (int i = 0; i < t; i++) {
			for (int j = 0; j < n; j++) {
				outToFile.print(seqS[i][j]);
			}
			outToFile.println();
		}
		outToFile.close();
		
		System.out.println("Done generating data set: " + filename);

	}//end of method main
	

	public static void generateMotif(int ll) {
		for (int i= 0; i < ll; i++) {
			double randomMotifKey = rn.nextDouble();
			selectBase(randomMotifKey);
			motif[i] = base;
		}
	}//end of method generateMotif
	
	public static void generateMutationSite(int ll, int dd) {
		//clean mutationSite[]
		for (int i = 0; i < mutationSite.length; i++) {
			mutationSite[i] = -1;
		}
		//generate mutation sites
		for (int i = 0; i < dd; i++) {
			int randomSite = rn.nextInt(ll);
			if (i == 0) {
				mutationSite[i] = randomSite;
			}
			else {
				for (int j = 0; j < i; j++) {
					if (randomSite == mutationSite[j]) {
						i--;
						break;
					}
					else {
						mutationSite[i] = randomSite;
					}
				}
			}
		}
		
	}//end of method generateMutationSite
	
	public static void selectBase(double randU) {
		if (randU >= 0.0 && randU <= 0.25) {
			base = 'A';
		}
		else if (randU > 0.25 && randU <= 0.50) {
			base = 'C';
		}
		else if (randU > 0.50 && randU <= 0.75) {
			base = 'G';
		}
		else if (randU > 0.75 && randU <= 1.00) {
			base = 'T';
		}
	}//end of method selectBase
	
	public static void generateSeq(int sequences, int LoS) {
		for (int i = 0; i < t; i++) {
			for (int j = 0; j < n; j++) {
				double randomBase = rn.nextDouble();
				selectBase(randomBase);
				seqS[i][j] = base;
			}
		}
	}//end of method generateSeq

}//end of class MF_DatasetGenerator