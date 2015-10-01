/**
  * EMS-GT.java
  * Solves the (l,d) planted motif problem.
  * Bit-array width is 32 bits (Integer type).
  *
  * @author Aia Sia, Julieta Nabos
  * @version 1.0 9/09/2015
  */

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class EMS_GT {

	static String inputFileHeader = "../datasets/";
	static String inputFileName = "";
	static final int MAX_INDICES = 120000000;
	static final int tPrime = 11;

	// from input file
	static int t, l, n, d;
	static int[] plantedAlignment;
	static String plantedMotif;
	static String foundMotifs;
	static char[][] seqS;

	static long mask;
	static int nl1;
	static int pmt;
	static int[] currNeighborhood;
	static int[] candidateMotifs;
	static long[][] lmerMappings;

//========================================================================================================
// MAIN METHOD
//========================================================================================================

	public static void main(String args[]) throws Exception {
		if( args.length > 0 ) {
			inputFileName = args[0];
		}

		readInput( inputFileHeader + inputFileName );

		// compute preliminaries
		mask = (((long) 1) << (2*l-2)) - 1;
		nl1 = n - l + 1;
		pmt = (int) (( (long)1 << (2*l) ) >>> 5); // (4^l)/32
		
		System.out.print("\n" + inputFileName);
		//System.out.println(plantedMotif);

		// start EMS-GT
		long sTime = System.nanoTime();
		collectCandidateMotifs();
		transformLmerSequences(tPrime);
		searchMotif();		
		long eTime = System.nanoTime();

		// compute runtime and used memory
		Runtime runtime = Runtime.getRuntime();
		long memUse = runtime.totalMemory() - runtime.freeMemory();
		runtime.gc();
		long memUseGC = runtime.totalMemory() - runtime.freeMemory();

		System.out.print(", " + (eTime - sTime) / 1000000000.0);
		System.out.print(", " + (eTime - sTime) / 1000000000.0 / 60.0);
		System.out.print(", " + memUse / 1024.0 / 1024.0);
		System.out.print(", " + memUseGC / 1024.0 / 1024.0);

		System.out.print(", " + plantedMotif + "," + foundMotifs);
		

	}

	public static void readInput(String filename) throws Exception {
		Scanner input = new Scanner (new File (filename ));
		t = input.nextInt();		
		n = input.nextInt();		
		l = input.nextInt();
		d = input.nextInt();
		plantedAlignment = new int[t];
		//plantedMotif = new char[l];
		seqS = new char[t][n];
		int a = 0;
		while (input.hasNextInt()) {	
			plantedAlignment[a] = input.nextInt();
			a++;
		}		
		plantedMotif = input.next().toUpperCase();
		//plantedMotif = motif.toCharArray();
		
		int n = 0;
		while (input.hasNext()) {
			String seq = input.next().toUpperCase();
			seqS[n] = seq.toCharArray();
			n++;
		}
		input.close();
	}

//========================================================================================================
// COLLECT CANDIDATE MOTIFS: intersect the d-neighborhoods of sequences 0 to tPrime.
//========================================================================================================

	public static void collectCandidateMotifs() throws Exception {
		generateNeigbhorhood(0);		
		candidateMotifs = currNeighborhood;
		for(int i=1; i < tPrime; i++) {
			generateNeigbhorhood(i);
			for(int j=0; j < pmt; j++)
				candidateMotifs[j] &= currNeighborhood[j];
		}
	}

	public static void generateNeigbhorhood(int s) throws Exception {
		// System.out.println("genNeigbhorhood(" + s + ")");
		currNeighborhood = new int[pmt];
		char[] currSeq = seqS[s];

		long mapping = 0;
		for(int i=0; i < l; i++) {
			char c = currSeq[i];
			int base = 0;
			switch(c) {
				case 'C': base=1; break;
				case 'G': base=2; break;
				case 'T': base=3; break;
			}
			mapping = (mapping << 2) + base;
		}
		addNeighbors(mapping, 0, d);

		/*if( s == 0 ) {
			System.out.println("Mapping: " + Long.toBinaryString(mapping) );
			System.out.println("Index: " + (mapping >>> 5));
			System.out.println("Offset: " + (mapping & 31));
			System.out.printf("\t\t\t0123456789ABCDEFGHIJKLMNOPQRSTUVW\n");
			for(int i = 0; i < pmt; i++) {
				if( currNeighborhood[i] != 0 ) {
					String row = Long.toBinaryString(currNeighborhood[i]);
					while( row.length() < 64 ) {
						row = "0" + row;
					}
					System.out.printf("%16d\t%s\n", i, row);
				}
			}
		}*/

		// System.out.println("addNeighbors(" + Long.toBinaryString(mapping) + "," + d + ")");	

		for(int i=l; i < n; i++) {
			char c = currSeq[i];
			int base = 0;
			switch(c) {
				case 'C': base=1; break;
				case 'G': base=2; break;
				case 'T': base=3; break;
			}
			mapping = ((mapping & mask) << 2) + base;
			addNeighbors(mapping, 0, d);
			// System.out.println("addNeighbors(" + Long.toBinaryString(mapping) + "," + d + ")");
		}
	}

	public static void addNeighbors(long mapping, int start, int allow_d) throws Exception {
		int shift=(l-start)*2;
		for(int k=start; k < l; ++k) {
			shift-=2;
			long alt1 = mapping ^ (((long) 1) << shift);
			long alt2 = mapping ^ (((long) 2) << shift);
			long alt3 = mapping ^ (((long) 3) << shift);

			int index = (int)(alt1 >>> 5);
			currNeighborhood[index] |= (1 << (alt1 & 31));
			
			index = (int)(alt2 >>> 5);
			currNeighborhood[index] |= (1 << (alt2 & 31));
			
			index = (int)(alt3 >>> 5);
			currNeighborhood[index] |= (1 << (alt3 & 31));
			
			if(allow_d > 1) {
				addNeighbors(alt1, k+1, allow_d-1);
				addNeighbors(alt2, k+1, allow_d-1);
				addNeighbors(alt3, k+1, allow_d-1);
			}
		}
	}

//========================================================================================================
// TRANSFORM L-MER SEQUENCES: create mappings for every l-mer in sequences tPrime+1 to t.
//========================================================================================================

	public static void transformLmerSequences(int tPrime) {
		lmerMappings = new long[t][nl1];
		for(int i=tPrime; i < t; i++) {
			char[] currSeq = seqS[i];
			long mapping = 0;

			for(int j=0; j < l; j++) {
				char c = currSeq[j];
				int base = 0;
				switch(c) {
					case 'C': base=1; break;
					case 'G': base=2; break;
					case 'T': base=3; break;
				}
				mapping = (mapping << 2) + base;
			}
			lmerMappings[i][0] = mapping;

			int k=0;
			for(int j=l; j < n;j++) {
				char c = currSeq[j];
				int base = 0;
				switch(c) {
					case 'C': base=1; break;
					case 'G': base=2; break;
					case 'T': base=3; break;
				}
				mapping = ((mapping & mask) << 2) + base;
				lmerMappings[i][++k]=mapping;
			}

		}
	}

//========================================================================================================
// SEARCH MOTIF: check each element of candidateMotifs if present in sequences tPrime+1 to t.
//========================================================================================================

	public static void searchMotif() throws Exception {
		int value, numMotifs = 0;
		foundMotifs = "";
		for(int i=0; i < pmt; i++) {
			if( (value = candidateMotifs[i]) == 0)
				continue;

			// System.out.println("Nonzero: candidateMotifs[" + i + "]\t= "
			// 	+ Long.toBinaryString(candidateMotifs[i]));
			long base = ((long) i) << 5;
			for(int j=0; j < 32; j++) {
				if( (value & 1) != 0) {
					long candidate = base + j;
					if( isMotif(candidate, tPrime) ) {
						foundMotifs += " " + decode(candidate, l);
						numMotifs++;
					}
				}
				value = value >>> 1;
			}
		}
	}

	public static boolean isMotif(long mapping, int tPrime) throws Exception {
		for(int i=tPrime; i < t; i++) {
			boolean found = false;
			for(int j=0; j < nl1; j++) {
				long lmer = lmerMappings[i][j];
				int hd = computeHD(mapping, lmer);
				if( hd <= d ) {
					found=true;
					break;
				}
			}
			if(!found) {
				return false;
			}
		}
		return true;
	}

	public static int computeHD(long lmer1, long lmer2) throws Exception {
		int distance=0;
		long result = lmer1 ^ lmer2;
		for(int i=0; i < l; i++) {
			if( (result & 3) != 0 )
				distance++;
			result = result >>> 2;
		}
		return distance;
	}

	public static String decode(long mapping, int strlen) throws Exception {
		String decoding = "";
		for(int i=0; i < strlen; i++) {
			int base = (int) mapping & 3;
			switch(base) {
				case 0: decoding = "A" + decoding; break;
				case 1: decoding = "C" + decoding; break;
				case 2: decoding = "G" + decoding; break;
				case 3: decoding = "T" + decoding; break;
			}
			mapping = mapping >>> 2;
		}
		return decoding;
	}

}