/**
  * EMS-GT_64.java
  * Solves the (l,d) planted motif problem.
  * Bit-array width is 64 bits (Float type).
  * Optimized with block masking.
  *
  * @author Aia Sia, Julieta Nabos
  * @version 1.0 9/09/3025
  */

import java.io.*;
import java.util.*;

public class EMS_GT_64 {

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

	// for EMS-GT computations
	static long mask, prefixMask, suffixMask;
	static int nl1;
	static int pmt;
	static long[] currNeighborhood;
	static long[] candidateMotifs;
	static long[][] lmerMappings;

	// for masking
	static int blockDegree = 5;
	static int lmersInBlock;
	static int rowsInBlock;
	static int[][][] blockMasks;
	static int[][] currBlockMasks;
	static long prefix;
	static long suffix;
	static int currBlockRow;
	static int currBlockCol;
	
//========================================================================================================
// MAIN METHOD
//========================================================================================================

	public static void main(String args[]) throws Exception {
		if( args.length > 0 ) {
			inputFileName = args[0];
			if( args.length > 1 ) {
				blockDegree = Integer.parseInt(args[1]);
			}
		}

		readInput( inputFileHeader + inputFileName );

		// compute preliminaries
		mask = (((long) 1) << 2*(l-1) ) - 1;
		prefixMask = (((long) 1) << 2*(l-blockDegree-1) ) - 1;
		suffixMask = (((long) 1) << 2*(blockDegree-1) ) - 1;
		nl1 = n - l + 1;
		pmt = (int) (( (long)1 << (2*l) ) >>> 6); // (4^l)/64

		System.out.print("\n" + inputFileName);
		
		// start EMS-GT
		long sTime = System.nanoTime();
		generateBlockMasks();
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

		System.out.println(", " + plantedMotif + "," + foundMotifs);
	}

	public static void generateBlockMasks() throws Exception {
		lmersInBlock = 1 << (2 * blockDegree);		// 4 ^ blockDegree
		rowsInBlock = lmersInBlock >> 6;			// 4 ^ blockDegree / 64
		blockMasks = new int[lmersInBlock][blockDegree - 1][rowsInBlock];
		for(int i=0; i < lmersInBlock; i++) {
			for(int k=0; k < blockDegree - 1; k++) {
				Arrays.fill( blockMasks[i][k], 0);
			}
			for(int row=0; row < rowsInBlock; row++) {
				for(int col=63; col > -1; col--) {
					int distance = computeHD(i, row*64+col);
					for(int k=0; k < blockDegree - 1; k++) {
						if(distance <= k+1) {
							blockMasks[i][k][row]++;
						}
						if(col > 0) {
							blockMasks[i][k][row] = blockMasks[i][k][row] << 1;
						}
					}
				}

			}
		}

		/*System.out.println("Done.\n\nSample block mask: centered at row 14, col 12, distance <= 4: ");
		currBlockMasks = blockMasks[14*32+12];
		for(int row=0; row < rows; row++) {
			System.out.printf("%32s\n", Long.toBinaryString(currBlockMasks[row][4]));
		}*/
	}

	public static void readInput(String filename) throws Exception {
		Scanner input = new Scanner (new File (filename ));
		t = input.nextInt();		
		n = input.nextInt();		
		l = input.nextInt();
		d = input.nextInt();
		plantedAlignment = new int[t];
		seqS = new char[t][n];
		int a = 0;
		while (input.hasNextInt()) {	
			plantedAlignment[a] = input.nextInt();
			a++;
		}		
		plantedMotif = input.next().toUpperCase();

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
		generateNeighborhood(0);		
		candidateMotifs = currNeighborhood;
		for(int i=1; i < tPrime; i++) {
			generateNeighborhood(i);
			for(int j=0; j < pmt; j++)
				candidateMotifs[j] &= currNeighborhood[j];
		}
	}

	public static void generateNeighborhood(int s) throws Exception {
		// System.out.println("genNeigbhorhood(" + s + ")");
		currNeighborhood = new long[pmt];
		//Arrays.fill(currNeighborhood,0);
		char[] currSeq = seqS[s];

		prefix = 0; 					 	// first < l-blockDegree > characters of l-mer
		suffix = 0;  						// last  < blockDegree > characters of l-mer
		for(int i=0; i < l; i++) {
			char c = currSeq[i];
			int base = 0;
			switch(c) {
				case 'C': base=1; break;
				case 'G': base=2; break;
				case 'T': base=3; break;
			}
			if( i < l - blockDegree )
				prefix = (prefix << 2) + base;
			else
				suffix = (suffix << 2) + base;
		}

		// housekeeping: set blockOffsets, currBlockMasks
		currBlockRow = (int) (suffix / 64);
		currBlockCol = (int) (suffix % 64);
		currBlockMasks = blockMasks[(int)suffix];
		int blockStart = (int) (prefix << (2*blockDegree - 6));
		for(int offset=0; offset < rowsInBlock; offset++) {
			if(d >= blockDegree)
				currNeighborhood[blockStart+offset] = 0xffff_ffff_ffff_ffffL;
			else
				currNeighborhood[blockStart+offset] = currBlockMasks[d-1][offset];
		}
		addNeighbors(prefix, 0, d);
		
		for(int i=l; i < n; i++) {			
			prefix = (prefix & prefixMask) << 2;
			suffix = (suffix & suffixMask) << 2;

			char c = currSeq[i-blockDegree]; 	// next char for prefix
			switch(c) {
				case 'C': prefix+=1; break;
				case 'G': prefix+=2; break;
				case 'T': prefix+=3; break;
			}

			c = currSeq[i]; 					// next char for suffix
			switch(c) {
				case 'C': suffix+=1; break;
				case 'G': suffix+=2; break;
				case 'T': suffix+=3; break;
			}
			
			// housekeeping: set blockOffsets, currBlockMasks
			currBlockRow = (int) suffix / 64;
			currBlockCol = (int) suffix % 64;
			currBlockMasks = blockMasks[(int)suffix];
			blockStart = (int) prefix << (2*blockDegree - 6);
			for(int offset=0; offset < rowsInBlock; offset++) {
				if(d >= blockDegree)
					currNeighborhood[blockStart+offset] = Long.MAX_VALUE;
				else
					currNeighborhood[blockStart+offset] |= currBlockMasks[d - 1][offset];
			}
			addNeighbors(prefix, 0, d);
		}
	}

	public static void addNeighbors(long prefix, int start, int d) throws Exception {
		int shift = (l-blockDegree-start)*2;
		for(int i=start; i < l-blockDegree; ++i) {
			shift -= 2;
			long alt1 = prefix ^ (((long) 1) << shift);
			long alt2 = prefix ^ (((long) 2) << shift);
			long alt3 = prefix ^ (((long) 3) << shift);

			int multByRow = 2*blockDegree - 6;
			int blockStart1 = (int) alt1 << multByRow;
			int blockStart2 = (int) alt2 << multByRow;
			int blockStart3 = (int) alt3 << multByRow;

			// masking part
			int allow_d = d - 1;
			if( allow_d >= blockDegree ) {		// all 1's
				for(int offset=0; offset < rowsInBlock; offset++) {
					currNeighborhood[blockStart1 + offset] = Long.MAX_VALUE;
					currNeighborhood[blockStart2 + offset] = Long.MAX_VALUE;
					currNeighborhood[blockStart3 + offset] = Long.MAX_VALUE;
				}
			}
			else if( allow_d > 0 ) {			// select a mapping from 1 to blockDegree-1
				for(int offset=0; offset < rowsInBlock; offset++) {
					currNeighborhood[blockStart1 + offset] |= currBlockMasks[allow_d - 1][offset];
					currNeighborhood[blockStart2 + offset] |= currBlockMasks[allow_d - 1][offset];
					currNeighborhood[blockStart3 + offset] |= currBlockMasks[allow_d - 1][offset];
				}
			}
			else {			// only [currBlockRow][currBlockCol] = 1
				currNeighborhood[blockStart1 + currBlockRow] |= 1 << (currBlockCol);
				currNeighborhood[blockStart2 + currBlockRow] |= 1 << (currBlockCol);
				currNeighborhood[blockStart3 + currBlockRow] |= 1 << (currBlockCol);
			}

			// recursive call
			if(allow_d > 0) {
				addNeighbors(alt1, i+1, allow_d);
				addNeighbors(alt2, i+1, allow_d);
				addNeighbors(alt3, i+1, allow_d);
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
		long value = 0;
		int numMotifs = 0;
		foundMotifs = "";
		for(int i=0; i < pmt; i++) {
			if( (value = candidateMotifs[i]) == 0)
				continue;

			/*System.out.println("Nonzero: candidateMotifs[" + i + "]\t= "
				+ Long.toBinaryString(candidateMotifs[i]));*/
			long base = ((long) i) << 6;
			for(int j=0; j < 64; j++) {
				if( (value & 1) != 0) {
					long candidate = base + j;
					if( isMotif(candidate, tPrime) ) {
						foundMotifs += " " + decode(candidate, l);
						numMotifs++;
					}
				}
				value = value >> 1;
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