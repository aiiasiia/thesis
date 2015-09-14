/*
 * CompatiblePairs.h
 *
 *  Created on: Nov 19, 2012
 *      Author: marius
 */

#ifndef COMPATIBLEPAIRS_H_
#define COMPATIBLEPAIRS_H_
#include "CompressedLmers.h"
#include "utils.h"
#include <cstring>

template<class indexType>
class CompatiblePairs {
public:
	CompatiblePairs(int nRows, int *rSize, indexType **rowItem, int L, int twoD,
			char *minLmerPtr, int nLmers, CompressedLmers *compressedLmers) :
			nLmers(nLmers) {
		if (compressedLmers != NULL)
			prepareCompatiblePairs<true>(nRows, rSize, rowItem, L, twoD,
					minLmerPtr, compressedLmers, nLmers);
		else
			prepareCompatiblePairs<false>(nRows, rSize, rowItem, L, twoD,
					minLmerPtr, compressedLmers, nLmers);
	}

	virtual ~CompatiblePairs() {
		deleteDataStructures();
	}

	int **pairOk;
private:

	template<bool areLmersCompressed>
	void precomputeCompatPairs(int *rowSize, indexType **rowItem, int n, int L,
			int twoD, char *minLmerPtr, CompressedLmers *compressedLmers) {
		for (int i = 0; i < n; ++i) {
			int nItemsA = rowSize[i];
			indexType *itemA = rowItem[i];
			for (int j = i + 1; j < n; ++j) {
				int nItemsB = rowSize[j];
				indexType *itemB = rowItem[j];
				for (int k = 0; k < nItemsA; ++k) {
					int aInd = itemA[k];
					int *p = pairOk[aInd];
					for (int l = 0; l < nItemsB; ++l) {
						int bInd = itemB[l];
						int h;
						if (areLmersCompressed) {
							h = compressedLmers->HamDistStoredLmers(aInd, bInd);
						} else {
							h = HamDist(minLmerPtr + aInd, minLmerPtr + bInd,
									L);
						}
						/*
						int h2 = HamDist(minLmerPtr + aInd, minLmerPtr + bInd,
								L);
						if (h != h2) {
							cout << "HORROR" << endl;
							cout << "HORROR" << endl;
							cout << "HORROR" << endl;
							cout << "HORROR" << endl;
							cout << "HORROR" << endl;
							cout << "HORROR" << endl;
							cout << "HORROR" << endl;
							cout << "HORROR" << endl;
						}*/
						if (h <= twoD) {
							setBit(p, bInd);
							setBit(pairOk[bInd], aInd);
						}
					}
				}
			}
		}
	}

	template<bool areLmersCompressed>
	void prepareCompatiblePairs(int n, int *rSize, indexType **rowItem, int L,
			int twoD, char *minLmerPtr, CompressedLmers *compressedLmer,
			int nLmers) {
		int bytes = sizeof(**pairOk);
		int bits = 8 * bytes;
		int nAllocLmersCol = (nLmers / bits) + 1;

		clock_t start = clock();
//		cout << "Pair size " << nLmers << " by " << nAllocLmersCol << endl;
		allocate(pairOk, nLmers, nAllocLmersCol);
		for (int i = 0; i < nLmers; ++i)
			memset(pairOk[i], 0, nAllocLmersCol * bytes);

		precomputeCompatPairs<areLmersCompressed>(rSize, rowItem, n, L, twoD,
				minLmerPtr, compressedLmer);
		clock_t end = clock();
		long time = (end - start) * 1000 / CLOCKS_PER_SEC;
		cout << "Time to preprocess pairs " << time << "ms" << endl;
	}

	void deleteDataStructures() {
		deAllocate(pairOk, nLmers);
	}

	int nLmers;
};

#endif /* COMPATIBLEPAIRS_H_ */
