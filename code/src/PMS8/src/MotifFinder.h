/*
 * MotifFinder.h
 *
 *  Created on: Nov 20, 2012
 *      Author: marius
 */

#ifndef MOTIFFINDER_H_
#define MOTIFFINDER_H_
#include "CompressedLmers.h"
#include "CompatiblePairs.h"
#include "StringSet.h"
#include "ClosestSub.h"
//#include <omp.h>

template<bool compressLmers, class indexType>
class MotifFinder {
public:
	MotifFinder(StringSet& s, MotifConfig *mc) :
			mc(mc) {
		minLmerPtr = s.memStart;

		compressedLmers = NULL;
		if (compressLmers)
			compressedLmers = new CompressedLmers(minLmerPtr, s.totalLength,
					mc->L, mc->sigmaLen);

		int nLmers = s.totalLength - mc->L + 1;
		rowSize = s.range;
		rowItem = createIndexMatrix(s);
		compatiblePairs = new CompatiblePairs<indexType>(mc->n, rowSize,
				rowItem, mc->L, 2 * mc->d, minLmerPtr, nLmers, compressedLmers);

		ms = new ClosestSub<compressLmers, indexType>(mc, rowSize, rowItem,
				minLmerPtr, compressedLmers, compatiblePairs);
	}

	void processIndex(int i) {
		ms->processIndex(i);
	}

	set<MyString> getMotifs() {
		return ms->getMotifs();
	}

	~MotifFinder() {
		if (compressLmers)
			delete compressedLmers;
		deAllocate(rowItem, mc->n);
		delete compatiblePairs;
		delete ms;
	}
private:
	indexType **createIndexMatrix(StringSet& s) {
		indexType **rowItem = new indexType*[s.n];
		for (int i = 0; i < s.n; ++i) {
			int nItems = s.range[i];
			indexType *item = new indexType[nItems];
			rowItem[i] = item;
			int l = s.s[i] - s.memStart;
			for (int c = 0; c < nItems; ++c)
				item[c] = l + c;
		}
		return rowItem;
	}

	MotifConfig *mc;
	char * minLmerPtr;
	int *rowSize;
	indexType **rowItem;
	CompatiblePairs<indexType> *compatiblePairs;
	ClosestSub<compressLmers, indexType> *ms;
	CompressedLmers *compressedLmers;
};

#endif /* MOTIFFINDER_H_ */
