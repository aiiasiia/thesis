/*
 * CompressedLmers.h
 *
 *  Created on: Nov 19, 2012
 *      Author: marius
 */

#ifndef COMPRESSEDLMERS_H_
#define COMPRESSEDLMERS_H_
#include "utils.h"

//template<class comprWordType, class comprPreprocType>
class CompressedLmers {
public:
	typedef uint16 comprWordType;
	typedef uint16 comprPreprocType;
	CompressedLmers(char *s, int sLength, int L, int sigmaLen) :
			L(L) {
		compressAllLmers(s, sLength, sigmaLen);
	}
	virtual ~CompressedLmers() {
		deleteCompressedLmers();
	}

	int HamDistStoredLmers(const int aInd, const int bInd) {
		comprWordType *ca = compressedLmer + aInd;
		comprWordType *cb = compressedLmer + bInd;
		int d = 0;
		int i = 0;
		comprWordType a = ca[0] ^ cb[0];
		for (; i + symbolsPerWord <= L;) {
			d += compressedMismatches[a];
			i += symbolsPerWord;
			a = ca[i] ^ cb[i];
		}
//		if (i < L)
			d += compressedMismatches[a >> blankBitsInLastLmerWord];
		return d;
	}

	int HamDistGeneratedVsStored(const comprWordType *ca,
			const int storedItemIndex) {
		comprWordType *cb = compressedLmer + storedItemIndex;
		int d = 0;
		int i = 0, j = 0;
		comprWordType a = ca[0] ^ cb[0];
		for (; i + symbolsPerWord <= L;) {
			d += compressedMismatches[a];
			++j;
			i += symbolsPerWord;
			a = ca[j] ^ cb[i];
		}
//		if (i < L)
			d += compressedMismatches[a >> blankBitsInLastLmerWord];
		return d;
	}

	void compressLmer(const char *l, comprWordType *cl) {
		int i = 0;
		for (; i + symbolsPerWord <= L; i += symbolsPerWord) {
			comprWordType w = 0;
			for (int j = 0; j < symbolsPerWord; ++j)
				w = (w << bitsPerSymbol) | *l++;
			*cl++ = w;
		}
		if (i < L) {
			comprWordType w = 0;
			for (; i < L; ++i)
				w = (w << bitsPerSymbol) | *l++;
			*cl = w << blankBitsInLastLmerWord;
		}
	}

private:

	int countNonZeroGroups(comprWordType w, comprWordType lastSymbolMask,
			int bitsPerSymbol, int symbolsPerWord) {
		int d = 0;
		for (int i = 0; i < symbolsPerWord; ++i) {
			d += (w & lastSymbolMask) != 0;
			w = w >> bitsPerSymbol;
		}
		return d;
	}

	void compressAllLmers(char *s, int sLength, int sigmaLen) {
		bitsPerSymbol = bitsFor(sigmaLen - 1);
		unsigned bitsPerWord = 8 * sizeof(comprWordType);
		symbolsPerWord = bitsPerWord / bitsPerSymbol;
		blankBitsInLastLmerWord = (symbolsPerWord - L % symbolsPerWord)
				* bitsPerSymbol;

		extractLastSymbols = new comprWordType[symbolsPerWord + 1];
		comprWordType w = ~0;
		for (int i = 0; i <= symbolsPerWord; ++i) {
			extractLastSymbols[i] = ~w;
			w = w << bitsPerSymbol;
		}

		w = 0;
		for (int i = 0; i < symbolsPerWord; ++i) {
			w = (w << bitsPerSymbol) | s[i];
		}
		comprWordType left0Mask = extractLastSymbols[symbolsPerWord];
		compressedLmer = new comprWordType[sLength];
		for (int i = 0, e = sLength - symbolsPerWord; i < e; ++i) {
			compressedLmer[i] = w;
			w = ((w << bitsPerSymbol) | s[i + symbolsPerWord]) & left0Mask;
		}
		for (int i = sLength - symbolsPerWord; i < sLength; ++i) {
			compressedLmer[i] = w;
			w = (w << bitsPerSymbol) & left0Mask;
		}

		int symbolsPerPreprocessedWord;
		if (8 * sizeof(comprPreprocType) < bitsPerWord)
			symbolsPerPreprocessedWord = (symbolsPerWord + 1) / 2;
		else
			symbolsPerPreprocessedWord = symbolsPerWord;
		int bitsPerPreprocessedWord = symbolsPerPreprocessedWord
				* bitsPerSymbol;
		int sz = 1 << bitsPerPreprocessedWord;
		compressedMismatches = new char[sz];
		comprWordType rightMask = extractLastSymbols[1];
		for (int i = 0; i < sz; ++i) {
			compressedMismatches[i] = countNonZeroGroups(i, rightMask,
					bitsPerSymbol, symbolsPerPreprocessedWord);
		}
	}

	void deleteCompressedLmers() {
		delete[] compressedLmer;
		delete[] extractLastSymbols;
		delete[] compressedMismatches;
	}

	const int L;
	int bitsPerSymbol;
	int symbolsPerWord;
	int blankBitsInLastLmerWord;
	comprWordType *compressedLmer;
	comprWordType *extractLastSymbols;
	char *compressedMismatches;
};

#endif /* COMPRESSEDLMERS_H_ */
