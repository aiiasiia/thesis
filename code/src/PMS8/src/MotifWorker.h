/*
 * MotifWorker.h
 *
 *  Created on: Nov 24, 2012
 *      Author: marius
 */

#ifndef MOTIFWORKER_H_
#define MOTIFWORKER_H_
#include <utility>
#include <algorithm>
#include <vector>
#include <set>
#include <string>
#include <cstring>
#include <cmath>
#include "utils.h"
#include "MyString.h"
#include "StringSet.h"
#include "MotifFinder.h"
using namespace std;

class MotifWorker {
	typedef uint32 indexType;
public:
	MotifWorker() :
			totalLen(0) {
		s = NULL;
		mf = NULL;
	}

	~MotifWorker() {
		if (s != NULL)
			delete s;
		if (mf != NULL)
			delete mf;
	}

	void init(int L, int d, int t = -1, int nPrime = -1) {
		mc.L = L;
		mc.d = d;
		mc.t = t;
		mc.nPrime = nPrime;
	}

	int *readAndEncodeInput(char *fileName, int& totalMsgWords) {
		vector<string> strings;
		int totalLen = parseInput(fileName, strings);
		if (totalLen < 0)
			return NULL;
		else {
			cout << "Read " << totalLen << " chars from file " << fileName
					<< endl;
			return encodeData(strings, totalMsgWords);
		}
	}

	void loadConfigFromBuffer(int *b) {
		// read alphabet
		char *bc = (char*) b;
		int n = b[0];
		int sigmaLen = b[n + 1];
		int pos = (n + 2) * sizeof(int);
		sigma = string(bc + pos, sigmaLen);
		pos += sigmaLen;

		// read config
		memcpy(&mc, bc + pos, sizeof(MotifConfig));
	}

	void printConfig() {
		cout << " Alphabet " << sigma << "  n=" << mc.n << " L=" << mc.L
				<< " d=" << mc.d << " vstride=" << mc.nPrime << endl;
	}

	void loadStringsFromBuffer(int *b) {
		s = decodeStrings(b, mc);
		mf = new MotifFinder<true, indexType>(*s, &mc);
	}

	int getNLmersInFirstString() {
		return s->range[0];
	}

	void process(int index) {
		mf->processIndex(index);
	}

	set<MyString> getMotifs() {
		return mf->getMotifs();
	}

	int *allocateMotifBuffer(int nMotifs, int& sz) {
		sz = 1 + ceil(((double) nMotifs * mc.L) / sizeof(int));
		return new int[sz];
	}

	int *encodeMotifs(set<MyString>& motifs, int nMotifs, int& requiredMem) {
		requiredMem = ceil(((double) nMotifs * mc.L) / sizeof(int)) + 1;
		int *buf = new int[requiredMem];
		buf[0] = nMotifs;
		char *bc = ((char *) buf) + sizeof(int);
		set<MyString>::iterator it = motifs.begin();
		for (int i = 0; i < nMotifs; ++i) {
			memcpy(bc, it->s, it->L);
			bc += mc.L;
			++it;
		}
		return buf;
	}

	int decodeMotifs(int *buf, set<MyString>& output) {
		int nm = buf[0];
		char *bc = ((char*) buf) + sizeof(int);
		for (int j = 0; j < nm; ++j, bc += mc.L) {
			output.insert(MyString(bc, mc.L));
		}
		return nm;
	}

	void printMotifs(set<MyString>& motifs) {
		std::cout << "Motifs:" << endl;
		for (set<MyString>::iterator it = motifs.begin(); it != motifs.end();
				++it) {
			std::cout << decodeString(it->s, mc.L, sigma)
					<< std::endl;
		}
	}

	int getScore(MyString& motif) {
		return 0;
	}

	static bool compare_pairs(pair<int, MyString> a, pair<int, MyString> b) {
		return a.first < b.first;
	}

	vector<MyString> rankMotifs(set<MyString>& motifs) {
		vector<MyString> rankedMotifs(motifs.begin(), motifs.end());
		return rankedMotifs;

		cout << "Ranking motifs" << endl;
		vector<pair<int, MyString> > pairs;
		for (set<MyString>::iterator it = motifs.begin(); it != motifs.end();
				++it) {
			MyString str = *it;
			int score = getScore(str);
			pairs.push_back(make_pair(score, str));
		}
		cout<<"Sorting paris" << endl;
		sort(pairs.begin(), pairs.end(), compare_pairs);
		cout<<"Sorted paris" << endl;

		vector<MyString> result;
		for (vector<pair<int, MyString> >::iterator it = pairs.begin();
				it != pairs.end(); ++it) {
			result.push_back(it->second);
		}
		cout<<"Collected paris " << result.size() << endl;
		return result;
	}

	void getMyChunk(int originalFirstStringLen, int myRank, int nProcessors,
			int& beginIndex, int& nLmers) {
		int totalLmers = (originalFirstStringLen - mc.L + 1);
		nLmers = totalLmers / nProcessors;
		beginIndex = myRank * nLmers;
		if (myRank < totalLmers % nProcessors) {
			nLmers++;
			beginIndex += myRank;
		} else
			beginIndex += totalLmers % nProcessors;
	}

	MotifConfig mc;
private:

	StringSet *decodeStrings(int *b, MotifConfig& mc) {
		int pos = (mc.n + 2) * sizeof(int) + mc.sigmaLen + sizeof(MotifConfig);
		char *bc = ((char*) b) + pos;
		return new StringSet(mc.n, mc.L, b + 1, bc);
	}

	void populateThresholds(int totalLen, MotifConfig& mc) {
		double avgStringLen = 1.0 * totalLen / mc.n;
//		for (int d = 1; d <= 50; ++d) {
//			mc.d = d;

		double logSigma = log2(mc.sigmaLen);
		double minStackSizeSquared = max(4.0,
				(2 * mc.d - 2) * logSigma - log2(avgStringLen) + 8);
		double t = sqrt(minStackSizeSquared);
		if (mc.t > 0)
			t = mc.t;
		mc.t = min(mc.n, (int) t);

		double nPrime = t + mc.n / 4 - log2(t);
		if (mc.nPrime > 0)
			nPrime = mc.nPrime;
		mc.nPrime = max(mc.t, min(mc.n, (int) nPrime));

		cout << "t=" << mc.t << " n'=" << mc.nPrime << endl;
//		}
		mc.firstThreshold = min(mc.nPrime - 1, 2);
		mc.bruteForceThreshold = 30;
	}

	int parseInput(char *fileName, vector<string>& strings) {
		if (readFasta(fileName, strings)) {
			// !!!!!!!!!!!!!!!!
//			strings.erase(strings.begin(),strings.begin()+100);
//			strings.erase(strings.begin()+750,strings.end());
			// !!!!!!!!!!!!!!!
			mc.n = strings.size();
			sigma = getAlphabet(strings);
			mc.sigmaLen = sigma.length();
			encodeStrings(strings, sigma);
			totalLen = sumLengths(strings);
			cout << "n=" << mc.n << " L=" << mc.L << " d=" << mc.d
					<< " totalStringLength=" << totalLen << endl;
			cout << "Alphabet " << sigma << endl;
			populateThresholds(totalLen, mc);
			return totalLen;
		} else {
			return -1;
		}
	}

	int *encodeData(vector<string>& strings, int& totalMsgWords) {
		totalMsgWords = mc.n + 3
				+ (sigma.length() + sizeof(MotifConfig) + totalLen)
						/ sizeof(int);
		int *b = new int[totalMsgWords];

		// write number of strings and string lengths
		b[0] = mc.n;
		for (int i = 0; i < mc.n; ++i)
			b[i + 1] = strings[i].length();

		// write alphabet
		b[mc.n + 1] = sigma.length();
		char *bc = (char*) b;
		int pos = (mc.n + 2) * sizeof(int);
		memcpy(bc + pos, sigma.c_str(), sigma.length());
		pos += sigma.length();

		// write configuration
		memcpy(bc + pos, (char*) &mc, sizeof(MotifConfig));
		pos += sizeof(MotifConfig);

		// write strings
		for (int j = 0; j < mc.n; ++j) {
			int len = strings[j].length();
			memcpy(bc + pos, strings[j].c_str(), len);
			pos += len;
		}
		return b;
	}

	StringSet *s;
	MotifFinder<true, indexType> *mf;
	string sigma;
	int totalLen;
};

#endif /* MOTIFWORKER_H_ */
