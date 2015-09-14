/*
 * StringSet.h
 *
 *  Created on: Jul 4, 2012
 *      Author: marius
 */

#ifndef STRINGSET_H_
#define STRINGSET_H_
#include <string>
#include <vector>
using namespace std;
#include "utils.h"

class StringSet {
public:
	StringSet(int n, int L, int *length, char *memStart) :
			externalMem(true), n(n), memStart(memStart), length(length) {
		s = new char*[n];
		char *m = memStart;
		range = new int[n];
		for (int i = 0; i < n; ++i) {
			s[i] = m;
			m += length[i];
			range[i] = max(0, length[i] - L + 1);
		}
		totalLength = m - memStart;
	}

	StringSet(vector<string>& strings) :
			externalMem(false) {
		n = strings.size();
		s = pack(strings);
		memStart = s[0];
		length = new int[n];
		totalLength = 0;
		for (int i = 0; i < n; ++i) {
			length[i] = strings[i].length();
			totalLength += length[i];
		}
		range = new int[n];
	}

	~StringSet() {
		if (!externalMem) {
			delete[] length;
			delete[] memStart;
		}
		delete[] range;
		delete[] s;
	}

	bool externalMem;
	int n;
	int totalLength;
	char *memStart;
	int *length;
	char **s;
	int *range;
};

#endif /* STRINGSET_H_ */
