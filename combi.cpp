#include "combi.h"
#include <stdlib.h>
#include <string.h>

CombiIterator::~CombiIterator() {
	free(combination);
	free(limit);
}

CombiIterator::CombiIterator(int S) {
	combination = (int*) calloc(S, sizeof(int));
	limit = (int*) calloc(S, sizeof(int));
	this->S = S;
	combination[S-1] = -1;
}

int CombiIterator::next_perm() {
	// look at each index
	for (i = S - 1; i >= 0; --i) {
		combination[i]++;
		if (combination[i] < limit[i])
			return 1;
		combination[i] = 0;
	}
	return 0;
}

void CombiIterator::reset() {
	memset(combination, 0, S * sizeof(int));
	combination[S-1] = -1;
}

int CombiIterator::has_next() {
	return next_perm();
}
