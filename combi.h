#ifndef COMBI_ITERATOR_H
#define COMBI_ITERATOR_H

class CombiIterator {
	private:
		int i, S; // Number of elements in combination
		int next_perm();

	public:
		int *combination;	// For holding a specific combination
		int *limit;			// For holding at limit[i] the number of possible values combination[i] can take on

		CombiIterator(int);
		~CombiIterator();		
		void reset();
		int has_next();
};

#endif

