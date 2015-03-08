#ifndef OUTBOUND_H
#define OUTBOUND_H

#include "present.h"
#include <set>

void test_outbound_trail(int, int, int);
void test_hull(int, int, int, int, int);
void kage();

class outbound {
private:	
	int in_pos, out_pos, inbound_rounds;
	int num_rounds;
	int sat_pairs;
	void verify_L_in();
	void handle_pair(u64, u64, u64);

public:
	std::set<std::pair<u64, key_state>> L_out, L_in;
	outbound(std::set<std::pair<u64, key_state>> &, int, int, int, int);
	
	void run_phase();
	
};

#endif
