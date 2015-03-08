#include "outbound.h"
#include "present.h"
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <random>
#include <inttypes.h>
#include <iostream>

#define DEBUG 1
#define VERIFY 1

outbound::outbound(std::set<std::pair<u64, key_state>> &start_pairs, int _ip, int _op, int _ib_rounds, int _num_rounds) {
	L_out = start_pairs;
	in_pos = _ip;
	out_pos = _op;
	inbound_rounds = _ib_rounds;
	num_rounds = _num_rounds;
	sat_pairs = 0;
#if DEBUG
	printf(">> Initialised outbound phase with #L_out = %d starting points\n\n", L_out.size());
#endif
}

// For verifying that the filtered set L_in indeed satisfies the
// inbound relation as well
void outbound::verify_L_in() {
	std::set<std::pair<u64, key_state>>::iterator it;

	u64 in_state, state, lk, rk;
	for (it = L_in.begin(); it != L_in.end(); ++it) {
		in_state = (*it).first;
		state = in_state;
		lk = (*it).second.left;
		rk = (*it).second.right;

		// Encrypt
		for (int r = 0; r < (inbound_rounds+num_rounds); ++r) {
			encrypt_one_round(state, lk);
			keysched_round(lk, rk, (r+1));
		}
		state ^= lk; // Whitening key

		// Check relation(s). This is hard coded for particular relations as of now.
		if ( ((in_state >> in_pos) & 0x1) != ((state >> out_pos) & 0x1) )
			printf("BAD: verification of pair of L_in failed  %d %d  state / state %016llx %016llx\n", 
				in_pos, out_pos, in_state, state);
	}
}

// This takes a pair from L_out that satisfies the outbound relation
// and determines the corresponding pair in L_in, by decryptiong
// for "inbound_rounds" rounds
void outbound::handle_pair(u64 state, u64 left, u64 right) {
	u64 RK[(inbound_rounds+1)][2];
	RK[inbound_rounds][0] = left;
	RK[inbound_rounds][1] = right;

	// backward keysched
	for (u64 r = inbound_rounds; r >= 1; --r) {
		left = RK[r][0];
		right = RK[r][1];

		inv_keysched_round(left, right, r);

		RK[r-1][0] = left;
		RK[r-1][1] = right;
	}
	for (int r = inbound_rounds; r >= 1; --r)
		decrypt_one_round(state, RK[(r-1)][0]);

	key_state in_pair;
	in_pair.left = RK[0][0];
	in_pair.right = RK[0][1];
	L_in.insert(std::make_pair(state, in_pair));
}

// Run the outbound phase
void outbound::run_phase() {
	std::set<std::pair<u64, key_state>>::iterator it;

#if DEBUG
	printf(">> Running outbound phase using %d = 2^{%f} starting points...\n", L_out.size(), log2((double)L_out.size()));
#endif	

	u64 in_state, state, lk, rk;
	for (it = L_out.begin(); it != L_out.end(); ++it) {
		// Set pair starting values
		in_state = (*it).first;
		state = in_state;
		lk = (*it).second.left;
		rk = (*it).second.right;

		// Encrypt
		for (int r = 0; r < num_rounds; ++r) {
			encrypt_one_round(state, lk);
			keysched_round(lk, rk, (r+inbound_rounds+1));
		}
		state ^= lk; // Final whitening key

		// Check outbound relation on single bit (positions provided as parameters)
		if (((in_state >> this->in_pos) & 0x1) == ((state >> this->out_pos) & 0x1)) {
			sat_pairs++;
#if VERIFY
			// Add the corresponding pair to L_in
			handle_pair(in_state, (*it).second.left, (*it).second.right);
#endif
		}
	}
#if DEBUG
	double correlation = 2.0f * (double) abs((L_out.size()/2)-sat_pairs) / L_out.size();
	printf(">> Finished. Found %d = 2^{%f} satisfying pairs, Pr(satisfied) = %f = 2^(%f)\n", 
		sat_pairs, log2((double)sat_pairs), 
		(double)sat_pairs / L_out.size(), log2((double)sat_pairs / L_out.size()));
	printf(">> Determined L_in with %d pairs\n\n", L_in.size());
	printf(">> This corresponds to a correlation of %f = 2^(%f)\n\n", correlation, log2(correlation));
#else
	printf("%d\n", sat_pairs);
#endif

#if VERIFY
	verify_L_in();
#endif
}

// This is for testing a particular weight-1 relation over some number of rounds,
// using e.g. fixed key or random keys, and random states.
void test_outbound_trail(int ip, int op, int rounds) {
	std::mt19937 rng;
	rng.seed(time(NULL));
	u64 in_state, state, lk, rk, s_lk, s_rk;
		
	std::uniform_int_distribution<u64> uni_dist(0, 0xFFFFFFFFFFFFFFFFull);
	
	for (int i = 0; i < 1000; ++i) {
		int _test_sat = 0;
		u64 runs = (0x1ull << 20);

		s_lk = uni_dist(rng);
		s_rk = uni_dist(rng);

		for (u64 i = 0; i < runs; ++i) {
			lk = s_lk;
			rk = s_rk;
			in_state = uni_dist(rng);
			state = in_state;

			// Encrypt
			for (u64 r = 0; r < rounds; ++r) {
				encrypt_one_round(state, lk);
				keysched_round(lk, rk, (r+7));
			}
			state ^= lk; // Final whitening key

			// Check relation
			if ( (((in_state >> ip) & 0x1) == ((state >> op) & 0x1)) )
				_test_sat++;
		}
		printf("%d\n", _test_sat);
	}
}

void test_hull(int inpos, int outpos, int rounds, int RUNS, int experiments) {
	// Seed RNG
	std::mt19937 rng;
	rng.seed(time(NULL));
	std::uniform_int_distribution<u64> uni_dist(0, 0xFFFFFFFFFFFFFFFFull);
	
	u64 m, in_state, out_state, k0_start, k1_start, k0, k1, count = 0;
	
	printf("# RUNS = %d\n\n", RUNS);
	
	// Make set of inputs
	int i, j, e;
	
	for (e = 0; e < experiments; ++e) {
		count = 0;
		m = uni_dist(rng);
	
		for (i = 0; i < RUNS; ++i) {
			in_state = m;
			k0 = uni_dist(rng);
			k1 = uni_dist(rng);
			#if PRESENT80
				k1_start &= 0xFFFFull;
			#endif
		
			out_state = in_state;
			
			// encrypt
			for (j = 0; j < rounds; ++j) {
				encrypt_one_round(out_state, k0);
				keysched_round(k0, k1, j+1);
			}
			
			// Check relation
			if ((((in_state >> inpos) & 0x1) == ((out_state >> outpos) & 0x1)))
				count++;
		}
		
		printf("%d\n", count);
	}
}


/**
* RANDOM NUMBER GENERATION
*/
#define __A 6364136223846793005ULL
#define __B 428856369ULL
unsigned long long __V = 8627181;
unsigned long long next_rand() {
	__V = __V*__A+__B;
	return __V;
}
#define random_u32() ((unsigned)((next_rand()>>32)))
#define random_u64() ((unsigned long long)random_u32()^next_rand())

using namespace std;

void kage() {	
	//~ std::mt19937 rng;
	//~ rng.seed(time(NULL));
	//~ std::uniform_int_distribution<u64> uni_dist(0, 0xFFFFFFFFFFFFFFFFull);
	
	int N = 1000;
	u64 S = 2619369;
	int ROUNDS = 6;
	
	double sqrtS = sqrt(S);
	cout << "S : " << S << "    sqrt(S) : " << sqrtS << endl;

	// Set up keys
	u64 k0vec[N], k1vec[N];
	u64 inner_counter[N];
	for (int i = 0; i < N; ++i) {
		inner_counter[i] = 0;
		
		k0vec[i] = random_u64();
		k1vec[i] = random_u64();
		
		#if PRESENT80
			k1 &= 0xFFFFull;
		#endif
	}
		//~ 
	//~ // Loop over inputs
	for (u64 k = 0; k < S; ++k) {
		// Random input
		u64 x = random_u64();
			
		// Loop over keys
		for (int i = 0; i < N; ++i) {
			// Re-set keys
			u64 ka = k0vec[i];
			u64 kb = k1vec[i];
			
			u64 y = x;
			
			// encrypt x to get y
			for (int j = 0; j < ROUNDS; ++j) {
				encrypt_one_round(y, ka);
				keysched_round(ka, kb, j+1);
			}
			
			// check relation
			if (((x >> 21) & 1) == ((y >> 21) & 1))
				inner_counter[i]++;
		}
	}
	for (int i = 0; i < N; ++i)
		cout << inner_counter[i] << ",";
	
	
	// Loop over keys k1...kN
	//~ for (int i = 0; i < N; ++i) {
		//~ u64 inner_counter = 0;
		//~ 
		//~ // Set random keys
		//~ const u64 k0 = random_u64();
		//~ const u64 k1 = random_u64();
		//~ 
		//~ for (u64 k = 0; k < S; ++k) {
			//~ // Re-set keys
			//~ u64 ka = k0;
			//~ u64 kb = k1;
			//~ 
			//~ // Random input
			//~ u64 x = random_u64();
			//~ u64 y = x;
			//~ 
			//~ // encrypt x to get y
			//~ for (int j = 0; j < ROUNDS; ++j) {
				//~ encrypt_one_round(y, ka);
				//~ keysched_round(ka, kb, j+1);
			//~ }
			//~ 
			//~ // check relation
			//~ if (((x >> 21) & 1) == ((y >> 21) & 1))
				//~ inner_counter++;
		//~ }
		//~ cout << inner_counter << endl;
	//~ }
}
