#include "present.h"
#include "inbound.h"
#include "combi.h"

#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <algorithm>
#include <vector>
#include <map>
#include <limits>
#include <math.h>

#define DEBUG 1
#define VERIFY 1

// This table corresponds to a generalised version of the T tables
// described in Appendix. More specifically, sbox_one_bit[a][b] holds
// inputs to the PRESENT S-Box s.t. the b-th least significant
u64 sbox_one_bit[2][4][8] =
{
	{ { 0x0ull, 0x2ull, 0x5ull, 0x6ull, 0x9ull, 0xBull, 0xCull, 0xFull }, // For bottom part
	{ 0x0ull, 0x1ull, 0x4ull, 0x5ull, 0x7ull, 0xBull, 0xCull, 0xEull },
	{ 0x3ull, 0x4ull, 0x5ull, 0x6ull, 0x8ull, 0xBull, 0xEull, 0xFull },
	{ 0x1ull, 0x2ull, 0x5ull, 0x8ull, 0xCull, 0xDull, 0xEull, 0xFull } },
	{ { 0x1ull, 0x3ull, 0x4ull, 0x7ull, 0x8ull, 0xAull, 0xDull, 0xEull },
	{ 0x2ull, 0x3ull, 0x6ull, 0x8ull, 0x9ull, 0xAull, 0xDull, 0xFull },
	{ 0x0ull, 0x1ull, 0x2ull, 0x7ull, 0x9ull, 0xAull, 0xCull, 0xDull },
	{ 0x0ull, 0x3ull, 0x4ull, 0x6ull, 0x7ull, 0x9ull, 0xAull, 0xBull } }
};

// Object used for genrating direct products
CombiIterator combi(16);

inbound::inbound() {
	srand(time(NULL));
	maskIn = (u64 *)calloc(NUM_ROUNDS, sizeof(u64));
	maskOut = (u64 *)calloc(NUM_ROUNDS, sizeof(u64));
	keyMaskLeft = (u64 *)calloc(32, sizeof(u64));
	keyMaskRight = (u64 *)calloc(32, sizeof(u64));
	trailSbox = (int *)calloc(NUM_ROUNDS, sizeof(int));
	trailPos = (int *)calloc(NUM_ROUNDS, sizeof(int));
}

inbound::~inbound() {
	free(maskIn);
	free(maskOut);
	free(keyMaskLeft);
	free(keyMaskRight);
	free(trailSbox);
	free(trailPos);
}

///////////////////////////////////////////////////////
//// HELPER FUNCTIONS		 					   ////
///////////////////////////////////////////////////////

// For finding the finding first set bit (ffs)
#if defined(_MSC_VER)
int __builtin_ffs(u64 x) {
	for (int p = 0; p < 4; ++p)
		if ((x >> p) & 0x1)
			return p+1;
	return 0;
}
#endif

// Considers each bit of an output mask, and marks 
//the input to active S-boxes with 0xF in the input mask
u64 inbound::get_input_mask(u64 mask) {
	u64 ret = 0x0ull;
	// Iterate over all bits and mark the active S-boxes with input mask 0xF
	for (int i = 0; i < 64; ++i) {
		if ((mask >> i) & 0x1)
			ret |= (0xFull << (4*(i/4)));
	}
	return ret;
}

///////////////////////////////////////////////////////
//// TESTING AND VALIDATION 					   ////
///////////////////////////////////////////////////////

// Test found set W
int inbound::verify_W_set() {
	std::set<u64>::iterator it;
	u64 x, kl, kr;
	// Inputs to E_in
	for (it = W.begin(); it != W.end(); ++it) {
		x = *it;
		kl = k0_val[0];
		kr = k0_val[0];
		x ^= kl;

		for (u64 round = 0x1; round <= NUM_ROUNDS; ++round) {
			encrypt_one_round(x, kl);
			
			if ( (((*it) >> trailPos[0]) & 0x1) != ((x >> trailPos[round-1]) & 0x1) ) {
				printf("problem after %d round\n", round);
				printf("*it = %016llx\n", *it);
				printf("x   = %016llx\n", x);
				return 0;
			}

			keysched_round(kl, kr, round);
		}
	}
	return 1;
}

void inbound::verify_W_layer(int round) {
	std::set<u64>::iterator it;
	u64 x;
	for (it = layerIn[round].begin(); it != layerIn[round].end(); ++it) {
		x = *it;
		x = sub_layer(x);
		x = perm_layer(x);

		if ((((*it) >> 21) & 0x1) != ((x >> 21) & 0x1))
			printf("fail!!! inbit = %x   outbit = %x\n",
				(((*it) >> 21) & 0x1), ((x >> 21) & 0x1));
	
	}
}

int inbound::verify_k0_constraints(u64 l, u64 r) {
	u64 RK[NUM_ROUNDS][2];
	RK[0][0] = l;
	RK[0][1] = r;
	u64 left, right;

	// Generate the rounds keys
	left =  RK[0][0];
	right = RK[0][1];
	for (u64 i = 0; i < NUM_ROUNDS-1; ++i) {
		keysched_round(left, right, (i+1));
		RK[i+1][0] = left;
		RK[i+1][1] = right;
	}

	// Check that they are ok for the required bits in top/bottom
	for (int i = 0; i < NUM_ROUNDS; ++i) {
		if ((keyMaskLeft[i] & RK[i][0]) != (keyMaskLeft[i] & round_keys[i][0]) ) {
			printf("%d : %016llX <=> %016llX\n", i, keyMaskLeft[i], RK[i][0]);
			return 0;
		}
	}
	return 1;
}

int inbound::check_k0_mask_val(u64 ml, u64 mr, u64 vl, u64 vr) {
	u64 l = ml & vl;
	u64 r = mr & vr;

	#if DEBUG
	printf(">> Testing matching key. Is it good? : %d\n", verify_k0_constraints(l, r));
	printf(">> Doing more tests : Checking if flipping any non-locked bit is ok. An error will be printed if it failed...\n\n");
	#endif
	int idx;
	int c1 = PRESENT80 ? 80 : 128;
	int c2 = PRESENT80 ? 16 : 64;

	for (int t = 0; t < c1; ++t) {
		idx = t;
		if (idx < c2) {
			// Check if bit t of right key mask is marked as non-fixed
			if (((mr >> idx) & 0x1) == 0) {
				if (!verify_k0_constraints(l, r ^ (0x1ull << idx)))
					printf("BAD ! idx = %d\n", idx);
			}
		}
		else {
			idx -= c2;
			// Check if bit t of left key mask is marked as non-fixed
			if (((ml >> idx) & 0x1) == 0) {
				if (!verify_k0_constraints(l ^ (0x1ull << idx), r))
					printf("BAD ! idx = %d\n", idx + c2);
			}
		}
	}
}

///////////////////////////////////////////////////////
//// INITIALIZATION FUNCTIONS 					   ////
///////////////////////////////////////////////////////

// Set up a specific trail for all levels. These masks denote the trail bit just 
// BEFORE an S-layer (considered in an outward-in fashion)
void inbound::setup_trails() {
	#if DEBUG
	printf(">> Setting up trails for %d levels...\n", NUM_ROUNDS);
	#endif
	for (int i = 0; i < NUM_ROUNDS; ++i) {
		trailPos[i] = 21;
		trailSbox[i] = trailPos[i] / 4;
		#if DEBUG
		printf(">> Round %d trail pos : %d, sbox : %d\n", i, trailPos[i], trailSbox[i]);
		#endif
	}
	#if DEBUG
	printf("\n");
	#endif
}

// Set up the input/output masks for the individual rounds
void inbound::setup_masks() {
	#if DEBUG
	printf(">> Setting up input/output masks for %d rounds...\n", NUM_ROUNDS);
	#endif

	maskOut[NUM_ROUNDS-1] = 0x1ull << trailPos[NUM_ROUNDS-1];
	maskIn[NUM_ROUNDS-1] = get_input_mask(maskOut[NUM_ROUNDS-1]);
	for (int i = NUM_ROUNDS-2; i >= 0; --i) {
		maskOut[i] = maskIn[i+1];
		maskIn[i] = get_input_mask(inv_perm_layer(maskOut[i]));
	}
	maskIn[0] = 0x1ull << trailPos[0];

#if DEBUG
	for (int i = 0; i < NUM_ROUNDS; ++i)
		printf(">> Round %d mask in / out: %016llX / %016llX\n", i, maskIn[i], maskOut[i]);
	printf("\n");
#endif
}

// Store masks for the key bits that are fixed in each round
void inbound::setup_fixed_key_bits() {
	for (int i = 0; i < NUM_ROUNDS; ++i) {
		keyMaskLeft[i] = maskIn[i];
		#if DEBUG
		printf(">> k_%d key mask : %016llX %016llX\n", i, keyMaskLeft[i], keyMaskRight[i]);
		#endif
	}
}

///////////////////////////////////////////////////////
//// FOR CONSTRUCTING W AND L_OUT				   ////
///////////////////////////////////////////////////////

// Construct inputs resp. outputs for a particular level
void inbound::do_layer(int level, u64 x, u64 mask) {
	std::vector<u64> sbox_combi[16];	// Hold the possible outputs of S-box j
	int s, p;

	// Iterate over S-boxes
	for (s = 0; s < 16; ++s) {
		// Check if S-box is even active
		if ((mask >> (4*s)) & 0xF) {
			// Iterate over the bits that are output of this S-box and see if they appear in the mask
			p = __builtin_ffs((mask >> (4*s)) & 0xF) - 1;	

			for (int i = 0; i < 8; ++i)
				sbox_combi[s].push_back(sbox_one_bit[(x >> (4*s + p)) & 0x1][p][i] << (s*4));
		}
		else // insert 0 for dummy
			sbox_combi[s].push_back(0x0ull);
	}

	// Pruning phase for keeping only those inputs to the trail S-box which satisfy the trail
	std::vector<u64> pruned;
	for (unsigned int i = 0; i < sbox_combi[trailSbox[level]].size(); ++i)
		if ((sbox_combi[trailSbox[level]][i] & (0x1ull << trailPos[level])) == (sub_layer(sbox_combi[trailSbox[level]][i]) & (0x1ull << trailPos[level])))
			pruned.push_back(sbox_combi[trailSbox[level]][i]);
	sbox_combi[trailSbox[level]] = pruned;
	
	// Combining phase. Create combination iterator and set limits
	combi.reset();
	for (int i = 0; i < 16; ++i)
		combi.limit[i] = sbox_combi[i].size();

	u64 y;
	int N = 0;
	while (combi.has_next() && (level > 0 || N < 50)) {
		y = 0;
		for (s = 0; s < 16; ++s)
			y ^= sbox_combi[s][combi.combination[s]];

		layerIn[level].insert(y);
		N++;
	}
}

void inbound::offset_W_layer_byKey(int round) {
	std::set<u64> tmpSet;
	std::set<u64>::iterator it;
	for (it = layerIn[round].begin(); it != layerIn[round].end(); ++it) {
		tmpSet.insert( (*it) ^ round_keys[round][0] );
	}
	layerIn[round] = tmpSet;
}

// Determine the set W
void inbound::determine_W() {
	std::set<u64>::iterator it;
	do_layer((NUM_ROUNDS-1), 0x0ull << trailPos[(NUM_ROUNDS-1)], maskOut[(NUM_ROUNDS-1)]);
	do_layer((NUM_ROUNDS-1), 0x1ull << trailPos[(NUM_ROUNDS-1)], maskOut[(NUM_ROUNDS-1)]);
	
	verify_W_layer(NUM_ROUNDS-1);
	offset_W_layer_byKey(NUM_ROUNDS-1);

	#if DEBUG
	printf(">> %d possible inputs to level %d found\n", layerIn[(NUM_ROUNDS-1)].size(), (NUM_ROUNDS-1));
	#endif
	for (int level = (NUM_ROUNDS-2); level >= 0; --level) {
		for (it = layerIn[level+1].begin(); it != layerIn[level+1].end(); ++it)
			do_layer(level, inv_perm_layer(*it), inv_perm_layer(maskOut[level]));

		#if DEBUG
		printf(">> %d possible inputs to level %d found\n", layerIn[level].size(), level);
		#endif

		verify_W_layer(level);
		if (level > 0)
			offset_W_layer_byKey(level);
	}
	W = layerIn[0];

	#if DEBUG
	printf(">> #W = %d = 2^{%f}\n\n", W.size(), log2((double)W.size()));
	#endif
}

void inbound::construct_L_out() {
	std::mt19937 rng;
	rng.seed(time(NULL));
	std::uniform_int_distribution<u64> uni_dist(0, 0xFFFFFFFFFFFFFFFFull);

	u64 randKL, randKR, left, right, x, y;
	std::set<u64>::iterator it;
	
	while (L_out.size() < L_OUT_SIZE) {
		randKL = uni_dist(rng) & (k0_mask[0] ^ 0xFFFFFFFFFFFFFFFFull);
		randKR = uni_dist(rng) & (k0_mask[1] ^ 0xFFFFFFFFFFFFFFFFull);

		if (PRESENT80)
			randKR &= 0xFFFFull;

		it = W.begin();
		while (it != W.end() && L_out.size() < L_OUT_SIZE) {
			x = *it;
			left = randKL;
			right = randKR;
			x ^= left;
					
			// forward keysched
			for (u64 r = 0x1; r <= NUM_ROUNDS; ++r) {
				y = x;
				encrypt_one_round(x, left);
				if ( ((y >> 21) & 0x1) != ((x >> 21) & 0x1) )
						printf("mistake in verification of inbound trail\n");
				keysched_round(left, right, r);
			}

			// Add to L_out
			key_state out_pair;
			out_pair.left = left;
			out_pair.right = right;
			L_out.insert(std::make_pair(x, out_pair));

			it++;
		}		
	}
}

///////////////////////////////////////////////////////
//// CODE RELATED TO DETERMINING CONSTRAINTS ON k3 ////
///////////////////////////////////////////////////////

void inbound::inv_ks_helper(u64 &ml, u64 &mr, u64 &vl, u64 &vr, u64 rc) {
#if DEBUG
	printf(" ks helper in: %016llx %016llx\n", ml, mr);
	printf("               %016llx %016llx\n", vl, vr);
#endif

	// add RC
	#if PRESENT80
	vl ^= (rc >> 1);
	vr ^= ((rc & 0x1) << 15);
	#else
	vl ^= (rc >> 2);
	vr ^= ((rc & 0x3) << 62);
	#endif

#if DEBUG
	printf(" after key add %016llx %016llx\n", ml, mr);
	printf("               %016llx %016llx\n", vl, vr);
#endif

	// check for active s-box
	if (ml & 0xF000000000000000ull) {
		printf("sbox active in keysched generating k_%d\nmask %016llx %04llx\nval  %016llx %04llx\n\n", 
			(rc-1), ml, mr, vl, vr);
	
		vl = (inv_sub_layer(ml & vl) & 0xF000000000000000ull) ^ (vl & 0x0FFFFFFFFFFFFFFFull);
		ml |= 0xF000000000000000ull;
	}
	#if PRESENT80 == 0
	// check for active s-box
	if (ml & 0x0F00000000000000ull) {
		printf("sbox active in keysched generating k_%d\nmask %016llx %04llx\nval  %016llx %04llx\n\n", 
			(rc-1), ml, mr, vl, vr);
	
		vl = (sub_layer(ml & vl) & 0x0F00000000000000ull) ^ (vl & 0xF0FFFFFFFFFFFFFFull);
		ml |= 0x0F00000000000000ull;
	}
	#endif

	// rotate
	inv_rotate_keys(ml, mr);
	inv_rotate_keys(vl, vr);

#if DEBUG
	printf(" after rotate  %016llx %016llx\n", ml, mr);
	printf("               %016llx %016llx\n", vl, vr);
#endif
}

void inbound::determine_k0_constraints() {
	u64 ml = 0; u64 mr = 0;
	u64 vl = 0; u64 vr = 0;
	u64 km[NUM_ROUNDS][2];
	u64 kv[NUM_ROUNDS][2];

#if DEBUG
	printf("\n>> Determining constraints on k_0 in inbound phase...\n");
#endif
	for (u64 round = NUM_ROUNDS; round > 0; --round) {
		if (ml & keyMaskLeft[round-1]) {
			printf("be careful!\n");
			//> ml        %016llx\n> ml&vl     %016llx\n> keyMaskLeft       %016llx\n> ml&vl&keyMaskLeft %016llx\n", 
			//	ml, (ml&vl), keyMaskLeft[round-1], (ml&vl&keyMaskLeft[round-1]));
		}
		ml |= keyMaskLeft[round-1];

		#if DEBUG
			printf("req on k_%d\nmask %016llx %016llx\nval  %016llx %016llx\n\n", round-1, ml, mr, vl, vr);
		#endif

		km[round-1][0] = ml;	km[round-1][1] = mr;
		kv[round-1][0] = vl;	kv[round-1][1] = vr;
		round_keys[round-1][0] = km[round-1][0] & kv[round-1][0];
		round_keys[round-1][1] = km[round-1][1] & kv[round-1][1];

		inv_ks_helper(ml, mr, vl, vr, round-1);
	}

	k0_mask[0] = km[0][0];	k0_mask[1] = km[0][1];
	k0_val[0]  = km[0][0] & kv[0][0];
	k0_val[1]  = km[0][1] & kv[0][1];

	check_k0_mask_val(km[0][0], km[0][1], kv[0][0], kv[0][1]);
}

// Wrapper for the search
void inbound::run_phase() {
	setup_trails();
	setup_masks();
	setup_fixed_key_bits();
	determine_k0_constraints();
	determine_W();
	//determine_W2();
	
#if VERIFY
	printf(">> Testing the W set... is it good ? : %d\n\n", verify_W_set());
#endif

	construct_L_out();
}