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

u64 sbox_one_bit[2][2][4][8] = {
{
	{ { 0x0ull, 0x2ull, 0x5ull, 0x6ull, 0x9ull, 0xBull, 0xCull, 0xFull }, // For bottom part
	{ 0x0ull, 0x1ull, 0x4ull, 0x5ull, 0x7ull, 0xBull, 0xCull, 0xEull },
	{ 0x3ull, 0x4ull, 0x5ull, 0x6ull, 0x8ull, 0xBull, 0xEull, 0xFull },
	{ 0x1ull, 0x2ull, 0x5ull, 0x8ull, 0xCull, 0xDull, 0xEull, 0xFull } },
	{ { 0x1ull, 0x3ull, 0x4ull, 0x7ull, 0x8ull, 0xAull, 0xDull, 0xEull },
	{ 0x2ull, 0x3ull, 0x6ull, 0x8ull, 0x9ull, 0xAull, 0xDull, 0xFull },
	{ 0x0ull, 0x1ull, 0x2ull, 0x7ull, 0x9ull, 0xAull, 0xCull, 0xDull },
	{ 0x0ull, 0x3ull, 0x4ull, 0x6ull, 0x7ull, 0x9ull, 0xAull, 0xBull } }
},
{
	{ { 0x1ull, 0x3ull, 0x4ull, 0x6ull, 0x9ull, 0xAull, 0xCull, 0xFull }, // For top part
	{ 0x0ull, 0x3ull, 0x4ull, 0x5ull, 0x7ull, 0x9ull, 0xCull, 0xEull },
	{ 0x3ull, 0x5ull, 0x6ull, 0x8ull, 0xBull, 0xCull, 0xEull, 0xFull },
	{ 0x0ull, 0x5ull, 0x6ull, 0x9ull, 0xAull, 0xBull, 0xCull, 0xDull } },
	{ { 0x0ull, 0x2ull, 0x5ull, 0x7ull, 0x8ull, 0xBull, 0xDull, 0xEull },
	{ 0x1ull, 0x2ull, 0x6ull, 0x8ull, 0xAull, 0xBull, 0xDull, 0xFull },
	{ 0x0ull, 0x1ull, 0x2ull, 0x4ull, 0x7ull, 0x9ull, 0xAull, 0xDull },
	{ 0x1ull, 0x2ull, 0x3ull, 0x4ull, 0x7ull, 0x8ull, 0xEull, 0xFull } }
} };

CombiIterator combi(16);

#define DEBUG 1
#define VERIFY 1

/**
	HELPER FUNCTIONS
**/
#if defined(_MSC_VER)
int __builtin_ffs(u64 x) {
	for (int p = 0; p < 4; ++p)
		if ((x >> p) & 0x1)
			return p+1;
	return 0;
}
#endif

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

// Considers each bit of an output mask, and marks the input to active S-boxes with 0xF in the input mask
u64 inbound::get_input_mask(u64 mask) {
	u64 ret = 0x0ull;
	// Iterate over all bits and mark the active S-boxes with input mask 0xF
	for (int i = 0; i < 64; ++i) {
		if ((mask >> i) & 0x1)
			ret |= (0xFull << (4*(i/4)));
	}
	return ret;
}

/**
	TESTING AND VALIDATION
**/

// Test if a value indeed works for following the trail down three levels AND up three levels
int inbound::test_value(u64 x, int top) {
	u64 y = x;

	if (top) {
		for (int r = 0; r < 3; ++r) {
			y = inv_perm_layer(y);
			y = inv_sub_layer(y);
		}
		if ( (maskIn[0] & y) && ((x >> trailPos[MATCH_ROUND]) & 0x1) )
			return 1;
		if ( ((maskIn[0] & y) == 0x0) && (((x >> trailPos[MATCH_ROUND]) & 0x1) == 0x0) )
			return 1;
	}
	else {
		for (int r = 0; r < 3; ++r) {
			y = sub_layer(y);
			y = perm_layer(y);
		}
		if ( ((x >> trailPos[MATCH_ROUND]) & 0x1) == ((y >> trailPos[NUM_ROUNDS-1]) & 0x1) )
			return 1;
	}

	return 0;
}

// Test found values
void inbound::test_level_values() {
#if DEBUG
	printf(">> Verifying that found values satisfy top/bottom 3 round trails...\n\n");
#endif
	std::set<u64>::iterator it;
	
	// Outputs from top
	for (it = W_top.begin(); it != W_top.end(); ++it)
		if (!test_value(*it, 1))
			printf(">> Value %016llX did not work for W_top\n", *it);

	// Inputs to bottom
	for (it = W_bottom.begin(); it != W_bottom.end(); ++it)
		if (!test_value(*it, 0))
			printf(">> Value %016llX did not work for W_bottom\n", *it);
}

int inbound::verify_K3_constraints(u64 l, u64 r) {
	u64 RK[NUM_ROUNDS][2];
	RK[MATCH_ROUND][0] = l;
	RK[MATCH_ROUND][1] = r;
	u64 left, right;

	// Generate the rounds keys
	left =  RK[MATCH_ROUND][0];
	right = RK[MATCH_ROUND][1];
	for (u64 i = 0; i < 3; ++i) {
		inv_keysched_round(left, right, (MATCH_ROUND-i));
		RK[MATCH_ROUND-i-1][0] = left;
		RK[MATCH_ROUND-i-1][1] = right;
	}

	left =  RK[MATCH_ROUND][0];
	right = RK[MATCH_ROUND][1];
	for (u64 i = 0; i < 2; ++i) {
		keysched_round(left, right, (MATCH_ROUND+1+i));
		RK[(MATCH_ROUND+1+i)][0] = left;
		RK[(MATCH_ROUND+1+i)][1] = right;
	}

	// Check that they are ok for the required bits in top/bottom
	for (int i = 0; i < NUM_ROUNDS; ++i) {
		if (keyMaskLeft[i] & RK[i][0]) {
			printf("%d : %016llX <=> %016llX\n", i, maskIn[i], RK[i][0]);
			return 0;
		}
	}
	return 1;
}

void inbound::verify_L_in() {
	u64 left, right, state;

	// Iterate over each pair in L_in
	for (std::set<std::pair<u64,key_state>>::iterator it = L_in.begin(); it != L_in.end(); ++it) {
		state = (*it).first;
		left = (*it).second.left;
		right = (*it).second.right;
		
		// Encrypt 6 rounds
		for (int r = 2; r < 8; ++r) {
			encrypt_one_round(state, left);			
			keysched_round(left, right, r);
		}

		// Check relation after 6 rounds
		if ( ((state >> trailPos[1]) & 0x1ull) != ((((*it).first) >> trailPos[6]) & 0x1ull) )
			printf("PROBLEM: trail was not verified for some element in L_in..\n");
	}
}

/**
   ATTACK-RELATED THINGS
**/

// Set up a specific trail for all levels. These masks denote the trail bit just BEFORE an S-layer (considered in an outward-in fashion)
void inbound::setup_trails() {
#if DEBUG
	printf(">> Setting up trails for %d levels...\n", NUM_ROUNDS);
#endif
	srand(time(NULL));
	for (int i = 0; i < NUM_ROUNDS; ++i)
		trailPos[i] = 21;

	for (int i = 0; i < NUM_ROUNDS; ++i)
		trailSbox[i] = trailPos[i] / 4;
#if DEBUG
	for (int i = 0; i < NUM_ROUNDS; ++i)
		printf(">> Round %d trail pos : %d, sbox : %d\n", i, trailPos[i], trailSbox[i]);
	printf("\n");
#endif
}

// Set up the input/output masks for the individual layers
void inbound::setup_masks() {
#if DEBUG
	printf(">> Setting up input/output masks for %d layers...\n", NUM_ROUNDS);
#endif
	maskIn[0] = (0x1ull << trailPos[0]);
	maskOut[0] = perm_layer(get_input_mask(maskIn[0]));
	for (int level = 1; level < 3; ++level) {
		maskIn[level] = maskOut[level-1];
		maskOut[level] = perm_layer(get_input_mask(maskIn[level]));
	}

	maskIn[0] = (0xFull << (trailSbox[0]*4)); // for m = 4 simultaneous trails

	maskOut[5] = (0x1ull << trailPos[5]);
	maskIn[5] = get_input_mask(maskOut[5]);
	for (int level = 4; level > 2; --level) {
		maskOut[level] = maskIn[level+1];
		maskIn[level] = get_input_mask(inv_perm_layer(maskOut[level]));
	}

	maskIn[3] = 0x0ull;
#if DEBUG
	for (int i = 0; i < NUM_ROUNDS; ++i)
		printf(">> Round %d mask in / out: %016llX / %016llX\n", i, maskIn[i], maskOut[i]);
	printf("\n");
#endif
}

void inbound::setup_fixed_key_bits() {
	for (int i = 0; i < NUM_ROUNDS; ++i)
		keyMaskLeft[i] = maskIn[i];

	//keyMaskLeft[0] |= 0x3FFFFFFFFF0FEE88ull;
	//keyMaskRight[0] = 0x73FFFFFFFFFFFFFFull;

#if DEBUG
	for (int i = 0; i < NUM_ROUNDS; ++i)
		printf(">> K_%d key mask : %016llX %016llX\n", i, keyMaskLeft[i], keyMaskRight[i]);
#endif
}

// Determine constraints on K_4
void inbound::determine_K3_constraints() {
#if DEBUG
	printf(">> Searching for the middle key 'locked bit positions' and their values...\n");
#endif
	// Make arrays for storing the round key masks and the values
	u64 kv[NUM_ROUNDS][2] = { 0 };
	u64 km[NUM_ROUNDS][2] = { 0 };
	u64 ml = 0x0ull;
	u64 mr = keyMaskRight[0];
	u64 vl = 0x0ull;
	u64 vr = 0x0ull;
	
	int p;
	for (u64 round = 0; round <= MATCH_ROUND; ++round) {
		// Check nibble 15
		if ((ml >> 60) & 0xF) {
			if (((keyMaskLeft[round] >> 60) & 0xF) == 0x0ull) {				
				p = __builtin_ffs((ml >> 60) & 0xF) - 1;
				
				ml |= (0xFull << 60);
				vl = (sbox_one_bit[1][0][p][0] << 60) ^ (vl & 0x0FFFFFFFFFFFFFFFull);
			}
			else
				printf("BAD\n");
		}
		// Check nibble 14
		if ((ml >> 56) & 0xF) {
			if (((keyMaskLeft[round] >> 56) & 0xF) == 0x0ull) {				
				p = __builtin_ffs((ml >> 56) & 0xF) - 1;				
				
				ml |= (0xFull << 56);
				vl = (sbox_one_bit[1][0][p][0] << 56) ^ (vl & 0xF0FFFFFFFFFFFFFFull);
			}
			else
				printf("BAD\n");
		}
		
		// Mark fixed bits for this level
		ml |= keyMaskLeft[round];
		
		// Modify value vectors to account for round constant
		vl ^= (round >> 2);
		vr ^= (round << 62);
		
		// Store for level
		km[round][0] |= ml;
		km[round][1] |= mr;
		kv[round][0] |= vl;
		kv[round][1] |= vr;
#if DEBUG
		printf("T # km[%d] : %016llX %016llX   wt : %d + %d\n", round, km[round][0], km[round][1], weight(km[round][0]), weight(km[round][1]));
#endif
		// Rotate for next level
		rotate_keys(ml, mr);
		rotate_keys(vl, vr);
	}
		
	// Bottom
	ml = 0x0ull;
	mr = 0x0ull;
	vl = 0x0ull;
	vr = 0x0ull;
	
	for (u64 round = (NUM_ROUNDS-1); round >= MATCH_ROUND; --round) {
		// Check nibble 15
		if ((keyMaskLeft[round] >> 60) & 0xF) {
			if (((ml >> 60) & 0xF) == 0x0ull) {				
				p = __builtin_ffs((keyMaskLeft[round] >> 60) & 0xF) - 1;
				
				ml |= (0xFull << 60);
				vl = (sbox_one_bit[1][0][p][0] << 60) ^ (vl & 0x0FFFFFFFFFFFFFFFull);
			}
			else
				printf("BAD\n");
		}
		// Check nibble 14
		if ((keyMaskLeft[round] >> 56) & 0xF) {
			if (((ml >> 56) & 0xF) == 0x0ull) {				
				p = __builtin_ffs((keyMaskLeft[round] >> 56) & 0xF) - 1;
				
				ml |= (0xFull << 56);
				vl = (sbox_one_bit[1][0][p][0] << 56) ^ (vl & 0xF0FFFFFFFFFFFFFFull);
			}
			else
				printf("BAD\n");
		}
		
		// Mark fixed bits for this level
		ml |= keyMaskLeft[round];
	
		// Modify value vectors to account for round constant
		vl ^= (round >> 2);
		vr ^= (round << 62);
	
		// Store for level
		km[round][0] |= ml;
		km[round][1] |= mr;
		kv[round][0] |= vl;
		kv[round][1] |= vr;
#if DEBUG
		printf("B # km[%d] : %016llX %016llX   wt : %d + %d\n", round, km[round][0], km[round][1], weight(km[round][0]), weight(km[round][1]));
#endif
		// Rotate for next level
		inv_rotate_keys(ml, mr);
		inv_rotate_keys(vl, vr);
	}
	
	// Set the midkey mask and value for the use of other methods
	K3_mask[0] = km[MATCH_ROUND][0];
	K3_mask[1] = km[MATCH_ROUND][1];
	K3_val[0] = kv[MATCH_ROUND][0] & K3_mask[0];
	K3_val[1] = kv[MATCH_ROUND][1] & K3_mask[1];
	
#if DEBUG
	printf(">> Mid key mask            : %016llX %016llX   wt : %d + %d\n", K3_mask[0], K3_mask[1], weight(K3_mask[0]), weight(K3_mask[1]));
	printf(">> AND'ed with mid key val : %016llX %016llX\n\n", K3_val[0], K3_val[1]);
#endif
#if VERIFY
	printf(">> Testing mid key. Is it good? : %d\n", verify_K3_constraints(K3_val[0], K3_val[1]));
	int idx;
	u64 l = K3_val[0];
	u64 r = K3_val[1];
	printf(">> Doing more tests : Checking if flipping any non-locked bit is ok. An error will be printed if it failed...\n\n");
	for (int t = 0; t < 128; ++t) {
		idx = t;
		if (idx < 64) {
			// Check if bit t of right key mask is marked as non-fixed
			if (((K3_mask[1] >> idx) & 0x1) == 0) {
				if (!verify_K3_constraints(l, r ^ (0x1ull << idx)))
					printf("BAD ! idx = %d\n", idx);
			}
		}
		else {
			idx -= 64;
			// Check if bit t of left key mask is marked as non-fixed
			if (((K3_mask[0] >> idx) & 0x1) == 0) {
				if (!verify_K3_constraints(l ^ (0x1ull << idx), r))
					printf("BAD ! idx = %d\n", idx + 64);
			}
		}
	}
#endif
}

// Construct inputs resp. outputs for a particular level
void inbound::do_layer(int level, u64 x, u64 mask) {
	std::vector<u64> sbox_combi[16];	// Hold the possible outputs of S-box j
	int s, p;
	int target_set = (level < MATCH_ROUND) ? 1 : 0;
	
	// Iterate over S-boxes
	for (s = 0; s < 16; ++s) {
		// Check if S-box is even active
		if ((mask >> (4*s)) & 0xF) {
			// Iterate over the bits that are output of this S-box and see if they appear in the mask
			p = __builtin_ffs((mask >> (4*s)) & 0xF) - 1;	

			for (int i = 0; i < 8; ++i)
				sbox_combi[s].push_back(sbox_one_bit[target_set][(x >> (4*s + p)) & 0x1][p][i] << (s*4));
		}
		else // insert 0 for dummy
			sbox_combi[s].push_back(0x0ull);
	}

	// Pruning phase for keeping only those inputs to the trail S-box which satisfy the trail
	std::vector<u64> pruned;	
	if (level < 4) {
		for (unsigned int i = 0; i < sbox_combi[trailSbox[level]].size(); ++i)
			if ((sbox_combi[trailSbox[level]][i] & (0x1ull << trailPos[level])) == (inv_sub_layer(sbox_combi[trailSbox[level]][i]) & (0x1ull << trailPos[level])))
				pruned.push_back(sbox_combi[trailSbox[level]][i]);
	}
	else {
		for (unsigned int i = 0; i < sbox_combi[trailSbox[level]].size(); ++i)
			if ((sbox_combi[trailSbox[level]][i] & (0x1ull << trailPos[level])) == (sub_layer(sbox_combi[trailSbox[level]][i]) & (0x1ull << trailPos[level])))
				pruned.push_back(sbox_combi[trailSbox[level]][i]);
	}
	sbox_combi[trailSbox[level]] = pruned;
	
	// Combining phase. Create combination iterator and set limits
	combi.reset();
	for (int i = 0; i < 16; ++i)
		combi.limit[i] = sbox_combi[i].size();

	u64 y;
	int N = 0;
	while (combi.has_next() && N < 100) {//(level < (MATCH_ROUND-1) || level > MATCH_ROUND || N < 15)) {
		y = 0;
		for (s = 0; s < 16; ++s)
			y ^= sbox_combi[s][combi.combination[s]];

		if (level < MATCH_ROUND)
			layerOut[level].insert(perm_layer(y));
		else
			layerIn[level].insert(y);

		N++;
	}
}

void inbound::handle_middle_match(u64 x, u64 _k4left, u64 _k4right) {
	// First, do 3 rounds of inverse key scheduling
	u64 left = _k4left;
	u64 right = _k4right;

	// forward keysched
	for (u64 r = (MATCH_ROUND+1); r <= NUM_ROUNDS; ++r) {
		encrypt_one_round(x, left);
		keysched_round(left, right, r);

#if VERIFY
		if ( (r < NUM_ROUNDS) && (maskIn[r] & left) )
			printf("Bad RK[%d] : %016llX    maskIn[%d] : %016llX\n", r, left, r, maskIn[r]);
#endif
	}
	
	// Add to L_out
	key_state out_pair;
	out_pair.left = left;
	out_pair.right = right;
	L_out.insert(std::make_pair(x, out_pair));
}

// Find the intersection between the outputs of level 3 (first apply the permutation layer!)
// and the input to layer 2
void inbound::match_in_the_middle() {
	L_out.clear();

	// Create two sets for the top resp. bottom part, where values are ANDED by the K3_mask
	std::set<u64>::iterator it;
	std::set<u64> list_top, list_bottom;
	std::map<u64, std::vector<u64> > map_top, map_bottom;
	u64 match_size = 0;
	u64 t, b, fully_det_k3_left, k3_right;
#if DEBUG
	printf(">> Trying to match W_top with W_bottom...\n");
#endif	
	// Make a set of W_bottom where the K4^l value has been added
	std::set<u64> W_bottom_translated;
	for (it = W_bottom.begin(); it != W_bottom.end(); ++it)
		W_bottom_translated.insert((*it) ^ K3_val[0]);

	for (it = W_top.begin(); it != W_top.end(); ++it) {
		list_top.insert((*it) & K3_mask[0]);
		map_top[((*it) & K3_mask[0])].push_back(*it);
	}
	for (it = W_bottom_translated.begin(); it != W_bottom_translated.end(); ++it) {
		list_bottom.insert((*it) & K3_mask[0]);
		map_bottom[((*it) & K3_mask[0])].push_back(*it);
	}
	
	std::set<u64> intersection;
	u64 x;
	set_intersection(list_top.begin(), list_top.end(), list_bottom.begin(), list_bottom.end(), inserter(intersection, intersection.begin()));
#if DEBUG
	printf(">> Computed intersection. Found %d elements\n", intersection.size());
#endif

	for (it = intersection.begin(); it != intersection.end(); ++it) {
		random_shuffle(map_top[*it].begin(), map_top[*it].end());
		random_shuffle(map_bottom[*it].begin(), map_bottom[*it].end());

		for (std::vector<u64>::iterator b_iter = map_bottom[(*it)].begin(); b_iter != map_bottom[(*it)].end(); ++b_iter) {
			for (std::vector<u64>::iterator t_iter = map_top[(*it)].begin(); t_iter != map_top[(*it)].end(); ++t_iter) {
				t = *t_iter;
				b = *b_iter ^ K3_val[0];
				
#if VERIFY
				// Check that t is indeed in W_top and b + K_4^l is indeed in W_bottom
				if (W_top.count(t) == 0)
					printf("THIS IS BAD:  t was not in W_top\n");
				if (W_bottom.count(b) == 0)
					printf("THIS IS BAD:  b was not in W_bottom\n");
#endif

				if ( (t & (0x1ull << trailPos[MATCH_ROUND])) == (b & (0x1ull << trailPos[MATCH_ROUND])) ) {
					fully_det_k3_left = t ^ b;
					
#if VERIFY
					if ((fully_det_k3_left & K3_mask[0]) != K3_val[0])
						printf("THIS IS BAD: key mismatch\n");
#endif
					x = rand();
					while ( (x != std::numeric_limits<u64>::max()) && (match_size < L_OUT_SIZE) ) {
						if ((x & K3_mask[1]) == K3_val[1]) {
							match_size++;
							handle_middle_match(t, fully_det_k3_left, x);	
							x++;
						}
						else
							x += (x & K3_mask[1]) ^ K3_val[1];
					}
				}
			}
		}
	}
#if DEBUG
	printf(">> Matching in the middle size : %d = 2^{%f}\n>> Size of L_out : %d = 2^{%f}\n\n", match_size, log2((double)match_size), 
		L_out.size(), log2((double)L_out.size()));
#endif
}

// Bottom 3 rounds
void inbound::determine_W_bottom() {	
	std::set<u64>::iterator it;
	do_layer((NUM_ROUNDS-1), 0x0ull << trailPos[(NUM_ROUNDS-1)], maskOut[(NUM_ROUNDS-1)]);
	do_layer((NUM_ROUNDS-1), 0x1ull << trailPos[(NUM_ROUNDS-1)], maskOut[(NUM_ROUNDS-1)]);
#if DEBUG
	printf(">> %d possible inputs to level %d found\n", layerIn[(NUM_ROUNDS-1)].size(), (NUM_ROUNDS-1));
#endif
	for (int level = (NUM_ROUNDS-2); level > (MATCH_ROUND-1); --level) {
		for (it = layerIn[level+1].begin(); it != layerIn[level+1].end(); ++it)
			do_layer(level, inv_perm_layer(*it), inv_perm_layer(maskOut[level]));
#if DEBUG
		printf(">> %d possible inputs to level %d found\n", layerIn[level].size(), level);
#endif
	}
	W_bottom = layerIn[MATCH_ROUND];
#if DEBUG
	printf(">> #W_bottom = %d = 2^{%f}\n\n", W_bottom.size(), log2((double)W_bottom.size()));
#endif
}

// Top 3 rounds
void inbound::determine_W_top() {
	std::set<u64>::iterator it;
	layerOut[0].insert( perm_layer( (0x2ull << (trailSbox[0]*4)) ) );
	layerOut[0].insert( perm_layer( (0xCull << (trailSbox[0]*4)) ) );
#if DEBUG
	printf(">> %d possible outputs from level %d found\n", layerOut[0].size(), 1);
#endif
	for (int level = 1; level < 3; ++level) {
		for (it = layerOut[level-1].begin(); it != layerOut[level-1].end(); ++it)
			do_layer(level, *it, maskIn[level]);
#if DEBUG
		printf(">> %d possible outputs from level %d found\n", layerOut[level].size(), level);
#endif
	}
	W_top = layerOut[(MATCH_ROUND-1)];
#if DEBUG
	printf(">> #W_top = %d = 2^{%f}\n\n", W_top.size(), log2((double)W_top.size()));
#endif
}

// Wrapper for the search
void inbound::run_phase() {
	setup_trails();
	setup_masks();
	setup_fixed_key_bits();
	determine_K3_constraints();

	determine_W_top();
	determine_W_bottom();
	
#if VERIFY 
	test_level_values();
#endif
}