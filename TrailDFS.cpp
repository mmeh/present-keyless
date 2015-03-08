#include "TrailDFS.h"

int TrailDFS::hamming_weight(u64 x) {
	int r = 0;
	while (x) {
		if (x & 0x1)
			r++;
		x >>= 1;
	}
	return r;
}

void TrailDFS::print_trail() {
	printf("TRAIL : ");
	for (vector<int>::iterator it = trail.begin(); it != trail.end(); ++it)
		printf("%d ", *it);
	printf("\n");
}

void TrailDFS::test_in_mask(u64 left, u64 right) {
	u64 t0, t1;

	int round = 0;
	for (round = 0; round < 5; ++round) {
		printf("r : %d  %016llX %016llX  wt : %d + %d\n", round, left, right, hamming_weight(left), hamming_weight(right));

		// rotate
		t0 = left;
		t1 = right;
		right = (t0 >> 3) ^ ((t1 & 0x7) << 61);
		left = (t1 >> 3) ^ ((t0 & 0x7) << 61);

		// Dirty 3 round counter bits
		left |= 0x1ull;
		right |= (0x3ull << 62);

		// Dirty 8 S-box bits
		left |= (0xFFull << 56);
	}

	printf("r : %d  %016llX %016llX  wt : %d + %d\n", round, left, right, hamming_weight(left), hamming_weight(right));
}

/**
 64-bit input mask "kmask" dictates by a 1 where a bit position is already under control,
 a 0 means it is not under control.
 Given "kmask", we compute the modified kmask which determines which bits are controlled
 in the key mask obtained by an inverse key schedule round. If the new key mask allows us 
 to control the desired key bits for that particular level.

 As we will stay below 7 rounds, the counter added in the key-schedule can be assumed at most 3 key bits.
*/
int TrailDFS::key_search(int level, int trail_sbox, u64 in_left_mask, u64 in_right_mask) {
	// Check if we have reached level 6 -- if so, we are done.
	if (level == 5) {
		printf(">> Level %d reached    wt : %d + %d    ", level, hamming_weight(in_left_mask), hamming_weight(in_right_mask));
		print_trail();
		//test_inkeymask(in_left_mask, in_right_mask);
		return 1;
	}

	// Temp. variables
	u64 my_left_mask, my_right_mask, nxt_left_mask, nxt_right_mask, t0, t1;
	int sbox, num_sbox_cands;

	// This is a hack. 
	// If trail_sbox == 16 then we allow to check on all possible trail
	// S-boxes above (16 of them), but if not, there are just 3 candidates, c.f. the array "trail_sbox_cand"
	// in the header file
	num_sbox_cands = (trail_sbox == 16) ? 16 : 3;

	// Iterate over possible trail S-boxes for the layer above
	// Depending on this, is the four key bits for this level, which we need to control
	for (int i = 0; i < num_sbox_cands; ++i) {
		sbox = trail_sbox_cand[trail_sbox][i];

		// Copy to new variables
		my_left_mask = in_left_mask;
		my_right_mask = in_right_mask;

		// Check against the 4 required key bits for that particular trail S-box
		if ((my_left_mask & levelMask[sbox]) == 0x0ull) {
			my_left_mask |= levelMask[sbox];
			trail.push_back(sbox);
			
			nxt_left_mask = my_left_mask;
			nxt_right_mask = my_right_mask;

			// Dirty 3 round counter bits
			nxt_left_mask |= 0x1ull;
			nxt_right_mask |= (0x3ull << 62);

			// Dirty 8 S-box bits
			nxt_left_mask |= (0xFFull << 56);

			// Do the shifting
			t0 = nxt_left_mask;
			t1 = nxt_right_mask;
			nxt_right_mask = (t1 >> 61) ^ (t0 << 3);
			nxt_left_mask = ((t1 & 0x1FFFFFFFFFFFFFFFull) << 3) ^ (t0 >> 61);

			//printf(">> Level %d => %d mask : %016llX %016llX wt : %d + %d  (sbox = %d)\n", level, level+1, nxt_left_mask, nxt_right_mask, hamming_weight(nxt_left_mask), hamming_weight(nxt_right_mask), sbox);
			key_search(level+1, sbox, nxt_left_mask, nxt_right_mask);

			trail.pop_back();
		}
	}
	
	return 0;
}

void TrailDFS::key_search_wrapper() {
	u64 left_mask = 0x0ull;
	u64 right_mask = 0x0ull;
	u64 t0, t1;

	for (int level = 0; level < 2; ++level) {
		// Dirty level key bits
		if (level == 0) {
			left_mask |= (0xFull << 20);
		}
		else if (level == 1) {
			left_mask |= (0xFFFFull << 16);
		}

		printf(">> Level %d mask : %016llX %016llX wt : %d + %d\n\n", level, left_mask, right_mask, hamming_weight(left_mask), hamming_weight(right_mask));

		// Dirty 3 round counter bits
		left_mask |= 0x1ull;
		right_mask |= (0x3ull << 62);

		// Dirty 8 S-box bits
		left_mask |= (0xFFull << 56);

		// Shift around
		t0 = left_mask;
		t1 = right_mask;
		right_mask = (t1 >> 61) ^ (t0 << 3);
		left_mask = ((t1 & 0x1FFFFFFFFFFFFFFFull) << 3) ^ (t0 >> 61);
	}

	key_search(2, 5, left_mask, right_mask);
}

void TrailDFS::wrapper2() {
	u64 left_mask = 0x0ull;
	u64 right_mask = 0x0ull;
	u64 t0, t1;

	int ohkuma_sboxes[4] = { 5, 6, 9, 10 };

	// Iterate over all possible trail S-boxes for level 0 (due to Ohkuma)
	for (int level_0_iter = 0; level_0_iter < 4; ++level_0_iter) {
		// Get S-box number and clear all masks and trails
		trail.push_back(ohkuma_sboxes[level_0_iter]);
		left_mask = 0x0ull;
		right_mask = 0x0ull;

		// LEVEL 0
		left_mask |= (0xFull << (4 * trail[0]));			// Key bits
		//left_mask |= 0x1ull;								// RC bits
		//right_mask |= (0x3ull << 62);						// RC bits
		//left_mask |= (0xFFull << 56);						// S-box bits
		t0 = left_mask;										// Shift
		t1 = right_mask;
		right_mask = (t1 >> 61) ^ (t0 << 3);
		left_mask = ((t1 & 0x1FFFFFFFFFFFFFFFull) << 3) ^ (t0 >> 61);

		for (int level_1_iter = 0; level_1_iter < 3; ++level_1_iter) {
			// Put level 1 trail S-box
			trail.push_back(trail_sbox_cand[trail[0]][level_1_iter]);

			// Work with a copy so we don't mess with the mask for level 0
			u64 lm_copy = left_mask;
			u64 rm_copy = right_mask;

			// LEVEL 1
			lm_copy |= (0xFFFFull << ((trail[0] % 4) * 16));	// Key bits
			//lm_copy |= 0x1ull;									// RC bits
			//rm_copy |= (0x3ull << 62);							// RC bits
			//lm_copy |= (0xFFull << 56);							// S-box bits
			t0 = lm_copy;										// Shift
			t1 = rm_copy;
			rm_copy = (t1 >> 61) ^ (t0 << 3);
			lm_copy = ((t1 & 0x1FFFFFFFFFFFFFFFull) << 3) ^ (t0 >> 61);

			for (int level_2_iter = 0; level_2_iter < 3; ++level_2_iter) {
				// Push trail S-box for level 2
				trail.push_back(trail_sbox_cand[trail[1]][level_2_iter]);

				printf(">> Wrapper calling level 2    masks : %016llX %016llX     %016llX %016llX\n", left_mask, right_mask, lm_copy, rm_copy);
				key_search(2, trail[2], lm_copy, rm_copy);

				// Pop level 2 trail S-box
				trail.pop_back();
			}

			// Pop level 1 trail S-box
			trail.pop_back();
		}

		// Pop level 0 trail S-box
		trail.pop_back();
	}

	for (vector<int>::iterator it = trail.begin(); it != trail.end(); ++it)
		printf("%d ", *it);
}