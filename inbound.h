#ifndef INBOUND_H
#define INBOUND_H

#include "present.h"
#include <set>

extern u64 sbox_one_bit[2][4][8];
const int NUM_ROUNDS = 3;
const int L_OUT_SIZE = 2*14605415;

class inbound {
private:
	///////////////////////////////////////////////////////
	//// DATA STRUCTURES							   ////
	///////////////////////////////////////////////////////
	std::set<u64> W; 

	std::set<u64> layerIn[NUM_ROUNDS], layerOut[NUM_ROUNDS];
	u64 *maskIn, *maskOut, *keyMaskLeft, *keyMaskRight;

	u64 round_keys[NUM_ROUNDS][2];

	u64 k0_mask[2];
	u64 k0_val[2];

	// For holding the definition of trail positions and the S-boxes they pass through
	int *trailSbox, *trailPos;
	
	///////////////////////////////////////////////////////
	//// FOR CONSTRUCTING W AND L_OUT				   ////
	///////////////////////////////////////////////////////
	void offset_W_layer_byKey(int);
	void determine_W();		
	u64 get_input_mask(u64);
	void do_layer(int, u64, u64);
	void construct_L_out();

	///////////////////////////////////////////////////////
	//// INITIALIZATION FUNCTIONS 					   ////
	///////////////////////////////////////////////////////
	void setup_trails();
	void setup_masks();
	void setup_fixed_key_bits();

	///////////////////////////////////////////////////////
	//// CODE RELATED TO DETERMINING CONSTRAINTS ON k3 ////
	///////////////////////////////////////////////////////
	void determine_k0_constraints();
	void inv_ks_helper(u64 &ml, u64 &mr, u64 &vl, u64 &vr, u64 rc);

	///////////////////////////////////////////////////////
	//// TESTING AND VALIDATION 					   ////
	///////////////////////////////////////////////////////
	int verify_W_set();
	int verify_k0_constraints(u64, u64);
	int check_k0_mask_val(u64, u64, u64, u64);
	int test_value(u64);
	void verify_W_layer(int);

public:	
	std::set<std::pair<u64, key_state>> S;
	std::set<std::pair<u64, key_state>> L_out;

	inbound();
	~inbound();

	void run_phase();
};

#endif // INBOUND_H
