#ifndef TRAILDFS_H
#define TRAILDFS_H

#include <stdint.h>
#include <set>
#include <vector>

using namespace std;

typedef uint64_t u64;

 int trail_sbox_cand[17][16] = { 
	 { 1, 2, 3 }, { 5, 6, 7 }, { 9, 10, 11 }, { 13, 14, 15 }, 
	 { 1, 2, 3 }, { 5, 6, 7 }, { 9, 10, 11 }, { 13, 14, 15 }, 
	 { 1, 2, 3 },  { 5, 6, 7 }, { 9, 10, 11 }, { 13, 14, 15 }, 
	 { 1, 2, 3 }, { 5, 6, 7 }, { 9, 10, 11 }, { 13, 14, 15 },
	 { 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 }
 };

u64 levelMask[16] = { 0x0001000100010001, 0x0002000200020002, 0x0004000400040004, 0x0008000800080008, 0x0010001000100010, 0x0020002000200020,
0x0040004000400040, 0x0080008000800080, 0x0100010001000100, 0x0200020002000200, 0x0400040004000400, 0x0800080008000800,
0x1000100010001000, 0x2000200020002000, 0x4000400040004000, 0x8000800080008000 };

class TrailDFS
{
private:
	int hamming_weight(u64);
	vector<int> trail;
	void print_trail();
	void test_in_mask(u64, u64);

public:
	int key_search(int, int, u64, u64);
	void key_search_wrapper();
	void wrapper2();
};

#endif