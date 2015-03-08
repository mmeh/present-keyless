#ifndef PRESENT_H
#define PRESENT_H

#include <stdint.h>

#define PRESENT80 0

typedef uint64_t u64;
typedef uint8_t  u8;

extern int P[64];
extern int P_inv[64];
extern u64 S[16];
extern u64 S_inv[16];

extern u64 S_compact[8][256];
extern u64 S_inv_compact[8][256];
extern u64 P_compact[8][256];
extern u64 P_inv_compact[8][256];

struct key_state {
	u64 left;
	u64 right;
};

extern bool operator<(const key_state &, const key_state &);

u64 sub_layer(u64);
u64 inv_sub_layer(u64);
u64 perm_layer(u64);
u64 inv_perm_layer(u64);

void encrypt_one_round(u64 &, u64);
void decrypt_one_round(u64 &, u64);

void dirty_rc_bits(u64 &, u64 &);
void dirty_sbox_bits(u64 &);

u8 parity(u64);
int weight(u64);

// PRESENT-128
void rotate_keys_128(u64 &, u64 &);
void inv_rotate_keys_128(u64 &, u64 &);
void inv_keysched_round_128(u64 &, u64 &, u64);
void keysched_round_128(u64 &, u64 &, u64);

// PRESENT-80
void rotate_keys_80(u64 &, u64 &);
void inv_rotate_keys_80(u64 &, u64 &);
void inv_keysched_round_80(u64 &, u64 &, u64);
void keysched_round_80(u64 &, u64 &, u64);

// WRAPPERS
void rotate_keys(u64 &, u64 &);
void inv_rotate_keys(u64 &, u64 &);
void inv_keysched_round(u64 &, u64 &, u64);
void keysched_round(u64 &, u64 &, u64);

#endif
