#ifndef VC_H
#define VC_H

#include <iostream>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/ulong_extras.h>
#include <flint/fmpz_mod.h>
#include <math.h>
#include <time.h>
#include <string>
#include <openssl/evp.h>
#include "util.h"
#include "Proof.h"


class Vc{
    public:
        int n1 = 128;                    // The number of row.
        int n2 = 128;                    // The number of coloum.
        int prime_length = 64;     // The lenth of different primes. (in bits)
        int block_size = 256;      // The size of vector block.
        int rsa_m = 2048;          // The size of rsa modulus.

        fmpz_t p;   // 1024-bit prime
        fmpz_t q;   // 1024-bit prime
        fmpz_t N;   // 2048-bit rsa modulus
        fmpz_mod_ctx_t ctxN;

        fmpz_t g;   // The generator of Z_N
        fmpz_t U;   // Accumulator of the primes

        fmpz*** data;             // The committed matrix
        fmpz_t commitment;        // The commitment of the matrix
        fmpz*** local_proof;
        fmpz** local_commitment;
        fmpz** hash_local_commitment;
        fmpz** global_proof;

        Vc();
        Vc(int e1, int e2);
        Vc(int e1, int e2, fmpz*** input_data);

        
        void compute_all_local_commitment();
        void hash_all_local_commitment();

        int verify_row(int x, int num_y, int* ys, fmpz** proofs);

        void open_row(int x, int num_y, int* ys, fmpz** result);
        void open_row(int x, int num_y, int* ys, fmpz** result, fmpz** Si);
        void test_open(int size);

        void update_all_proof(int x, int y, fmpz_t delta);
        void prepare_update(int row_idx, int col_idx);

        fmpz** __aggregate(fmpz** vec, int* pos1, int size1, fmpz** proof1, 
                            int* pos2, int size2, fmpz** proof2);


        // Interface for outer call
        void generate_commitment();
        Proof open_single(int row_idx, int col_idx);
        Proof open_multiple(int* rows, int* cols, int size);

        int verify(Proof* proof, int* rows, int* cols, int size);
        int verify(Proof* proof, int k, int l);

        Proof aggregate(int* x1, int* y1, int size_1, Proof* p1,
                        int* x2, int* y2, int size_2, Proof* p2);
        Proof agg_many_to_one(int* x, int* y, int size, Proof* proofs);
};

#endif