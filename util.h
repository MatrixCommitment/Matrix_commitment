#ifndef UTIL_H
#define UTIL_H

#include <iostream>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/ulong_extras.h>
#include <flint/fmpz_mod.h>
#include <math.h>
#include <time.h>
#include <string>
#include <openssl/evp.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
extern "C" {
#include <relic/relic.h>
}


/* Hash a large integer or group element to a 256-bit number. */
void hash_to_fmpz(fmpz_t num, fmpz_t retval);

/* Maps each integer to different primes. */
ulong generate_prime(int idx);

/* Generate the rsa commitment for a single vector. */
void compute_local_commitment(fmpz** vector, int size, fmpz_t result, fmpz_t g, fmpz_mod_ctx_t ctxN);
void compute_local_commitment(fmpz** vector, fmpz** si, int size, fmpz_t result, fmpz_t g, fmpz_mod_ctx_t ctxN);

int belongs_to(int arr[], int n, int target);

void selectRandomNumbers(int* randomNumbers, int n, int m);

int count_x(int* rows, int size);

int process_points(int* rows, int* cols, int size, int* xs, int* num_y, int** ys);

void select_points(int n, int m, int* x, int* y);

/* Calculate n Si in O(nlogn) time. */
fmpz** compute_si(fmpz_t g, int left, int right, fmpz_mod_ctx_t ctx);
fmpz** compute_si(fmpz_t g, int* index, int size, int left, int right, fmpz_mod_ctx_t ctx);

fmpz** merge(fmpz** left_result, fmpz** right_result, int left, int mid, int right);


bool is_generator(fmpz_t g, fmpz_t p, fmpz_t q);


void generate_secure_prime(fmpz_t p, flint_rand_t rnd);


void find_generator(fmpz_t g, fmpz_t p, fmpz_t q);


void shamir_trick(fmpz_t result, fmpz_t r_x, fmpz_t r_y, fmpz_t x, fmpz_t y, fmpz_mod_ctx_t ctx);


int* mergeAndSort(int* array1, int* array2, int size1, int size2);


void fmpz2bn(fmpz_t in, bn_t out);


void bn2fmpz(fmpz_t out, bn_t in);


void mergeAndSortPoints(int* x1, int* y1, int size_1, int* x2, int* y2, int size_2, int* x, int* y);


void hash1(int i, ep_t c, int* xs, int size, ep_t* ci, fmpz_t retval);
void hash2(int i, int size, ep_t* ci, int** ys, fmpz_t retval);
void hash3(int j, ep_t ci, int* xs, int size, fmpz_t retval);

#endif