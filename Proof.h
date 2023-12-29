#ifndef PROOF_H
#define PROOF_H

#include <iostream>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/ulong_extras.h>
#include <flint/fmpz_mod.h>
#include "util.h"

class Proof
{
public:
    int num_row;
    fmpz*** local_proofs;     /* m*2 matrix */
    fmpz** local_commitment;  /* m*1 array */
    fmpz** global_proof;      /* 1*2 array */
    Proof();
    Proof(int rows);
    ~Proof();
    void set_local_proofs(fmpz_t si, fmpz_t lambda);
    void set_local_proofs(fmpz*** proofs);
    void set_local_commitment(fmpz_t cmt);
    void set_local_commitment(int* rows, fmpz** cmt);
    void set_global_proof(fmpz_t si, fmpz_t lambda);


};





#endif