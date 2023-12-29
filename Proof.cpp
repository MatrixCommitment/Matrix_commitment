#include "Proof.h"


Proof::Proof(){
    int rows = 1;
    num_row = rows;
    local_proofs = FLINT_ARRAY_ALLOC(rows, fmpz**);
    local_commitment = FLINT_ARRAY_ALLOC(rows, fmpz*);
    for (size_t i = 0; i < rows; i++){
        local_proofs[i] = FLINT_ARRAY_ALLOC(2, fmpz*);
        local_proofs[i][0] = FLINT_ARRAY_ALLOC(1, fmpz);
        local_proofs[i][1] = FLINT_ARRAY_ALLOC(1, fmpz);
        local_commitment[i] = FLINT_ARRAY_ALLOC(1, fmpz);
        fmpz_init(local_proofs[i][0]);
        fmpz_init(local_proofs[i][1]);
        fmpz_init(local_commitment[i]);
    }
    global_proof = FLINT_ARRAY_ALLOC(2, fmpz*);
    global_proof[0] = FLINT_ARRAY_ALLOC(1, fmpz);
    fmpz_init(global_proof[0]);
    global_proof[1] = FLINT_ARRAY_ALLOC(1, fmpz);
    fmpz_init(global_proof[1]);
}


Proof::Proof(int rows)
{
    num_row = rows;
    local_proofs = FLINT_ARRAY_ALLOC(rows, fmpz**);
    local_commitment = FLINT_ARRAY_ALLOC(rows, fmpz*);
    for (size_t i = 0; i < rows; i++){
        local_proofs[i] = FLINT_ARRAY_ALLOC(2, fmpz*);
        local_proofs[i][0] = FLINT_ARRAY_ALLOC(1, fmpz);
        local_proofs[i][1] = FLINT_ARRAY_ALLOC(1, fmpz);
        local_commitment[i] = FLINT_ARRAY_ALLOC(1, fmpz);
        fmpz_init(local_proofs[i][0]);
        fmpz_init(local_proofs[i][1]);
        fmpz_init(local_commitment[i]);
    }
    global_proof = FLINT_ARRAY_ALLOC(2, fmpz*);
    global_proof[0] = FLINT_ARRAY_ALLOC(1, fmpz);
    fmpz_init(global_proof[0]);
    global_proof[1] = FLINT_ARRAY_ALLOC(1, fmpz);
    fmpz_init(global_proof[1]);
}

Proof::~Proof()
{
}


/* Set local proof for single element. */
void Proof::set_local_proofs(fmpz_t si, fmpz_t lambda)
{
     fmpz_init_set(local_proofs[0][0], si);
     fmpz_init_set(local_proofs[0][1], lambda);
}

void Proof::set_local_proofs(fmpz*** proofs){
    for (int i = 0; i < num_row; i++){
        fmpz** row_proof = proofs[i];
        fmpz_set(local_proofs[i][0], row_proof[0]);
        fmpz_set(local_proofs[i][1], row_proof[1]);
    }
    return;
}


/* Set local commitment for single element. */
void Proof::set_local_commitment(fmpz_t cmt)
{
    fmpz_init_set(local_commitment[0], cmt);
}

void Proof::set_local_commitment(int* rows, fmpz** cmt){
    for (int i = 0; i < num_row; i++){
        fmpz_set(local_commitment[i], cmt[rows[i]]);
    }
    return;
}

void Proof::set_global_proof(fmpz_t si, fmpz_t lambda)
{
    fmpz_init_set(global_proof[0], si);
    fmpz_init_set(global_proof[1], lambda);
}


