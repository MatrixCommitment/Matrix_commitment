#include <iostream>
#include <flint/fmpz_mod.h>
#include <math.h>
#include <time.h>
#include "Vc.h"
#include "Proof.h"


Vc::Vc(){
    fmpz_init(p);
    fmpz_init(q);
    fmpz_set_str(p, "89884656743115795386465259539451236680898848947115328636715040578866337902750481566354238661203768010560056939935696678829394884407208311246423715319737062188883946712432742638151109800623047059726541476042502884419075341171231440736956555270413618581675255342293149119973622969239858152417678164812112069763", 10);
    fmpz_set_str(q, "89884656743115795386465259539451236680898848947115328636715040578866337902750481566354238661203768010560056939935696678829394884407208311246423715319737062188883946712432742638151109800623047059726541476042502884419075341171231440736956555270413618581675255342293149119973622969239858152436124908885821620397", 10);
    fmpz_init(N);
    fmpz_mul(N,p,q);

    // generate random data
    flint_rand_t rnd;
	flint_randinit(rnd);
    data = FLINT_ARRAY_ALLOC(n1, fmpz**);
	for (int i = 0; i < n1; ++i){
		data[i] = FLINT_ARRAY_ALLOC(n2, fmpz*);
        for (int j = 0; j < n2; j++){
            data[i][j] = FLINT_ARRAY_ALLOC(1, fmpz);
            fmpz_init(data[i][j]);
            fmpz_randbits(data[i][j], rnd, 256);
        }
	}

    // find g and calculate U
    fmpz_mod_ctx_init(ctxN, N);
    fmpz_init(g);
    fmpz_set_str(g,"3",10);
    fmpz_init(U);
    fmpz_init_set(U,g);
    for (int i = 0; i < n1; i++){
        fmpz_mod_pow_ui(U,U,generate_prime(i),ctxN);
    }
}

Vc::Vc(int e1, int e2){
    n1 = pow(2,e1);
    n2 = pow(2,e2);
    fmpz_init(p);
    fmpz_init(q);
    // fmpz_set_str(p, "89884656743115795386465259539451236680898848947115328636715040578866337902750481566354238661203768010560056939935696678829394884407208311246423715319737062188883946712432742638151109800623047059726541476042502884419075341171231440736956555270413618581675255342293149119973622969239858152417678164812112069763", 10);
    // fmpz_set_str(q, "89884656743115795386465259539451236680898848947115328636715040578866337902750481566354238661203768010560056939935696678829394884407208311246423715319737062188883946712432742638151109800623047059726541476042502884419075341171231440736956555270413618581675255342293149119973622969239858152436124908885821620397", 10);
    fmpz_set_str(p, "1205156213460516294290058303014157056456046623972844475679837519532628695795901600334542512053673024831724383140444002393931208489397479162484806493945387325727606669690812612385391038958840749838422771568693910028798672928952299554730693561049753982498907820671150338814736677640808714205897081983892935185184484554610795971527116005781379225040289793925450496857446141738323315590757531902436687591130253123496418949352985506262921662200616493428502380169659067", 10);
    fmpz_set_str(q, "1205156213460516294290058303014157056456046623972844475679837519532628695795901600334542512053673024831724383140444002393931208489397479162484806493945387325727606669690812612385391038958840749838422771568693910028798672928952299554730693561049753982498907820671150338814736677640808714205897081983892935185184484554610795971527116005781379225040289793925450496857446141738323315590757531902436687591130253123496418949352985506262921662200616511875246453879212529", 10);
    
    fmpz_init(N);
    fmpz_mul(N,p,q);

    // generate random data
    flint_rand_t rnd;
	flint_randinit(rnd);
    data = FLINT_ARRAY_ALLOC(n1, fmpz**);
	for (int i = 0; i < n1; ++i){
		data[i] = FLINT_ARRAY_ALLOC(n2, fmpz*);
        for (int j = 0; j < n2; j++){
            data[i][j] = FLINT_ARRAY_ALLOC(1, fmpz);
            fmpz_init(data[i][j]);
            fmpz_randbits(data[i][j], rnd, 256);
        }
	}

    // find g and calculate U
    fmpz_mod_ctx_init(ctxN, N);
    fmpz_init(g);
    fmpz_set_str(g,"3",10);
    fmpz_init(U);
    fmpz_init_set(U,g);
    for (int i = 0; i < n1; i++){
        fmpz_mod_pow_ui(U,U,generate_prime(i),ctxN);
    }
}

Vc::Vc(int e1, int e2, fmpz*** input_data){
    n1 = pow(2,e1);
    n2 = pow(2,e2);
    fmpz_init(p);
    fmpz_init(q);
    // fmpz_set_str(p, "89884656743115795386465259539451236680898848947115328636715040578866337902750481566354238661203768010560056939935696678829394884407208311246423715319737062188883946712432742638151109800623047059726541476042502884419075341171231440736956555270413618581675255342293149119973622969239858152417678164812112069763", 10);
    // fmpz_set_str(q, "89884656743115795386465259539451236680898848947115328636715040578866337902750481566354238661203768010560056939935696678829394884407208311246423715319737062188883946712432742638151109800623047059726541476042502884419075341171231440736956555270413618581675255342293149119973622969239858152436124908885821620397", 10);
    fmpz_set_str(p, "1205156213460516294290058303014157056456046623972844475679837519532628695795901600334542512053673024831724383140444002393931208489397479162484806493945387325727606669690812612385391038958840749838422771568693910028798672928952299554730693561049753982498907820671150338814736677640808714205897081983892935185184484554610795971527116005781379225040289793925450496857446141738323315590757531902436687591130253123496418949352985506262921662200616493428502380169659067", 10);
    fmpz_set_str(q, "1205156213460516294290058303014157056456046623972844475679837519532628695795901600334542512053673024831724383140444002393931208489397479162484806493945387325727606669690812612385391038958840749838422771568693910028798672928952299554730693561049753982498907820671150338814736677640808714205897081983892935185184484554610795971527116005781379225040289793925450496857446141738323315590757531902436687591130253123496418949352985506262921662200616511875246453879212529", 10);
    fmpz_init(N);
    fmpz_mul(N,p,q);

    // generate random data
    flint_rand_t rnd;
	flint_randinit(rnd);
    data = input_data;

    // find g and calculate U
    fmpz_mod_ctx_init(ctxN, N);
    fmpz_init(g);
    fmpz_set_str(g,"3",10);
    fmpz_init(U);
    fmpz_init_set(U,g);
    for (int i = 0; i < n1; i++){
        fmpz_mod_pow_ui(U,U,generate_prime(i),ctxN);
    }
}


void Vc::compute_all_local_commitment(){
    local_commitment = FLINT_ARRAY_ALLOC(n1, fmpz*);
    fmpz** si = compute_si(g, 0, n2, ctxN);
    for (int i=0; i < n1; i++){
        local_commitment[i] = FLINT_ARRAY_ALLOC(1, fmpz);
        // compute_local_commitment(data[i], n2, local_commitment[i],g, ctxN);
        compute_local_commitment(data[i], si, n2, local_commitment[i],g, ctxN);
    }
}


void Vc::hash_all_local_commitment(){
    hash_local_commitment = FLINT_ARRAY_ALLOC(n1, fmpz*);
    for (int i=0; i< n1; i++){
        hash_local_commitment[i] = FLINT_ARRAY_ALLOC(1, fmpz);
        fmpz_init(hash_local_commitment[i]);
        hash_to_fmpz(local_commitment[i], hash_local_commitment[i]);
    }
}


void Vc::generate_commitment(){

    compute_all_local_commitment();
    
    hash_all_local_commitment();

    fmpz_init(commitment);

    compute_local_commitment(hash_local_commitment, n1, commitment, g, ctxN);
    // std::cout << "Commitment generated." << std::endl;
}


Proof Vc::open_single(int row_idx, int col_idx){
    Proof result(1);
    fmpz_t lambda, si;
    fmpz_init(lambda);
    fmpz_set_si(lambda,1);
    fmpz_init_set(si,g);
    for (int i = 0; i < n2; i++){
        if (i == col_idx)
            continue;
        fmpz_t sk;
        fmpz_init_set(sk,g);
        for (size_t l = 0; l < n2; l++){
            if (l == col_idx || l == i)
                continue;
            fmpz_mod_pow_ui(sk,sk,generate_prime(l), ctxN);
        }
        fmpz_mod_pow_fmpz(sk,sk,data[row_idx][i],ctxN);
        fmpz_mod_mul(lambda, lambda, sk, ctxN);
        fmpz_mod_pow_ui(si, si, generate_prime(i), ctxN);
    }
    result.set_local_proofs(si, lambda);
    // Local proof done;

    // local_commitment[row_idx]
    result.set_local_commitment(local_commitment[row_idx]);
    // Local commitment done.

    fmpz_t Lambda, Si;
    fmpz_init(Lambda);
    fmpz_set_si(Lambda,1);
    fmpz_init_set(Si,g);
    for (int i = 0; i < n1; i++){
        if (i == row_idx)
            continue;
        fmpz_t sk;
        fmpz_init_set(sk,g);
        for (size_t l = 0; l < n2; l++){
            if (l == row_idx || l == i)
                continue;
            fmpz_mod_pow_ui(sk,sk,generate_prime(l), ctxN);
        }
        fmpz_mod_pow_fmpz(sk,sk,hash_local_commitment[i],ctxN);
        fmpz_mod_mul(Lambda, Lambda, sk, ctxN);
        fmpz_mod_pow_ui(Si, Si, generate_prime(i), ctxN);
    }
    result.set_global_proof(Si, Lambda);

    return result;
}

Proof Vc::open_multiple(int* rows, int* cols, int size){
    int num_x = count_x(rows, size);
    int* xs = new int[num_x];
    int* num_y = new int[num_x];
    int** ys = new int*[num_x];
    num_x = process_points(rows, cols, size, xs, num_y, ys);
  
    Proof result(num_x);
    
    fmpz*** local_proofs = FLINT_ARRAY_ALLOC(num_x,fmpz**);
    for (int i = 0; i < num_x; i++){
        local_proofs[i] = FLINT_ARRAY_ALLOC(2, fmpz*);
        local_proofs[i][0] = FLINT_ARRAY_ALLOC(1, fmpz);
        fmpz_init(local_proofs[i][0]);
        local_proofs[i][1] = FLINT_ARRAY_ALLOC(1, fmpz);
        fmpz_init(local_proofs[i][1]);

        open_row(xs[i],num_y[i], ys[i], local_proofs[i]);
    }
    result.set_local_proofs(local_proofs);

    result.set_local_commitment(xs, local_commitment);
    
    // std::cout << "here===" << std::endl;

    /* Calculate global proofs. */
    /* Calculate Si first. */
    fmpz_t Si;
    fmpz_init_set(Si,g);
    for (int i = 0; i < n1; i++){
        if (belongs_to(xs, i, i) > 0)
            continue;
        fmpz_mod_pow_ui(Si, Si, generate_prime(i), ctxN);
    }

    /* Calculate lambda next. */
    fmpz_t Lambda;
    fmpz_init(Lambda);
    fmpz_set_si(Lambda, 1);
    
    if (num_x == n1){
        fmpz_t s;
        fmpz_init_set(s, g);
        result.set_global_proof(s, Lambda);
        return result;
    }

    int com_size = n1 - num_x;
    int* index = new int[com_size];
    int count = 0;
    for (int i = 0; i < n1; i++){
        if (belongs_to(xs, num_x, i)==-1){
            index[count] = i;
            count++;
        }        
    }

    
    fmpz** sj = compute_si(g, index, com_size, 0, com_size, ctxN);
    fmpz** s_i = compute_si(g, 0, n1, ctxN);
    for (int i = 0; i < com_size; i++){
        fmpz_t temp;
        fmpz_init_set(temp, sj[i]);
        compute_local_commitment(data[index[i]], s_i, n2, local_commitment[i], g, ctxN);
        fmpz_mod_pow_fmpz(temp, temp, hash_local_commitment[index[i]], ctxN);
        fmpz_mod_mul(Lambda, Lambda, temp, ctxN);
    }
    
    result.set_global_proof(Si, Lambda);
    return result;
}

void Vc::open_row(int x, int num_y, int* ys, fmpz** result){
    /* Calculate Si first. */
    // std::cout << "open row begin" << std::endl;
    fmpz_t si;
    fmpz_init_set(si,g);
    for (int i = 0; i < n2; i++){
        if (belongs_to(ys, num_y, i) > 0)
            continue;
        fmpz_mod_pow_ui(si, si, generate_prime(i), ctxN);
    }

    /* Calculate lambda next. */
    fmpz_t lambda;
    fmpz_init(lambda);
    fmpz_set_si(lambda, 1);

    int com_size = n2 - num_y;
    int* index = new int[com_size];
    int count = 0;
    for (int i = 0; i < n2; i++){
        if (belongs_to(ys, num_y, i)==-1){
            index[count] = i;
            count++;
        }        
    }

    fmpz** sj = compute_si(g, index, com_size, 0, com_size, ctxN);
    for (int i = 0; i < com_size; i++){
        fmpz_t temp;
        fmpz_init_set(temp, sj[i]);
        fmpz_mod_pow_fmpz(temp, temp, data[x][index[i]], ctxN);
        fmpz_mod_mul(lambda, lambda, temp, ctxN);
    }

    fmpz_set(result[0], si);
    fmpz_set(result[1], lambda);

    // std::cout << "open row finished" << std::endl;
    return;
}


void Vc::test_open(int size){
    int* rows = new int[size];
    int* cols = new int[size];
    select_points(n1, size, rows, cols);
    
    Proof pi = open_multiple(rows, cols, size);

    verify(&pi, rows, cols, size);
}

int Vc::verify(Proof* proof, int* rows, int* cols, int size){
    int num_x = count_x(rows, size);
    int* xs = new int[num_x];
    int* num_y = new int[num_x];
    int** ys = new int*[num_x];
    num_x = process_points(rows, cols, size, xs, num_y, ys);

    for (int i = 0; i < num_x; i++){
        verify_row(xs[i], num_y[i], ys[i], proof->local_proofs[i]);
    }

    fmpz_t v1, v2;
    fmpz_init_set(v1, proof->global_proof[0]);
    fmpz_init_set(v2, proof->global_proof[1]);
    for (int i = 0; i < n1; i++){
        ulong pi = generate_prime(xs[i]);
        fmpz_mod_pow_ui(v1, v1, pi, ctxN);
        fmpz_mod_pow_ui(v2, v2, pi, ctxN);
    }
    if (fmpz_equal(v1, U) == 0){
        return 0;
    }
    fmpz** si = compute_si(proof->global_proof[0], xs, num_x, 0, num_x, ctxN);
    for (int i = 0; i < n1; i++){
        fmpz_mod_pow_fmpz(si[i], si[i], hash_local_commitment[xs[i]], ctxN);
        fmpz_mod_mul(v2, v2, si[i], ctxN);
    }
    if (fmpz_equal(v2, commitment) == 0){
        return 0;
    }
    return 1;
}

int Vc::verify(Proof* proof, int k, int l){
    fmpz_t v1, v2;
    fmpz_init_set(v1, proof->local_proofs[0][0]);
    fmpz_init_set(v2, proof->local_proofs[0][1]);

    fmpz_mod_pow_ui(v1, v1, generate_prime(l), ctxN);
    fmpz_equal(U, v1);

    fmpz_mod_pow_ui(v2, v2, generate_prime(l), ctxN);
    fmpz_t z1, z2;
    fmpz_init_set_si(z1,1);
    for (int i = 0; i < n2; i++){
        if (i == k || i == l)
            continue;
        fmpz_mod_pow_ui(z1, z1, generate_prime(i), ctxN);
    }
    fmpz_init_set(z2, z1);
    fmpz_mod_pow_ui(z1, z1, generate_prime(k), ctxN);
    fmpz_mod_pow_ui(z2, z2, generate_prime(l), ctxN);

    fmpz_mod_pow_fmpz(z1, z1, data[k][l], ctxN);
    fmpz_mod_mul_fmpz(v2, v2, z1, ctxN);
    fmpz_equal(proof->local_commitment[0], v2);

    // Verify the second half
    fmpz_mod_pow_ui(v1, proof->global_proof[0], generate_prime(l), ctxN);
    fmpz_equal(U, v1);
    fmpz_mod_pow_ui(v2, proof->global_proof[1], generate_prime(l), ctxN);
    fmpz_mod_pow_fmpz(z2, z2, hash_local_commitment[k], ctxN);
    fmpz_mod_mul_fmpz(v2, v2, z2, ctxN);
    fmpz_equal(commitment, v2);
    return 1;
}



int Vc::verify_row(int x, int num_y, int* ys, fmpz** proofs){
    fmpz_t v1, v2;
    fmpz_init_set(v1, proofs[0]);
    fmpz_init_set(v2, proofs[1]);
    for (int i = 0; i < n2; i++){
        ulong pi = generate_prime(ys[i]);
        fmpz_mod_pow_ui(v1, v1, pi, ctxN);
        fmpz_mod_pow_ui(v2, v2, pi, ctxN);
    }
    if (fmpz_equal(v1, U) == 0){
        return 0;
    }
    fmpz** si = compute_si(proofs[0], ys, num_y, 0, num_y, ctxN);
    for (int i = 0; i < n2; i++){
        fmpz_mod_pow_fmpz(si[i], si[i], data[x][ys[i]], ctxN);
        fmpz_mod_mul(v2, v2, si[i], ctxN);
    }
    if (fmpz_equal(v2, local_commitment[x]) == 0){
        return 0;
    }
    return 1;
}

fmpz** Vc::__aggregate(fmpz** vec, int* pos1, int size1, fmpz** proof1, 
                            int* pos2, int size2, fmpz** proof2){
    fmpz_t e1, e2;
    fmpz_init_set_ui(e1, 1);
    fmpz_init_set_si(e2, 1);
    for (int i = 0; i < size1; i++){
        fmpz_mul_si(e1, e1, generate_prime(pos1[i]));
    }
    for (int i = 0; i < size2; i++){
        fmpz_mul_si(e2, e2, generate_prime(pos2[i]));
    }


    fmpz** result = FLINT_ARRAY_ALLOC(2, fmpz*);
    result[0] = FLINT_ARRAY_ALLOC(1, fmpz);
    fmpz_init(result[0]);
    result[1] = FLINT_ARRAY_ALLOC(1, fmpz);
    fmpz_init(result[1]);

    shamir_trick(result[0], proof1[0], proof2[0], e1, e2, ctxN);

    fmpz** phi_j = compute_si(result[0], pos2, size2, 0, size2, ctxN);
    fmpz** psi_i = compute_si(result[0], pos1, size1, 0, size1, ctxN);

    fmpz_t den_1, den_2;
    fmpz_init_set_si(den_1, 1);
    fmpz_init_set_si(den_2, 1);
    for (int j = 0; j < size2; j++){
        fmpz_t temp;
        fmpz_init(temp);
        fmpz_mod_pow_fmpz(temp, phi_j[j], vec[pos2[j]], ctxN);
        fmpz_mod_mul(den_1, den_1, temp, ctxN);
    }
    for (int j = 0; j < size1; j++){
        fmpz_t temp;
        fmpz_init(temp);
        fmpz_mod_pow_fmpz(temp, psi_i[j], vec[pos1[j]], ctxN);
        fmpz_mod_mul(den_1, den_1, temp, ctxN);
    }
    fmpz_mod_inv(den_1, den_1, ctxN);
    fmpz_mod_inv(den_2, den_2, ctxN);
    fmpz_mod_mul_fmpz(den_1, den_1, proof1[1], ctxN);
    fmpz_mod_mul_fmpz(den_2, den_2, proof2[1], ctxN);
    shamir_trick(result[1], den_1, den_2, e1, e2, ctxN);
    return result;
}


Proof Vc::aggregate(int *x1, int *y1, int size_1, Proof *p1, 
                    int *x2, int *y2, int size_2, Proof *p2){
    // Processing input:
    int num_x1 = count_x(x1, size_1);
    int* xs1 = new int[num_x1];
    int* num_y1 = new int[num_x1];
    int** ys1 = new int*[num_x1];
    num_x1 = process_points(x1, y1, size_1, xs1, num_y1, ys1);

    int num_x2 = count_x(x2, size_2);
    int* xs2 = new int[num_x2];
    int* num_y2 = new int[num_x2];
    int** ys2 = new int*[num_x2];
    num_x2 = process_points(x2, y2, size_2, xs2, num_y2, ys2);

    int size = size_1 + size_2;
    int* x = new int[size];
    int* y = new int[size];
    mergeAndSortPoints(x1, y1, size_1, x2, y2, size_2, x, y);

    int num_x = count_x(x, size);
    int* xs = new int[num_x];
    int* num_y = new int[num_x];
    int** ys = new int*[num_x];
    num_x = process_points(x, y, size, xs, num_y, ys);

    Proof retval(num_x);
    fmpz*** local_proofs = FLINT_ARRAY_ALLOC(num_x, fmpz**);
    for (int i = 0; i < num_x; i++){
        local_proofs[i] = FLINT_ARRAY_ALLOC(2, fmpz*);
        local_proofs[i][0] = FLINT_ARRAY_ALLOC(1, fmpz);
        local_proofs[i][1] = FLINT_ARRAY_ALLOC(1, fmpz);
        fmpz_init(local_proofs[i][0]);
        fmpz_init(local_proofs[i][1]);
    }


    // Aggregation for local proofs
    int i = 0;
    int j = 0;
    int k = 0;
    while (1){
        int cur_x = xs[k];
        if (cur_x == xs1[i] && cur_x != xs2[j]){
            // directly copy the proof.
            fmpz_set(local_proofs[k][0], p1->local_proofs[i][0]);
            fmpz_set(local_proofs[k][1], p1->local_proofs[i][1]);
            fmpz_set(retval.local_commitment[k], p1->local_commitment[i]);
            i++;
            k++;
        } else if (cur_x != xs1[i] && cur_x == xs2[j]){
            // Copy the proof here
            fmpz_set(local_proofs[k][0], p2->local_proofs[j][0]);
            fmpz_set(local_proofs[k][1], p2->local_proofs[j][1]);
            fmpz_set(retval.local_commitment[k], p2->local_commitment[j]);
            j++; 
            k++;
        } else /*(cur_x == xs1[i] == xs2[j])*/{
            /* Aggregate them */
            fmpz** new_pf = __aggregate(data[cur_x], ys1[i], num_y1[i], p1->local_proofs[i], ys2[j], num_y2[j], p2->local_proofs[j]);
            fmpz_set(local_proofs[k][0], new_pf[0]);
            fmpz_set(local_proofs[k][1], new_pf[1]);
            fmpz_set(retval.local_commitment[k], p1->local_commitment[i]);
            i++;
            j++;
            k++;
        }
        if (k == num_x){
            break;
        } 
    }
    retval.set_local_proofs(local_proofs);

    // std::cout << "local proof aggregated." << std::endl;

    // Now aggregate the global proof
    // fmpz** new_global_pf = __aggregate(hash_local_commitment, xs1, num_x1, p1->global_proof, xs2, num_x2, p2->global_proof);
    // fmpz_set(retval.global_proof[0], new_global_pf[0]);
    // fmpz_set(retval.global_proof[1], new_global_pf[1]);
    fmpz_set(retval.global_proof[0], p1->global_proof[0]);
    fmpz_set(retval.global_proof[1], p1->global_proof[1]);

    return retval;
}


Proof Vc::agg_many_to_one(int* x, int* y, int size, Proof* proofs){
    if (size == 1){
        return proofs[0];
    }

    int size_1 = size/2;
    int size_2 = size - size_1;

    Proof p1 = agg_many_to_one(x, y, size_1, proofs);
    Proof p2 = agg_many_to_one(x+size_1, y+size_1, size_2, proofs+size_1);

    Proof p = aggregate(x, y, size_1, &p1, x+size_1, y+size_1, size_2, &p2);

    return p;
}