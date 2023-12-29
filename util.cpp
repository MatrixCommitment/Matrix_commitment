#include "util.h"

ulong generate_prime(int idx){
    // if (DEBUG){
    //     return (ulong)(idx);
    // }
    unsigned char result[32] = {0};
    std::string input = std::__cxx11::to_string(idx);
    EVP_Digest(input.c_str(), input.size(), result, NULL, EVP_sha256(), NULL);
    ulong* res = (ulong*)result;    
    ulong retval = n_nextprime(*res, 1);
    // std::cout << retval << std::endl;
    return retval;
}

// ulong generate_prime(int idx){
//     unsigned char result[32] = {0};
//     std::string input_1 = std::__cxx11::to_string(idx);
//     for (int i = 1; ;i++){
//         std::string input_2 = std::__cxx11::to_string(i);
//         std::string input = input_1 + input_2;
//         EVP_Digest(input.c_str(), input.size(), result, NULL, EVP_sha256(), NULL);
//         ulong* res = (ulong*)result;
//         if (n_is_prime(*res)){
//             return *res;
//         }
//     }
// }


void hash_to_fmpz(fmpz_t num, fmpz_t retval){
	char* sinput = fmpz_get_str(NULL, 16, num);
	std::string input = sinput;
	int length = input.size();
	unsigned char result[32] = {0};
	EVP_Digest(input.c_str(), length, result, NULL, EVP_sha256(), NULL);
	unsigned int * uptr = (unsigned int *)result;
	long int arr[8] = {0};
	for (int i = 0; i < 8; ++i){
		arr[i] = (long)*uptr;
		uptr ++;
	}
	fmpz_set_ui(retval,0);
	for (int i = 0; i < 8; ++i){
		fmpz_t local;
		fmpz_init(local);
		fmpz_set_si(local, arr[i]);
		fmpz_mul_2exp(local, local, 32*(7-i));
		fmpz_add(retval, retval, local);
		fmpz_clear(local);
	}
}

void compute_local_commitment(fmpz** vector, int size, fmpz_t result, fmpz_t g, fmpz_mod_ctx_t ctxN){
    fmpz_t retval;
    fmpz_init(retval);
    fmpz_set_ui(retval,1);
    for (int i=0; i < size; i++){
        // Compute s_i
        fmpz_t si;
        fmpz_init_set(si,g);
        for (int j = 0; j < size; j++){
            if (j==i)
                continue;
            fmpz_mod_pow_ui(si,si,generate_prime(j),ctxN);
        }
        fmpz_mod_pow_fmpz(si,si,vector[i], ctxN);
        fmpz_mod_mul(retval, retval, si, ctxN);
    }
    fmpz_init_set(result, retval);
}

void compute_local_commitment(fmpz** vector, fmpz** si, int size, fmpz_t result, fmpz_t g, fmpz_mod_ctx_t ctxN){
    fmpz_t retval;
    fmpz_init(retval);
    fmpz_set_ui(retval,1);
    for (int i=0; i < size; i++){
        // Compute s_i
        fmpz_t Si;
        fmpz_init_set(Si,si[i]);
        fmpz_mod_pow_fmpz(Si,Si,vector[i], ctxN);
        fmpz_mod_mul(retval, retval, Si, ctxN);
    }
    fmpz_init_set(result, retval);
}


int belongs_to(int arr[], int n, int target){
    int l = 0, r = n - 1;
    while (l <= r) {
        int mid = l + (r - l) / 2;
        if (arr[mid] == target) {
            return mid;
        } else if (arr[mid] < target) {
            l = mid + 1;
        } else {
            r = mid - 1;
        }
    }
    return -1;
}


// 生成 [0, n-1] 之间的随机整数
int generateRandomNumber(int n)
{
    return rand() % n;
}

// 选择大小为m的不重复随机数
void selectRandomNumbers(int* randomNumbers, int n, int m)
{

    // 初始化随机数数组，全部置为-1
    for (int i = 0; i < m; i++) {
        randomNumbers[i] = -1;
    }

    // 选出m个不一样的随机数
    int i = 0;
    while (i < m) {
        int randomNumber = generateRandomNumber(n);
        // 遍历随机数数组，判断当前数是否重复
        int j = 0;
        for (; j < i; j++) {
            if (randomNumbers[j] == randomNumber) {
                break;
            }
        }
        // 如果不重复，把数添加到随机数数组中
        if (j == i) {
            randomNumbers[i] = randomNumber;
            i++;
        }
    }

    // 对随机数数组进行排序
    for (int i = 0; i < m - 1; i++) {
        for (int j = i + 1; j < m; j++) {
            if (randomNumbers[i] > randomNumbers[j]) {
                int temp = randomNumbers[i];
                randomNumbers[i] = randomNumbers[j];
                randomNumbers[j] = temp;
            }
        }
    }

    return;
}


int count_x(int* rows, int size){
    int num_x = 1;
    for (int i = 1; i < size; i++){
        if (rows[i]!=rows[i-1])
            num_x ++;
    }
    return num_x;
}


int process_points(int* rows, int* cols, int size, int* xs, int* num_y, int** ys){
    int num_x = 1;
    for (int i = 1; i < size; i++){
        if (rows[i]!=rows[i-1])
            num_x ++;
    }
    xs[0] = rows[0];
    num_y[0] = 1;
    int index = 0;
    for (int i = 1; i < size; i++){
        if (rows[i]!=rows[i-1]){
            index ++;
            xs[index] = rows[i];
            num_y[index] = 0;
        }
        num_y[index] ++;
    }
    
    int total = 0;
    for (int i = 0; i < num_x; i++){
        ys[i] = new int[num_y[i]];
        for (int j = 0; j < num_y[i]; j++){
            ys[i][j] = cols[total];
            total++;
        }
    }
    // for (int i = 0; i < num_x; i++){
    //         std::cout << xs[i] << " " << num_y[i] << std::endl;
    // }
    return num_x;
}


void select_points(int n, int m, int* x, int* y) {
    int i, j, k;
    int selected[m], temp;
    srand(time(NULL));
    for (i = 0; i < m; i++) {
        selected[i] = rand() % (n * n);
        // std::cout << selected[i] << std::endl;
        for (j = 0; j < i; j++) {
            if (selected[i] == selected[j]) {
                i--;
                break;
            }
        }
    }
    for (i = 0; i < m; i++) {
        x[i] = selected[i] % n;
        y[i] = selected[i] / n;
    }
    for (i = 0; i < m - 1; i++) {
        for (j = i + 1; j < m; j++) {
            if (x[i] > x[j]) {
                temp = x[i];
                x[i] = x[j];
                x[j] = temp;
                temp = y[i];
                y[i] = y[j];
                y[j] = temp;
            } else if (x[i] == x[j] && y[i] > y[j]) {
                temp = y[i];
                y[i] = y[j];
                y[j] = temp;
            }
        }
    }
}


/* Left inclusive right exclusive*/
fmpz** compute_si(fmpz_t g, int left, int right, fmpz_mod_ctx_t ctx){
    if (right == left + 1){
        fmpz** retval = FLINT_ARRAY_ALLOC(1, fmpz*);
        retval[0] = FLINT_ARRAY_ALLOC(1, fmpz);
        fmpz_init_set(retval[0], g);
        return retval; 
    }

    int mid = (left + right) / 2;

    fmpz_t g_left;
    fmpz_init_set(g_left, g);
    for (int i = left; i < mid; i++){
        fmpz_mod_pow_ui(g_left, g_left, generate_prime(i), ctx);
    }
    fmpz_t g_right;
    fmpz_init_set(g_right, g);
    for (int i = mid; i < right; i++){
        fmpz_mod_pow_ui(g_right, g_right, generate_prime(i), ctx);
    }
    fmpz** left_result = compute_si(g_right, left, mid, ctx);
    fmpz** right_result = compute_si(g_left, mid, right, ctx);

    return merge(left_result, right_result, left, mid, right);
}


/* Left inclusive right exclusive
    for non-consecutive indices. */
fmpz** compute_si(fmpz_t g, int* index, int size, int left, int right, fmpz_mod_ctx_t ctx){
    if (right == left + 1){
        fmpz** retval = FLINT_ARRAY_ALLOC(1, fmpz*);
        retval[0] = FLINT_ARRAY_ALLOC(1, fmpz);
        fmpz_init_set(retval[0], g);
        return retval;
    }

    int mid = (left + right) / 2;

    fmpz_t g_left;
    fmpz_init_set(g_left, g);
    for (int i = left; i < mid; i++){
        fmpz_mod_pow_ui(g_left, g_left, generate_prime(index[i]), ctx);
    }
    fmpz_t g_right;
    fmpz_init_set(g_right, g);
    for (int i = mid; i < right; i++){
        fmpz_mod_pow_ui(g_right, g_right, generate_prime(index[i]), ctx);
    }
    fmpz** left_result = compute_si(g_right, index, size, left, mid, ctx);
    fmpz** right_result = compute_si(g_left, index, size, mid, right, ctx);

    return merge(left_result, right_result, left, mid, right);
}


fmpz** merge(fmpz** left_result, fmpz** right_result, int left, int mid, int right){
    int size = right - left;
    int left_size = mid - left;
    int right_size = right - mid;
    
    // std::cout << left << "-" << mid << "-" << right << std::endl;
    fmpz** retval = FLINT_ARRAY_ALLOC(size, fmpz*);
    for (int i = 0; i < size; i++){
        retval[i] = FLINT_ARRAY_ALLOC(1, fmpz);
        if (i < left_size){
            fmpz_init_set(retval[i], left_result[i]);
        } else {
            fmpz_init_set(retval[i], right_result[i-left_size]);
        }
    }
    // std::cout << "finished" << std::endl;
    return retval;
}

bool is_generator(fmpz_t g, fmpz_t p, fmpz_t q){
    fmpz_t phi, p1,q1;
    fmpz_init(phi);
    fmpz_init(p1);
    fmpz_init(q1);
    fmpz_sub_si(p1, p, 1);
    fmpz_sub_si(q1, q, 1);
    fmpz_mul(phi, p1, q1);

    std::cout << "Here" << std::endl;
    fmpz_factor_t fac;
    fmpz_factor_init(fac);
    fmpz_factor(fac, phi);
    std::cout << "There" << std::endl;

    fmpz_t N;
    fmpz_init(N);
    fmpz_mul(N,p,q);

    fmpz_mod_ctx_t ctx;
    fmpz_mod_ctx_init(ctx, N);

    for (int i = 0; i < fac->num; i++){
        fmpz_t p;
        fmpz_init_set(p, fac->p+i);
        fmpz_t exp;
        fmpz_init(exp);
        fmpz_divexact(exp, phi, p);

        fmpz_t temp;
        fmpz_init(temp);
        fmpz_mod_pow_fmpz(temp, g, exp, ctx);
        if (fmpz_equal_si(temp, 1)){
            return false;
        }
    }
    return true;
}


void generate_secure_prime(fmpz_t p, flint_rand_t rnd){
    fmpz_t result;
    fmpz_init(result);
    while (1){
        fmpz_randprime(result, rnd, 1024, 1);
        fmpz_t t;
        fmpz_init(t);
        fmpz_sub_ui(t, result,1);
        fmpz_divexact_ui(t,t,2);
        // fmpz_print(result);
        // std::cout << "\n";
        if (fmpz_is_prime(t)){
            break;
        }
    }
    fmpz_set(p, result);
}


void find_generator(fmpz_t g, fmpz_t p, fmpz_t q){
    
    fmpz_t result;
    fmpz_init(result);
    fmpz_set_ui(result, 1);
    
    fmpz_t N;
    fmpz_init(N);
    fmpz_mul(N,p,q);
    
    fmpz_mod_ctx_t ctx;
    fmpz_mod_ctx_init(ctx, N);
    fmpz_t phi, p1,q1;
    fmpz_init(phi);
    fmpz_init(p1);
    fmpz_init(q1);
    fmpz_sub_si(p1, p, 1);
    fmpz_sub_si(q1, q, 1);
    fmpz_mul(phi, p1, q1);
    fmpz_divexact_ui(p1,p1,2);
    fmpz_divexact_ui(q1,q1,2);

    
    while (1){
        fmpz_add_ui(result, result, 1);
        fmpz_print(result);
        std::cout << "\n";

        if (fmpz_equal(result, N)){
            break;
        }

        fmpz_t r;
        fmpz_init(r);
        fmpz_mod(r,result,p);
        if (fmpz_equal_si(r, 0)){
            continue;
        }
        fmpz_mod(r,result,q);
        if (fmpz_equal_si(r, 0)){
            continue;
        }

        fmpz_t exp;
        fmpz_init(exp);
        fmpz_divexact(exp, phi, p1);

        fmpz_t temp;
        fmpz_init(temp);
        fmpz_mod_pow_fmpz(temp, result, exp, ctx);
        if (fmpz_equal_si(temp, 1)){
            continue;
        }

        fmpz_divexact(exp, phi, q1);
        fmpz_mod_pow_fmpz(temp, result, exp, ctx);
        if (fmpz_equal_si(temp, 1)){
            continue;
        }

        fmpz_divexact_ui(exp, phi, 2);
        fmpz_mod_pow_fmpz(temp, result, exp, ctx);
        if (fmpz_equal_si(temp, 1)){
            continue;
        }

        fmpz_set(g, result);
        break;
    }
}


void shamir_trick(fmpz_t result, fmpz_t r_x, fmpz_t r_y, fmpz_t x, fmpz_t y, fmpz_mod_ctx_t ctx){
    fmpz_t d,a,b;
    fmpz_init(d);
    fmpz_init(a);
    fmpz_init(b);
    fmpz_xgcd_canonical_bezout(d,a,b,x,y);
    fmpz_mod_pow_fmpz(result, r_x, b, ctx);
    fmpz_t temp;
    fmpz_init(temp);
    fmpz_mod_pow_fmpz(temp, r_y, a, ctx);
    fmpz_mod_mul(result, result, temp, ctx);
    return;
}


int* mergeAndSort(int* array1, int* array2, int size1, int size2) {
    int* mergedArray = new int[size1+size2];
    int i = 0, j = 0, k = 0;

    while (i < size1 && j < size2) {
        if (array1[i] < array2[j])
            mergedArray[k++] = array1[i++];
        else 
            mergedArray[k++] = array2[j++];
    }

    while (i < size1)
        mergedArray[k++] = array1[i++];
    
    while (j < size2)
        mergedArray[k++] = array2[j++];
    
    return mergedArray;
}


void fmpz2bn(fmpz_t in, bn_t out){
    bn_new(out);
    size_t size = fmpz_sizeinbase(in, 10);
    char* buffer = new char[size+2];
    freopen("/dev/null", "a", stdout);
    setbuf(stdout, buffer);
    fmpz_print(in);
    freopen ("/dev/tty", "a", stdout);
    bn_read_str(out, buffer, size+2, 10);
    fcloseall();
}


void bn2fmpz(fmpz_t out, bn_t in){
    int size = bn_size_str(in, 10);
    char *in_char = new char[size];
    bn_write_str(in_char, size, in, 10);
    fmpz_set_str(out, in_char, 10);
}


void mergeAndSortPoints(int* x1, int* y1, int size_1, int* x2, int* y2, int size_2, int* x, int* y) {
    std::vector<std::pair<int, int>> points;
    
    for (int i = 0; i < size_1; ++i) {
        points.emplace_back(x1[i], y1[i]);
    }
    
    for (int i = 0; i < size_2; ++i) {
        points.emplace_back(x2[i], y2[i]);
    }
    
    std::sort(points.begin(), points.end());
    
    for (int i = 0; i < size_1 + size_2; ++i) {
        x[i] = points[i].first;
        y[i] = points[i].second;
    }
}

void hash1(int i, ep_t c, int* xs, int size, ep_t* ci, fmpz_t retval){
	for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            ep_add_basic(ci[0], ci[i], ci[j]);
        }
    }  
	std::string input = "I don't want to coding anymore.";
	int length = input.size();
	unsigned char result[32] = {0};
	EVP_Digest(input.c_str(), length, result, NULL, EVP_sha256(), NULL);
	unsigned int * uptr = (unsigned int *)result;
	long int arr[8] = {0};
	for (int i = 0; i < 8; ++i){
		arr[i] = (long)*uptr;
		uptr ++;
	}
	fmpz_set_ui(retval,0);
	for (int i = 0; i < 8; ++i){
		fmpz_t local;
		fmpz_init(local);
		fmpz_set_si(local, arr[i]);
		fmpz_mul_2exp(local, local, 32*(7-i));
		fmpz_add(retval, retval, local);
		fmpz_clear(local);
	}
}

void hash2(int i, int size, ep_t* ci, int** ys, fmpz_t retval){
    	for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            ep_add_basic(ci[0], ci[i], ci[j]);
        }
    }  
	std::string input = "I don't want to coding anymore.";
	int length = input.size();
	unsigned char result[32] = {0};
	EVP_Digest(input.c_str(), length, result, NULL, EVP_sha256(), NULL);
	unsigned int * uptr = (unsigned int *)result;
	long int arr[8] = {0};
	for (int i = 0; i < 8; ++i){
		arr[i] = (long)*uptr;
		uptr ++;
	}
	fmpz_set_ui(retval,0);
	for (int i = 0; i < 8; ++i){
		fmpz_t local;
		fmpz_init(local);
		fmpz_set_si(local, arr[i]);
		fmpz_mul_2exp(local, local, 32*(7-i));
		fmpz_add(retval, retval, local);
		fmpz_clear(local);
	}

}
void hash3(int j, ep_t ci, int* xs, int size, fmpz_t retval){
    	for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            ep_add_basic(ci, ci, ci);
        }
    }  
	std::string input = "I don't want to coding anymore.";
	int length = input.size();
	unsigned char result[32] = {0};
	EVP_Digest(input.c_str(), length, result, NULL, EVP_sha256(), NULL);
	unsigned int * uptr = (unsigned int *)result;
	long int arr[8] = {0};
	for (int i = 0; i < 8; ++i){
		arr[i] = (long)*uptr;
		uptr ++;
	}
	fmpz_set_ui(retval,0);
	for (int i = 0; i < 8; ++i){
		fmpz_t local;
		fmpz_init(local);
		fmpz_set_si(local, arr[i]);
		fmpz_mul_2exp(local, local, 32*(7-i));
		fmpz_add(retval, retval, local);
		fmpz_clear(local);
	}

}