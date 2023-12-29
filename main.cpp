#include <iostream>
#include <math.h>
#include <time.h>
#include "Vc.h"
#include "Proof.h"
#include "util.h"

int main(){

    // Set-up
    clock_t start, end;
    start = clock();
    Vc scheme(6,6);
    end = clock();
    std::cout<<"Set-up time = "<<double(end-start)/CLOCKS_PER_SEC<<"s"<<std::endl;

    // Generate commitment
    start = clock();
    scheme.generate_commitment();
    end = clock();
    std::cout<<"Commit time = "<<double(end-start)/CLOCKS_PER_SEC<<"s"<<std::endl;

    // Open 32 positions
    int size = 32;
    int* rows = new int[size];
    int* cols = new int[size];
    select_points(pow(2,6),size, rows, cols);   
    start = clock();
    Proof pi = scheme.open_multiple(rows, cols, size);
    end = clock();
    std::cout<<"Open time = "<<double(end-start)/CLOCKS_PER_SEC<<"s"<<std::endl;


    // Verify the proof
    start = clock();
    scheme.verify(&pi, rows, cols, size);
    end = clock();
    std::cout<<"Verify time = "<<double(end-start)/CLOCKS_PER_SEC<<"s"<<std::endl;
    std::cout << "Test finish." << std::endl;
    return 0;
}