#include <iostream>
#include <cstdio>
#include "omp.h"

struct csr{
    int N, NT;
    int* IAT;
    int* JA;
    double* coef;
};

int alloc_csr(csr &A, const int N, const int NT){
    A.N = N; A.NT = NT;
    A.IAT = (int*) malloc( (N+1) * sizeof(int));
    A.JA = (int*) malloc( (NT) * sizeof(int));
    A.coef = (double*) malloc( (NT) * sizeof(double));
    if (A.IAT == nullptr) return 1;
    if (A.JA == nullptr) return 2;
    if (A.coef == nullptr) return 3;
    return 0;
}

int dealloc_csr( csr &A){
    free(A.coef);
    free(A.IAT);
    free(A.JA);
    A.N = 0;
    A.NT=0;
    return 0;
}

csr sum_csr(const csr &A,const csr &B){
    // somma matrici sparse con lo stesso pattern
    csr C;
    int ierr = alloc_csr(C, A.N, A.NT);
    // if (ierr!=0){}
    double* coef_a = A.coef; // riferimenti diretti o handle
    double* coef_b = B.coef;
    double* coef_c = C.coef;
    #pragma omp parallel for num_threads(4) 
    {
        for (int i = 0; i< A.NT; i++){
            coef_c[i] = coef_a[i]+ coef_b[i];
        }
    }
    return C;
}