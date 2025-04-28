#ifndef myblas_h
#define myblas_h

#include <cstdio>
#include <stdlib.h>

struct csr{
    int nr;
    int nt;
    int* iat;
    int* ja;
    double* coef;
};

void sumvec( double* v, double* u, double* res, int nn);
void daxpy( double* v, double* u, double a, int nn);
double scalarprod(double* v, double* u, int nn);
void matvecprod(double** A, double* v, double* res, int nc, int nr);
void sparsematvecprod(csr& A, double* v, double* res);
void matmatprod(double** A, double** B, double** C, int m, int n, int p);
void transpose(double** A, double** At, int nr, int nc);

#endif