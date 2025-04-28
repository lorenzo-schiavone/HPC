#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <omp.h>
#include "my_blas.h"

// Compute res = v + u.
void sumvec(double* v, double* u, double* res, int nn) {
    int np = 4;
    int bsize = nn / np;
    #pragma omp parallel num_threads(np)
    {
        int myid = omp_get_thread_num();
        int istart = myid * bsize;
        int iend = (myid + 1) * bsize;
        if (myid == np - 1)
            iend = nn;
        for (int i = istart; i < iend; i++) {
            res[i] = v[i] + u[i];
        }
    }
}

// Compute v = v + a * u.
void daxpy(double* v, double* u, double a, int nn) {
    int np = 4;
    int bsize = nn / np;
    #pragma omp parallel num_threads(np)
    {
        int myid = omp_get_thread_num();
        int istart = myid * bsize;
        int iend = (myid + 1) * bsize;
        if (myid == np - 1)
            iend = nn;
        for (int i = istart; i < iend; i++) {
            v[i] += a * u[i];
        }
    }
}

// Compute the dot (scalar) product of vectors v and u.
double scalarprod(double* v, double* u, int nn) {
    double alpha = 0.0;
    #pragma omp parallel for reduction(+:alpha)
    for (int i = 0; i < nn; i++) {
        alpha += v[i] * u[i];
    }
    return alpha;
}

// Matrix-vector product: res = A * v.
// A has dimensions nr x nc.
void matvecprod(double** A, double* v, double* res, int nc, int nr) {
    for (int i = 0; i < nr; i++) {
        res[i] = 0.0;
        for (int j = 0; j < nc; j++) {
            res[i] += A[i][j] * v[j];
        }
    }
}

// Matrix-matrix product: C = A * B.
// A: m x n, B: n x p, C: m x p.
void matmatprod(double** A, double** B, double** C, int m, int n, int p) {
    #pragma omp parallel for
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < p; j++) {
            C[i][j] = 0.0;
            for (int k = 0; k < n; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

// Transpose the matrix A (nr x nc) into At (nc x nr).
void transpose(double** A, double** At, int nr, int nc) {
    for (int i = 0; i < nr; i++) {
        for (int j = 0; j < nc; j++) {
            At[j][i] = A[i][j];
        }
    }
}

// Sparse matrix-vector product for a matrix stored in CSR format.
void sparsematvecprod(csr& A, double* v, double* res) {
    int* iat = A.iat;
    int* ja = A.ja;
    double* coef = A.coef;
    int nr = A.nr;
    #pragma omp parallel for
    for (int i = 0; i < nr; i++) {
        res[i] = 0.0;
        for (int j = iat[i]; j < iat[i + 1]; j++) {
            res[i] += coef[j] * v[ja[j]];
        }
    }
}
