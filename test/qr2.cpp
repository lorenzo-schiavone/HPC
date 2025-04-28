#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "my_blas.h"

void qr(double** A, double** Q, double** R, int m, int n) {
    // QR factorization of A^T (which is m x n matrix)
    for (int i = 0; i < m; i++) {
        Q[0][i] = A[0][i]; // First column of Q = first row of A
    }
    
    R[0][0] = sqrt(scalarprod(Q[0], Q[0], m));
    for (int i = 0; i < m; i++) {
        Q[0][i] /= R[0][0];
    }
    
    for (int i = 1; i < n; i++) { // Process each column of A^T
        for (int j = 0; j < m; j++) {
            Q[i][j] = A[i][j]; // Copy i-th row of A
        }
        
        for (int j = 0; j < i; j++) {
            R[j][i] = scalarprod(Q[j], Q[i], m);
            daxpy(Q[i], Q[j], -R[j][i], m);
        }
        
        R[i][i] = sqrt(scalarprod(Q[i], Q[i], m));
        for (int j = 0; j < m; j++) {
            Q[i][j] /= R[i][i];
        }
    }
}

int main() {
    int nr = 10, nc = 10; // Original A dimensions
    
    // Allocate A (nr x nc)
    double** A = new double*[nr];
    A[0] = new double[nr*nc];
    for (int i = 1; i < nr; i++) A[i] = A[i-1] + nc;
    
    // Initialize A randomly
    for (int i = 0; i < nr; i++) 
        for (int j = 0; j < nc; j++) 
            A[i][j] = (double)rand()/RAND_MAX;
    
    // Allocate Q (nc x nc)
    double** Q = new double*[nc];
    Q[0] = new double[nc*nc];
    for (int j = 1; j < nc; j++) Q[j] = Q[j-1] + nc;
    
    // Allocate R (nc x nr)
    double** R = new double*[nc];
    R[0] = new double[nc*nr];
    for (int j = 1; j < nc; j++) R[j] = R[j-1] + nr;
    
    qr(A, Q, R, nc, nr); // Pass dimensions of A^T
    
    // Verify Q*R - A^T
    double** C = new double*[nc];
    C[0] = new double[nc*nr];
    for (int j = 1; j < nc; j++) C[j] = C[j-1] + nr;
    
    matmatprod(Q, R, C, nc, nc, nr);
    
    printf("Q*R - A^T:\n");
    for (int i = 0; i < nc; i++) {
        for (int j = 0; j < nr; j++) {
            printf("% .2e ", C[i][j] - A[i][j]);
        }
        printf("\n");
    }
    
    // Cleanup
    delete[] A[0]; delete[] A;
    delete[] Q[0]; delete[] Q;
    delete[] R[0]; delete[] R;
    delete[] C[0]; delete[] C;
    
    return 0;
}