#include <cstdio>
#include <stdlib.h>
#include <cmath>
#include "my_blas.h"

/*
QR FACTORIZATION FOR SMALL MATRICES: column based
*/
void qr(double** A, double** Q, double** R, int nr, int nc){
    double norm;
    
    // First column
    for (int i = 0; i < nr; i++)
        Q[0][i] = A[0][i];
    
    R[0][0] = sqrt(scalarprod(Q[0], Q[0], nr));
    
    for (int i = 0; i < nr; i++)
        Q[0][i] /= R[0][0];
    
    // Remaining columns
    for (int i = 1; i < nc; i++){
        for (int j = 0; j < nr; j++)
            Q[i][j] = A[i][j];  // copy A[i] into Q[i]

        for (int j = 0; j < i; j++){
            R[j][i] = scalarprod(Q[j], Q[i], nr);
            daxpy(Q[i], Q[j], -R[j][i], nr);  // Q[i] -= R[j][i] * Q[j]
        }

        R[i][i] = sqrt(scalarprod(Q[i], Q[i], nr));
        for (int j = 0; j < nr; j++)
            Q[i][j] /= R[i][i];
    }
}

int main(int argc, const char* argv[]){
    int nr=10;
    int nc=10;
    // ALLOC A
    double** A = (double **) malloc(nr *sizeof(double*)); // one per row
    double* Abuf;
    Abuf = (double*) malloc ( nr*nc * sizeof(double));
    for (int i=0;i<nr;i++){
        A[i] = &Abuf[i*nc];
    }

    // INITIALIZE RANDOMLY
    double* Ai;
    for (int i=0; i<nr; i++){
        Ai = A[i];
        for (int j=0;j<nc; j++){
            Ai[j]=double(rand())/ RAND_MAX;
        }
    } 

    // ALLOC At
    double** At = (double **) malloc(nc *sizeof(double*)); // one per row
    double* Atbuf;
    Atbuf = (double*) malloc ( nr*nc * sizeof(double));
    for (int i=0;i<nc;i++){
        At[i] = &Atbuf[i*nr];
    }

    transpose(A,At, nr, nc);

    // PREALLOC Q R
    double** R = (double **) malloc(nr *sizeof(double*)); // one per row
    double* Rbuf;
    Rbuf = (double*) malloc ( nr*nc * sizeof(double));
    for (int i=0;i<nr;i++){
        R[i] = &Rbuf[i*nc];
    }

    double** Q = (double **) malloc(nr *sizeof(double*)); // one per row
    double* Qbuf;
    Qbuf = (double*) malloc ( nr*nr * sizeof(double));
    for (int i=0;i<nr;i++){
        Q[i] = &Qbuf[i*nr];
    }
    // TEST QR = I and print
    qr(At,Q,R,nr,nc);

    // PREALLOC RES
    double** C = (double **) malloc(nr *sizeof(double*)); // one per row
    double* Cbuf;
    Cbuf = (double*) malloc ( nr*nc * sizeof(double));
    for (int i=0;i<nr;i++){
        C[i] = &Cbuf[i*nc];
    }

    matmatprod(Q,R,C,nr,nc);

    // print C
    // Print difference C - A
    printf("Q * R - A:\n");
    for (int i = 0; i < nr; i++) {
        for (int j = 0; j < nc; j++) {
            printf("% .3e ", C[i][j] - A[j][i]);
        }
        printf("\n");
    }

    
    return 0;
}