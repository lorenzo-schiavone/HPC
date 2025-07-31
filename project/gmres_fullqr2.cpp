#include <stdlib.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <omp.h>
#include <math.h>

#include "qr.h"

void matcsrvecprod(int nrows, int* iat, int* ja, double* coef, double* v, double* x, int np) {
    #pragma omp parallel for num_threads(np)
    for (int i = 0; i < nrows; i++) {
        x[i] = 0.0;
        for (int j = iat[i]; j < iat[i+1]; j++) {
            x[i] += coef[j] * v[ja[j]];
        }
    }
}

void gmres(int nrows, int* iat, int* ja, double* coef, double* rhs, double tol, int maxit, int np, double* x) {
    double* diag = (double*)malloc(nrows * sizeof(double));
    double* Mb = (double*)malloc(nrows * sizeof(double));
    double beta = 0.;
    
    // Initialize diag, x, Mb
    #pragma omp parallel for reduction(+:beta) num_threads(np)
    for (int i = 0; i < nrows; i++) {
        double d = 1.0;
        x[i] = 0.0;
        
        for (int j = iat[i]; j < iat[i+1]; j++) {
            if (ja[j] == i) {
                d = coef[j];
                break;
            }
        }
        diag[i] = d;
        Mb[i] = rhs[i] / d;
        beta += Mb[i] * Mb[i];
    }

    beta = sqrt(beta);
    if(beta < tol){
        memset(x, 0, nrows * sizeof(double));
        free(diag);
        free(Mb);
        return;
    }

    double* resvec = (double*)calloc(maxit, sizeof(double));
    double** V = (double**)malloc(maxit * sizeof(double*));
    double* Vbuff = (double*)malloc(nrows * maxit * sizeof(double));
    for (int i = 0; i < maxit; i++) V[i] = &Vbuff[nrows * i];

    int nrowH = maxit + 1;  
    int ncolH = maxit;      
    double** H = (double**) malloc(ncolH * sizeof(double*));
    double*  Hbuff = (double*) calloc(nrowH * ncolH, sizeof(double));
    for(int j = 0; j < ncolH; j++){
        H[j] = &Hbuff[j * nrowH];
    }

    double* vnew = (double*)malloc(nrows * sizeof(double));
    double* vold; 
    double* Vj;

    double** Q; double** R;
    double hj = 0.;
    double hnew = 0.;
    double exit_cond;
    int it_exit = 0;
    double* y;
    // Initialize V[0]
    #pragma omp parallel for num_threads(np) 
    for (int i = 0; i < nrows; i++){
        V[0][i] = Mb[i] / beta;
    } 
    vold = V[0];
    resvec[0] = beta;
    exit_cond = tol * beta;
    for(int it = 1; it <= maxit; it++){
        // Matrix-vector product: vnew = A * vold / diag
        
        #pragma omp parallel for num_threads(np)
        for (int i = 0; i < nrows; i++) {
            vnew[i] = 0.0;
            for (int j = iat[i]; j < iat[i+1]; j++) {
                vnew[i] += coef[j] * vold[ja[j]];
            }
            vnew[i] /= diag[i];
        }

        // Orthogonalization
        for (int j = 0; j < it; j++) {
            Vj = V[j];
            hj = 0.0;
            #pragma omp parallel for reduction(+:hj) num_threads(np)
            for (int i = 0; i < nrows; i++){
                hj += Vj[i] * vnew[i];
            } 
            H[it-1][j] = hj;
            
            #pragma omp parallel for num_threads(np)
            for (int i = 0; i < nrows; i++){
                vnew[i] -= hj * Vj[i];
            } 
        }

        // Compute hnew
        hnew = 0.0;
        #pragma omp parallel for reduction(+:hnew) num_threads(np)
        for (int i = 0; i < nrows; i++) {
            hnew += vnew[i] * vnew[i];
        }
        hnew = sqrt(hnew);
        if (hnew < 1e-15) {
            printf("happy breakdown!\n");
            qr(H, it, it, &Q, &R);
            it_exit = it;
            break;
        }
        #pragma omp parallel for num_threads(np)
        for (int i = 0; i < nrows; i++) {
            V[it][i] = vnew[i] / hnew;
        }

        H[it-1][it] = hnew;
        qr(H, it + 1, it, &Q, &R);
        resvec[it] = fabs(beta * Q[it][0]);
        if(resvec[it] <= exit_cond){
            it_exit = it;
            break;
        }
        vold = V[it]; 
    }
    if(it_exit == 0){it_exit = maxit; // no convergence
        printf("No convergence: res= %e\n", resvec[it_exit]);
    }
    else{
        printf("Number of iterations: %u\n", it_exit);
    }

    // Solve Ry = beta Q^T e1
    double* finalrhs = (double*)malloc(it_exit * sizeof(double));
    for (int i = 0; i < it_exit; i++) finalrhs[i] = beta * Q[i][0];
    y = (double*)malloc(it_exit * sizeof(double));
    for (int i = it_exit - 1; i >= 0; i--) {
        y[i] = finalrhs[i];
        for (int j = i + 1; j < it_exit; j++) y[i] -= R[j][i] * y[j];
        y[i] /= R[i][i];
    }

    // Update solution x
    #pragma omp parallel for num_threads(np)
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < it_exit; j++) {
            x[i] += y[j] * V[j][i];
        }
    }

    // printf("resvec: \n");
    // for (int i= 0; i<it; i++){
    //     printf("%e ", resvec[i]);
    // }
    // printf("\n");
    free(finalrhs);
    free(y);
    free(diag);
    free(Mb);
    free(resvec);
    free(V);
    free(Vbuff);
    free(H);
    free(Hbuff);
    free(vnew);
}