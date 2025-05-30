#include <stdlib.h>
#include <cstdlib>
#include <fstream>
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
    double* resvec = (double*)malloc(maxit * sizeof(double));
    for (int i = 0; i < maxit; i++) resvec[i] = 0.;

    double** V = (double**)malloc(maxit * sizeof(double*));
    double* Vbuff = (double*)malloc(nrows * maxit * sizeof(double));
    for (int i = 0; i < maxit; i++) V[i] = &Vbuff[nrows * i];

    double** H = (double**)malloc(maxit * sizeof(double*));
    double* Hbuff = (double*)malloc(maxit * maxit * sizeof(double));
    for (int i = 0; i < maxit; i++) H[i] = &Hbuff[maxit * i];

    double* vnew = (double*)malloc(nrows * sizeof(double));
    double* vold = (double*)malloc(nrows * sizeof(double));
    double* Vj;

    double** Q; double** R;
    double hj = 0.;
    double hnew = 0.;
    int it = 0;
    bool flag = false;
    double exit_cond;
    double* y;

    // Initialize diag, x, Mb
    #pragma omp parallel num_threads(np)
    {
        #pragma omp for
        for (int i = 0; i < nrows; i++) {
            diag[i] = 1.0;
            x[i] = 0.0;
            Mb[i] = rhs[i];
            for (int j = iat[i]; j < iat[i+1]; j++) {
                if (ja[j] == i) {
                    diag[i] = coef[j];
                    Mb[i] /= diag[i];
                    break;
                }
            }
        }

        // Compute beta
        #pragma omp for reduction(+:beta) 
        for (int i = 0; i < nrows; i++) {
            beta += Mb[i] * Mb[i];
            }
        #pragma omp single
        {
            beta = sqrt(beta);
        }
        // Initialize V[0]
        #pragma omp for 
        for (int i = 0; i < nrows; i++){
            V[0][i] = Mb[i] / beta;
            } 
    }
    vold = V[0];
    resvec[0] = beta;
    exit_cond = tol * beta;

    if (beta > tol) {
        while ((it < maxit - 1) && (resvec[it] > exit_cond) && (!flag)) {
            it++;
            // Matrix-vector product: vnew = A * vold / diag
            #pragma omp parallel num_threads(np)
            {
                #pragma omp for
                for (int i = 0; i < nrows; i++) {
                    vnew[i] = 0.0;
                    for (int j = iat[i]; j < iat[i+1]; j++) {
                        vnew[i] += coef[j] * vold[ja[j]];
                    }
                    vnew[i] /= diag[i];
                }

                // Orthogonalization
                for (int j = 0; j < it; j++) {
                    
                    #pragma omp single
                    {
                    Vj = V[j];
                    hj = 0.0;
                    }
                    #pragma omp for reduction(+:hj)
                    for (int i = 0; i < nrows; i++){
                        hj += Vj[i] * vnew[i];
                    } 
                    H[it-1][j] = hj;

                    #pragma omp for
                    for (int i = 0; i < nrows; i++){
                        vnew[i] -= hj * Vj[i];
                    } 
                }

                // Compute hnew
                #pragma omp single
                {
                    hnew = 0.0;
                }
                #pragma omp for reduction(+:hnew)
                for (int i = 0; i < nrows; i++) {
                    hnew += vnew[i] * vnew[i];
                }
                #pragma omp single
                {
                    hnew = sqrt(hnew);
                    if (hnew < 1e-10) {
                    printf("happy breakdown!\n");
                    qr(H, it, it, &Q, &R);
                    flag = true;
                    }
                } 
                if(!flag) {
                    #pragma omp for 
                    for (int i = 0; i < nrows; i++) {
                        V[it][i] = vnew[i] / hnew;
                    }
                    #pragma omp single
                    {
                        H[it-1][it] = hnew;
                        qr(H, it + 1, it, &Q, &R);
                        resvec[it] = fabs(beta * Q[it][0]);
                        vold = V[it]; 
                    }
                }
            }
        }

        // Solve Ry = beta Q^T e1
        double* finalrhs = (double*)malloc(it * sizeof(double));
        for (int i = 0; i < it; i++) finalrhs[i] = beta * Q[i][0];
        y = (double*)malloc(it * sizeof(double));
        for (int i = it - 1; i >= 0; i--) {
            y[i] = finalrhs[i];
            for (int j = i + 1; j < it; j++) y[i] -= R[j][i] * y[j];
            y[i] /= R[i][i];
        }

        // Update solution x
        #pragma omp parallel for num_threads(np)
        for (int i = 0; i < nrows; i++) {
            for (int j = 0; j < it; j++) {
                x[i] += y[j] * V[j][i];
            }
        }
        // free(finalrhs);
        // free(y);
    }

    // free(diag);
    // free(Mb);
    // free(resvec);
    // free(V);
    // free(Vbuff);
    // free(H);
    // free(Hbuff);
    // free(vnew);
    // free(vold);
    return;
}