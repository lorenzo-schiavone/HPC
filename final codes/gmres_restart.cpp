// gmres_restart.c
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "qr.h"

// prodotto matrice‑vettore CSR: x = A * v
void matcsrvecprod(int nrows, int* iat, int* ja, double* coef,
                   double* v, double* x, int np) {
    #pragma omp parallel for num_threads(np)
    for (int i = 0; i < nrows; i++) {
        double sum = 0.0;
        for (int j = iat[i]; j < iat[i+1]; j++) {
            sum += coef[j] * v[ja[j]];
        }
        x[i] = sum;
    }
}

// GMRES with restart: solve A x = rhs, tol, maxit totale, restart dopo m iterazioni
void gmres(int nrows, int* iat, int* ja, double* coef,
           double* rhs, double tol, int maxit, int restart,
           int np, double* x) {
    // allocazioni preliminari
    double *diag = (double *)malloc(nrows * sizeof(double));
    double *Mb   = (double *)malloc(nrows * sizeof(double));
    double *r    = (double *)malloc(nrows * sizeof(double));
    // spazio per il sottospazio di Krylov di dimensione restart
    double **V    = (double **)malloc((restart+1) * sizeof(double*));
    double *Vbuff = (double *) malloc(nrows * (restart+1) * sizeof(double));
    for(int i = 0; i < restart+1; i++) V[i] = &Vbuff[nrows*i];
    // matrice di Hessenberg (restart+1)×restart
    double **H    = (double **)malloc(restart * sizeof(double*));
    double *Hbuff = (double *)calloc((restart+1)*restart, sizeof(double));
    for(int j = 0; j < restart; j++) H[j] = &Hbuff[j*(restart+1)];
    double *resvec = (double *) malloc((restart+1) * sizeof(double));
    // buffer per vnew e soluzione intermedia
    double *vnew = (double *) malloc(nrows * sizeof(double));
    // variabili per QR
    double **Q = NULL, **R = NULL;
    int total_it = 0;
    int inner_it = 0;
    double exit_tol;

    // inizializza x a zero
    #pragma omp parallel for num_threads(np)
    for (int i = 0; i < nrows; i++) x[i] = 0.0;

    // precondizionatore diagonale: Mb = M⁻¹ rhs
    #pragma omp parallel for num_threads(np)
    for (int i = 0; i < nrows; i++) {
        double d = 1.0;
        for (int j = iat[i]; j < iat[i+1]; j++)
            if (ja[j] == i) { d = coef[j]; break; }
        diag[i] = d;
    }

    while (total_it < maxit) {
        // printf("Outer cycle %u\n", total_it);
        // --- 1) calcola r = rhs - A*x, applica precondizionatore M⁻¹
        matcsrvecprod(nrows, iat, ja, coef, x, r, np);
        #pragma omp parallel for num_threads(np)
        for (int i = 0; i < nrows; i++)
            Mb[i] = (rhs[i] - r[i]) / diag[i];

        // --- 2) norma beta = ||Mb||
        double beta = 0.0;
        #pragma omp parallel for reduction(+:beta) num_threads(np)
        for (int i = 0; i < nrows; i++) beta += Mb[i]*Mb[i];
        beta = sqrt(beta);
        if (!(total_it)){
            exit_tol = tol * beta;
            if (beta < tol) {
                break;
            }
        }
        // --- 3) inizializza V[0]
        #pragma omp parallel for num_threads(np)
        for (int i = 0; i < nrows; i++)
            V[0][i] = Mb[i] / beta;

        resvec[0] = beta;
        int it_exit = 0;
    
        // --- 4) ciclo interno GMRES fino a restart
        for (int it = 1; it <= restart; it++) {
            // 4.1 A * V[it-1] → vnew
            matcsrvecprod(nrows, iat, ja, coef, V[it-1], vnew, np);
            #pragma omp parallel for num_threads(np)
            for (int i = 0; i < nrows; i++)
                vnew[i] /= diag[i];

            // 4.2 Ortogonalizzazione contro V[0..it-1]
            for (int j = 0; j < it; j++) {
                double hj = 0.0;
                #pragma omp parallel for reduction(+:hj) num_threads(np)
                for (int i = 0; i < nrows; i++)
                    hj += V[j][i] * vnew[i];
                H[it-1][j] = hj;
                #pragma omp parallel for num_threads(np)
                for (int i = 0; i < nrows; i++)
                    vnew[i] -= hj * V[j][i];
            }

            // 4.3 nuovo h
            double hnew = 0.0;
            #pragma omp parallel for reduction(+:hnew) num_threads(np)
            for (int i = 0; i < nrows; i++) hnew += vnew[i]*vnew[i];
            hnew = sqrt(hnew);

            if (hnew < 1e-15) {
                it_exit = it;
                break;  // happy breakdown
            }
            H[it-1][it] = hnew;
            #pragma omp parallel for num_threads(np)
            for (int i = 0; i < nrows; i++)
                V[it][i] = vnew[i] / hnew;

            // 4.4 QR e residuo
            qr(H, it+1, it, &Q, &R);
            resvec[it] = fabs(beta * Q[it][0]);
            if (resvec[it] <= exit_tol) {
                it_exit = it;
                break;
            }
        }
        if (it_exit == 0) it_exit = restart; // no convergence

        // --- 5) risolvi R y = beta * Q^T e1
        double *finalrhs = (double *)malloc(it_exit * sizeof(double));
        for (int i = 0; i < it_exit; i++) finalrhs[i] = beta * Q[i][0];
        double *y = (double *)malloc(it_exit * sizeof(double));
        for (int i = it_exit-1; i >= 0; i--) {
            y[i] = finalrhs[i];
            for (int j = i+1; j < it_exit; j++)
                y[i] -= R[j][i] * y[j];
            y[i] /= R[i][i];
        }

        // --- 6) aggiorna soluzione x += V[:,0:it_exit] * y
        #pragma omp parallel for num_threads(np)
        for (int i = 0; i < nrows; i++) {
            double sum = 0.0;
            for (int j = 0; j < it_exit; j++)
                sum += V[j][i] * y[j];
            x[i] += sum;
        }

        free(finalrhs);
        free(y);

        total_it++;
        inner_it=it_exit;
        if (it_exit < restart){
            break;  // convergenza interna
        }
    }

    printf("Total Iteration: %u [cycle %u]\n", inner_it, total_it);
    printf("Final Residue: %1.4e\n", resvec[inner_it]);
    // libera memoria
    free(diag);
    free(Mb);
    free(r);
    free(Vbuff);
    free(V);
    free(Hbuff);
    free(H);
    free(resvec);
}



