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

// flat H: H_ij = H_at(H, ncolH, i, j ) - reference so that we can read and write
inline double &H_at(double *H_flat, int ncolH, int i, int j){
    return H_flat[i * ncolH + j];
}

// GMRES with restart: solve A x = rhs, tol, maxit numero di outer cycle, restart dopo restart iterazioni
void gmres(int nrows, int* iat, int* ja, double* coef,
           double* rhs, double tol, int maxit, int restart,
           int np, double* x) {
    // alloc
    double *diag = (double *)malloc(nrows * sizeof(double));
    double *Mb   = (double *)malloc(nrows * sizeof(double));
    double *r    = (double *)malloc(nrows * sizeof(double));
    // Krylov Basis
    double **V    = (double **)malloc((restart+1) * sizeof(double*));
    double *Vbuff = (double *) malloc(nrows * (restart+1) * sizeof(double));
    for(int i = 0; i < restart+1; i++) V[i] = &Vbuff[nrows*i];
    // flat H 
    int nrowH = restart + 1, ncolH = restart;
    double *H   = (double*) calloc(nrowH * ncolH, sizeof(double));
    
    double *vnew = (double *) malloc(nrows * sizeof(double));
    // givens rotation
    double *cs     = (double*) malloc(restart * sizeof(double));
    double *sn     = (double*) malloc(restart * sizeof(double));
    double *s      = (double*) calloc(restart+1, sizeof(double));
    
    // ouptut info
    int total_it = 0;
    int inner_it = 0;
    double exit_tol;
    double *resvec = (double *) malloc((restart+1) * sizeof(double));
    
    // set x to 0
    memset(x,0,nrows * sizeof(double)); // if we want an initial guess it is enough to comment this

    #pragma omp parallel for num_threads(np)
    for (int i = 0; i < nrows; i++) {
        double d = 1.0;
        for (int j = iat[i]; j < iat[i+1]; j++)
            if (ja[j] == i) { d = coef[j]; break; }
        diag[i] = d;
    }

    while (total_it < maxit) {
        // rhs
        memset(s, 0, (restart+1)*sizeof(double));
        memset(cs, 0, (restart)*sizeof(double));
        memset(sn, 0, (restart)*sizeof(double));        

        matcsrvecprod(nrows, iat, ja, coef, x, r, np);
        #pragma omp parallel for num_threads(np)
        for (int i = 0; i < nrows; i++)
            Mb[i] = (rhs[i] - r[i]) / diag[i];

        double beta = 0.0;
        #pragma omp parallel for reduction(+:beta) num_threads(np)
        for (int i = 0; i < nrows; i++) beta += Mb[i]*Mb[i];

        beta = sqrt(beta);
        if (!(total_it)){ //first iteration
            exit_tol = tol * beta;
            if (beta < tol) {
                break;
            }
        }
        #pragma omp parallel for num_threads(np)
        for (int i = 0; i < nrows; i++)
            V[0][i] = Mb[i] / beta;

        s[0]      = beta;
        resvec[0] = beta;
        int it_exit = 0;
    
        for (int it = 1; it <= restart; it++) {
            int j = it - 1;
            // A * V[it-1] -> vnew
            #pragma omp parallel for num_threads(np)
            for(int i = 0; i < nrows; i++){
                double sum = 0.0;
                for(int k = iat[i]; k < iat[i+1]; k++){
                    sum += coef[k] * V[j][ja[k]];
                }
                vnew[i] = sum / diag[i]; // M^-1 vnew 
            }

            for(int k = 0; k <= j; k++){
                double dot = 0.0;
                #pragma omp parallel for reduction(+:dot) num_threads(np)
                for(int i = 0; i < nrows; i++)
                    dot += V[k][i] * vnew[i];
                H_at(H,ncolH, k,   j) = dot;      // H(k,j)
                #pragma omp parallel for num_threads(np)
                for(int i = 0; i < nrows; i++)
                    vnew[i] -= dot * V[k][i];
            }
            // norm vnew
            double vnorm = 0.0;
            #pragma omp parallel for reduction(+:vnorm) num_threads(np)
            for(int i = 0; i < nrows; i++){
                vnorm += vnew[i]*vnew[i];
            }
            vnorm = sqrt(vnorm);
            H_at(H,ncolH, j+1, j) = vnorm;       // H(j+1, j)

            if(vnorm < 1e-15){
                it_exit = it; //happy breakdown
                break;
            }
            #pragma omp parallel for num_threads(np)
            for (int i = 0; i < nrows; i++)
                V[it][i] = vnew[i] / vnorm;

            // previous givens rotations to H(0:j+1, j)
            for(int k = 0; k < j; k++){
                double h_kj   = H_at(H,ncolH, k,   j);
                double h_k1j  = H_at(H,ncolH, k+1, j);
                double t      = cs[k]*h_kj + sn[k]*h_k1j;
                double t1     = -sn[k]*h_kj + cs[k]*h_k1j;
                H_at(H,ncolH, k,   j) = t;
                H_at(H,ncolH, k+1, j) = t1;
            }
            // new cs sn to zero out H(j+1,j)
            double hjj  = H_at(H,ncolH, j,   j),
                hj1j = H_at(H,ncolH, j+1, j),
                rho  = sqrt(hjj*hjj + hj1j*hj1j); 
            cs[j] =  hjj/rho;
            sn[j] =  hj1j/rho;
            H_at(H,ncolH, j,   j) = rho;
            H_at(H,ncolH, j+1, j) = 0.0;

            // G_j^T ( G_{j-1}^T ....)e_1
            double s_j   = s[j];
            double s_j1  = s[j+1];
            s[j]   =  cs[j]*s_j + sn[j]*s_j1;
            s[j+1] = -sn[j]*s_j + cs[j]*s_j1;

            resvec[it] = fabs(s[j+1]);

            // exit check
            if(resvec[it] <= exit_tol){
                it_exit = it;
                break;
            }
        }
        if (it_exit == 0) it_exit = restart; // no convergence

        double *y = (double *)malloc(it_exit * sizeof(double));
        for (int i = it_exit-1; i >= 0; i--) {
            double sum = 0.0;
            for(int k = i+1; k < it_exit; k++)
                sum += H_at(H,ncolH, i, k) * y[k];
            y[i] = (s[i] - sum) / H_at(H,ncolH, i, i);
        }

        // x += V[:,0:it_exit] * y
        #pragma omp parallel for num_threads(np)
        for (int i = 0; i < nrows; i++) {
            double sum = 0.0;
            for (int j = 0; j < it_exit; j++)
                sum += V[j][i] * y[j];
            x[i] += sum;
        }

        free(y);

        total_it++;
        inner_it=it_exit;
        if (it_exit < restart){
            break;  // convergenza interna
        }
    }

    printf("Total Iteration: %u [cycle %u]\n", inner_it, total_it);
    printf("Final Residue: %1.4e\n", resvec[inner_it]);
    
    free(diag);
    free(Mb);
    free(r);
    free(Vbuff);
    free(V);
    free(H);
    free(resvec);
    free(vnew);
    free(s);
    free(sn);
    free(cs);
}



