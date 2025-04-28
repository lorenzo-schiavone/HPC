#include <cstdio>
#include <stdlib.h>
#include <omp.h>
#include "my_blas.h"

void sumvec( double* v, double* u, double* res, int nn){
    // double* res = (double*) malloc (nn*sizeof(double));
    int np=4;
    int bsize = nn / np;

    #pragma omp parallel num_threads(np)
    {
        int myid = omp_get_thread_num();
        int istart = myid*bsize;
        int iend = (myid +1)*bsize; 
        if (myid == np-1){
            iend = nn;
        }
        for (int i=istart; i<iend; i++){
            res[i]= v[i] + u[i];
        } 
    }

}

void daxpy( double* v, double* u, double a, int nn){
    // v = v + a * x
    int np= 4;
    int bsize = nn / np;
    #pragma omp parallel num_threads(np)
    {
        int myid = omp_get_thread_num();
        int istart = myid*bsize;
        int iend = (myid +1)*bsize; 
        if (myid == np-1){
            iend = nn;
        }
        for (int i=istart; i<iend; i++){
            v[i]+= a * u[i];
        } 
    }
}

double scalarprod(double* v, double* u, int nn){
    double alpha = 0.;
    int np=4;
    double* ridv = (double*) malloc(np*sizeof(double)); // uno per processore usato
    int bsize = nn/np;
    #pragma omp parallel num_threads(np)// riduzione scelta da lui. possiamo essere un po' più precisi
    {
        int myid = omp_get_thread_num();
        int istart = myid*bsize;
        int iend = (myid +1)*bsize; 
        if (myid == np-1){
            iend = nn;
        }
        // #pragma omp for schedule(static)
        double alpha_loc=0.;
        for (int i=istart; i<iend; i++){
            alpha_loc+=v[i]*u[i];
        } 
        ridv[myid] = alpha_loc;
        
    }
    for (int j=0; j<np;j++){ 
        alpha+=ridv[j];
    }
    free(ridv);
    return alpha;
}

void matvecprod(double** A, double* v, double* res, int nc, int nr){
    double* Ai;
    for (int i =0; i < nr; i++){
        // crea handle
        Ai = A[i];
        res[i]=0;
        for (int j =0; j< nc; j++){
            res[i] +=Ai[j] * v[j];
        }
    }
}

// void matmatprod(double** A, double** B, double** C, int nr, int nc){
//     double* Ci;
//     double* Ai;
//     for (int i=0;i < nr; i++){
//         Ci = C[i];
//         Ai = A[i];
//         for (int j=0;j< nc; j++){
//             Ci[j] = 0.0; 
//             for (int k=0; k<nc; k++){
//                 Ci[j] += Ai[k] * B[k][j];
//             }
//         }
//     }
// }

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

void transpose(double** A, double** At, int nr, int nc){
    double* Ai;
    for (int i=0;i<nr;i++){
        Ai = A[i];
        for (int j=0;j<nc;j++){
            At[j][i]= Ai[j];
        } 
    }
}

void sparsematvecprod(csr& A, double* v, double* res){
    int *iat = A.iat;
    int*ja = A.ja;
    double *coef = A.coef;
    int nr = A.nr;
    double *x = (double*) malloc(nr*sizeof(double));
    int np=4;
    #pragma omp parallel for num_threads(np)
    for (int i = 0; i < nr; i++){
        res[i] = 0.0;
        for (int j = iat[i]; j < iat[i+1]; j++){
            res[i] += coef[j] * v[ja[j]];
        }
    }
}

