#include <stdlib.h>
#include <cstdlib>
#include <fstream>
#include <omp.h>
#include <math.h>

#include "qr.h"

void printMatrix(double** mat, int nrows, int ncols, const char* name) {
    printf("%s:\n", name);
    for (int i = 0; i < nrows; ++i) {
        for (int j = 0; j < ncols; ++j) {
            printf("%2.2f ", mat[j][i]);
        }
        printf("\n");
    }
}

double scalarProd(double* v, double* w, int n, int np){
    double alpha = 0.; // alpha inizzializzato con l'elemento neutro dell'operatore
    #pragma omp parallel num_threads(np) reduction(+:alpha) // riduzione scelta da lui. possiamo essere un po' più precisi
    {
        #pragma omp for schedule(static)
        for (int i=0; i<n; i++){
            alpha+=v[i]*w[i];
        } 
    }
    return alpha;
}

double norm(double* rhs, int n, int np){
    return sqrt(scalarProd(rhs, rhs, n, np));
}

void matcsrvecprod(int nn,int* iat,int* ja,double* coef,double* v, double* x, int np){
    #pragma omp parallel for num_threads(np)
    for (int i = 0; i < nn; i++){
        x[i] = 0.0;
        for (int j = iat[i]; j < iat[i+1]; j++){
            x[i] += coef[j] * v[ja[j]];
        }
    }
}

void daxyps(int nn, double* vnew, double* w, double alpha, int np){
    // substitute vnew with vnew + alpha *w 
    #pragma omp parallel for num_threads(np)
    for (int i = 0; i < nn; i++){
        vnew[i]+= alpha * w[i];
    }
}

void gmres(int n, int* iat, int* ja, double* coef, double* rhs, double tol, int maxit, int np, double* x){
    
    // initial estimate is 0
    // double* x = (double*) malloc (n*sizeof(double));
    // for (int i=0; i<n; i++){
    //     x[i]=0.;
    // }
    double beta = norm(rhs, n, np);

    double* resvec = (double*) malloc (maxit*sizeof(double));
    resvec[0] = beta;
    for (int i=1; i<maxit; i++){
        resvec[i]=0.;
    }
    double exit_cond = tol * beta;
    if (resvec[0]< tol){ // i don't like this how should it be?
        return;
    }
    
    double** V = (double**) malloc(maxit*sizeof(double*));
    double* Vbuff = (double*) malloc(n*maxit*sizeof(double));
    int k=0;
    for (int i=0;i<maxit; i++){
        V[i] = &Vbuff[n*k];
        k++;
    }
    // first col initialization
    for (int i=0; i<n;i++){
        V[0][i] = rhs[i]/beta;
    }

    double** H = (double**) malloc(maxit*sizeof(double*));
    double* Hbuff = (double*) malloc(maxit*maxit*sizeof(double));
    // double** Q = (double**) malloc(maxit*sizeof(double*));
    // double* Qbuff = (double*) malloc(maxit*maxit*sizeof(double));
    // double** R = (double**) malloc(maxit*sizeof(double*));
    // double* Rbuff = (double*) malloc(maxit*maxit*sizeof(double));
    for (int i=0;i<maxit; i++){
        H[i] = &Hbuff[maxit*i];
        // Q[i] = &Qbuff[maxit*k];
        // R[i] = &Rbuff[maxit*k];
    }


    double* vnew = (double*) malloc (n*sizeof(double));
    for (int i=0; i<n; i++){
        vnew[i]=0.;
    }

    // printf("v0: \n");
    // for (int i=0; i<n; i++){
    //     printf("%f ", V[0][i]);

    // }
    // printf("\n");
    double** Q; double**R;
    double hj = 0.;
    int it_number = 0;
    printf("exit_cond: %f\n", exit_cond);
    // broken it_number and qr factorization :(
    while ((resvec[it_number] > exit_cond) && (it_number < maxit)){
        it_number ++;
        matcsrvecprod(n, iat, ja, coef, V[it_number-1], vnew, np); // in place on (old) vnew
        
        // printf("vnew: \n");
        // for (int i=0; i<n; i++){
        //     printf("%f ", vnew[i]);

        // }
        // printf("\n");

        for (int j = 0; j< it_number; j++){
            hj = scalarProd(vnew, V[j],n,np);
            daxyps(n, vnew, V[j], - hj, np); // in place for vnew
            H[it_number-1][j]=hj;   
        }

        hj = norm(vnew, n, np);
        if (hj < 1e-10){
            printf("happy breakdown!\n");
            qr(H, it_number-1, it_number-1, &Q, &R);
            break;
        }
        else {// should i put an else?
            for (int i=0; i<n; i++){
                V[it_number][i] = vnew[i]/hj; 
            }
            H[it_number-1][it_number] = hj;
        }
        qr(H, it_number+1, it_number, &Q, &R);
        // it_number ++;
        resvec[it_number] = abs(beta*Q[it_number-1][0]);

        printf("resvec[%u]: ", it_number);
        printf("%f\n", resvec[it_number] );
    }
    printf("exited from while loop!\n");
    printMatrix(Q, it_number,it_number, "Q");
    printMatrix(R, it_number,it_number-1, "R");

    double* finalrhs = (double*) malloc ((it_number-1)*sizeof(double));
    for (int i=0;i<it_number-1;i++){
        finalrhs[i] = Q[i][0]*beta;
    }
    // check again
    double* y = (double*) malloc((it_number-1)*sizeof(double));
    for (int i = it_number-1; i >= 0; i--) {
        y[i] = finalrhs[i];
        for (int j = i+1; j < it_number; j++) {
            y[i] -= R[i][j] * y[j];
        }
        y[i] /= R[i][i];
    }

    for (int j = 0; j < it_number-1; j++) {
        daxyps(n, x, V[j], y[j], np);
    }
    return;
}
