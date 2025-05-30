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

double scalarProd(double* v, double* w, int nrows, int np){
    double alpha = 0.; // alpha inizzializzato con l'elemento neutro dell'operatore
    #pragma omp parallel num_threads(np) reduction(+:alpha) // riduzione scelta da lui. possiamo essere un po' più precisi
    {
        #pragma omp for schedule(static)
        for (int i=0; i<nrows; i++){
            alpha+=v[i]*w[i];
        } 
    }
    return alpha;
}

double norm(double* v, int nrows, int np){
    return sqrt(scalarProd(v, v, nrows, np));
}

// coef*v -> x
void matcsrvecprod(int nrows,int* iat,int* ja,double* coef,double* v, double* x, int np){
    #pragma omp parallel for num_threads(np)
    for (int i = 0; i < nrows; i++){
        x[i] = 0.0;
        for (int j = iat[i]; j < iat[i+1]; j++){
            x[i] += coef[j] * v[ja[j]];
        }
    }
}

// substitute v with v + alpha * w 
void daxyps(int nrows, double* v, double* w, double alpha, int np){
    #pragma omp parallel for num_threads(np)
    for (int i = 0; i < nrows; i++){
        v[i]+= alpha * w[i];
    }
}

void gmres(int nrows, int* iat, int* ja, double* coef, double* rhs, double tol, int maxit, int np, double* x){
    
    // initial estimate is 0
    // x = (double*) malloc (nrows*sizeof(double));
    for (int i=0; i<nrows; i++){
        x[i]=0.;
    }

    double* diag = (double*) malloc(nrows * sizeof(double));
    # pragma omp parallel for num_threads(np)
    for (int i = 0; i < nrows; i++) {
        diag[i] = 1.0;
        for (int j = iat[i]; j < iat[i+1]; j++) {
            if (ja[j] == i) { 
                diag[i] = coef[j];
                break;
            }
        }
    }
    double* Mb = (double*) malloc(nrows * sizeof(double));
    #pragma parallel for num_threads(np)
    for (int i = 0; i < nrows; i++) {
        Mb[i] = rhs[i] / diag[i];
    }
    double beta = norm(Mb, nrows, np); 
    
    // double beta = norm(Mb, nrows, np);

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
    double* Vbuff = (double*) malloc(nrows*maxit*sizeof(double));
    for (int i=0;i<maxit; i++){
        V[i] = &Vbuff[nrows*i];
    }
    // first col initialization
    for (int i=0; i<nrows;i++){
        V[0][i] = Mb[i]/beta;
    }

    double** H = (double**) malloc(maxit*sizeof(double*));
    double* Hbuff = (double*) malloc(maxit*maxit*sizeof(double));
    for (int i=0;i<maxit; i++){
        H[i] = &Hbuff[maxit*i];
    }


    double* vnew = (double*) malloc (nrows*sizeof(double));
    for (int i=0; i<nrows; i++){
        vnew[i]=0.;
    }

    double** Q; double**R;
    double hj = 0.;
    double hnew =0.;
    int it = 0;
    bool flag = false;
    // printf("exit_cond: %f\n", exit_cond);
    while ((resvec[it] > exit_cond) && (it < maxit-1)){
        it ++;
        matcsrvecprod(nrows, iat, ja, coef, V[it-1], vnew, np); // in place on (old) vnew
        #pragma omp parallel for num_threads(np)
        for (int i = 0; i < nrows; i++) {
            vnew[i] /= diag[i]; // vnew = M^{-1} * A * V[it-1]
        }

        for (int j = 0; j< it; j++){
            hj = scalarProd(vnew, V[j],nrows,np);
            daxyps(nrows, vnew, V[j], - hj, np); // in place for vnew
            H[it-1][j]=hj;   
        }

        hnew = norm(vnew, nrows, np);
        if (hnew < 1e-10){
            printf("happy breakdown!\n");
            qr(H, it, it, &Q, &R);
            flag = true;
            break;
        }
        else {// should i put an else?
            for (int i=0; i<nrows; i++){
                V[it][i] = vnew[i]/hnew; 
            }
            H[it-1][it] = hnew;
        
            qr(H, it+1, it, &Q, &R);
            // it ++;
            resvec[it] = fabs(beta*Q[it][0]);

            // printf("resvec[%u]: ", it);
            // printf("%f\n", resvec[it] );
        }
    }
    // printf("exited from while loop!\n");
    if (flag){
        // printMatrix(Q, it-1,it-1, "Q");
        // printMatrix(R, it-1,it-1, "R");
    } else{
        // printMatrix(Q, it,it, "Q");
        // printMatrix(R, it,it-1, "R");
        }
        
        double* finalrhs = (double*) malloc ((it-1)*sizeof(double));
        for (int i=0;i<it-1;i++){
            finalrhs[i] = Q[i][0]*beta; // prima riga * beta ((beta *Q^Te1))
        }
        // I have to solve Ry = beta Q^Te1

        double* y = (double*) malloc((it-1)*sizeof(double));
        for (int i = it-2; i >= 0; i--) {
            y[i] = finalrhs[i];
            for (int j = i+1; j < it-1; j++) {
                y[i] -= R[j][i] * y[j];
            }
            y[i] /= R[i][i];
        }
        for (int j = 0; j < it-1; j++) {
            daxyps(nrows, x, V[j], y[j], np);
        }
    return;
}
