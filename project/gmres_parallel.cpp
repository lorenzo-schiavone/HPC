#include <stdlib.h>
#include <cstdlib>
#include <fstream>
#include <omp.h>
#include <math.h>

#include "qr.h"

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

void gmres(int nrows, int* iat, int* ja, double* coef, double* rhs, double tol, int maxit, int np, double* x){

    // inizializzazione variabili
    double* diag = (double*) malloc(nrows * sizeof(double));
    double* Mb = (double*) malloc(nrows * sizeof(double));
    double beta = 0.;
    double* resvec = (double*) malloc (maxit*sizeof(double));
    for (int i=0; i<maxit; i++){
        resvec[i]=0.;
    }
    // base arnoldi
    double** V = (double**) malloc(maxit*sizeof(double*));
    double* Vbuff = (double*) malloc(nrows*maxit*sizeof(double));
    for (int i=0;i<maxit; i++){
        V[i] = &Vbuff[nrows*i];
    }
    // hessenberg matrix
    double** H = (double**) malloc(maxit*sizeof(double*));
    double* Hbuff = (double*) malloc(maxit*maxit*sizeof(double));
    for (int i=0;i<maxit; i++){
        H[i] = &Hbuff[maxit*i];
    }
    // temp vector
    double* vnew = (double*) malloc (nrows*sizeof(double));
    double* vold = (double*) malloc (nrows*sizeof(double));

    // puntatori Q,R
    double** Q; double**R;
    double hj = 0.; // per tenere prodotti scalari
    double hnew =0.; // per tenere norme
    int it = 0;
    bool flag = false;
    double exit_cond;
    double* y;

    // inizio regione parallela
    # pragma omp parallel num_threads(np)
    {
        # pragma omp for 
        for (int i = 0; i < nrows; i++) {
            diag[i] = 1.0; x[i]=0.; // initialize sol vector x and diag
            vnew[i] = 0.;
            Mb[i]= rhs[i]; // initialize Mb as rhs, if diag[i]= 1 do nothing
            for (int j = iat[i]; j < iat[i+1]; j++) {
                if (ja[j] == i) {  // if there is a diagonal term put it in diag
                    diag[i] = coef[j];
                    Mb[i]/= diag[i];
                    break;
                }
            }
        }
        // calcolo beta: beta = norm(Mb, nrows, np); 
        #pragma omp for reduction(+:beta)
        for (int i=0; i<nrows; i++){
            beta+=Mb[i]*Mb[i];
        }
        # pragma omp single
        {
            beta = sqrt(beta);
        }
        resvec[0] = beta;
        exit_cond = tol * beta;
        if (beta > tol){
            // first column of the basis
            # pragma omp for
            for (int i=0; i<nrows;i++){
                V[0][i] = Mb[i]/beta;
            }
            # pragma omp single
            {
                vold = V[0];
            }

            while ((resvec[it] > exit_cond) && (it < maxit-1) && (!flag)){
                # pragma omp single
                {
                    it++;
                }
                // vnew = AV[it-1] // vold is supposed to contain V[it-1]
                #pragma omp for 
                for (int i = 0; i < nrows; i++){
                    vnew[i] = 0.0;
                    for (int j = iat[i]; j < iat[i+1]; j++){
                        vnew[i] += coef[j] * vold[ja[j]];
                    }
                    vnew[i] /= diag[i];
                }

                // remove orthogonal component from previous vector basis
                
                for (int j = 0; j< it; j++){
                    vold = V[j];
                    // scalar product
                    hj = 0.;
                    # pragma omp for reduction(+:hj)
                    for (int i=0; i<nrows; i++){
                        hj+=vold[i]*vnew[i];
                    }

                    // daxyps
                    # pragma omp for 
                    for (int i = 0; i < nrows; i++){
                        vnew[i]-= (hj * vold[i]);
                    }
                    H[it-1][j]=hj; 
                    #pragma omp barrier
                }
                
                // calcolo norma di vnew una volta tolte le componenti linearmente dipendenti ai vettori precedenti
                hnew = 0.;
                #pragma omp for reduction(+:hnew)
                for (int i=0; i<nrows; i++){
                    hnew+=vnew[i]*vnew[i];
                }
                # pragma omp single
                {
                    hnew = sqrt(hnew);
                    if (hnew < 1e-10){
                        printf("happy breakdown!\n");
                        qr(H, it, it, &Q, &R);
                        flag = true;
                    }
                }
                if (!flag){
                    # pragma omp for
                    for (int i=0; i<nrows; i++){
                        V[it][i] = vnew[i]/hnew; 
                    }
                    #pragma omp single
                    {
                        H[it-1][it] = hnew;
                        qr(H, it+1, it, &Q, &R);
                        resvec[it] = fabs(beta*Q[it][0]);
                        vold = V[it];
                        printf("resvec[%u]: ", it);
                        printf("%f\n", resvec[it]);
                    }
                }
            }
            # pragma omp single
            {
                double* finalrhs = (double*) malloc ((it-1)*sizeof(double));
                for (int i=0;i<it-1;i++){
                    finalrhs[i] = Q[i][0]*beta; // prima riga * beta ((beta *Q^Te1))
                }

                // I have to solve Ry = beta Q^Te1
                y = (double*) malloc((it-1)*sizeof(double));
                for (int i = it-2; i >= 0; i--) {
                    y[i] = finalrhs[i];
                    for (int j = i+1; j < it-1; j++) {
                        y[i] -= R[j][i] * y[j];
                    }
                    y[i] /= R[i][i];
                }
            }
            
            // sommo su x
            for (int j = 0; j < it-1; j++) {
                vold = V[j];
                # pragma omp for
                for (int i = 0; i < nrows; i++){
                    x[i]-= (y[j] * vold[i]);
                }
            }
        }
    }          
    return;
}
