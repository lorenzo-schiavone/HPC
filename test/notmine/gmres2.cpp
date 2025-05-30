#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

/* Prototipi */
double scalarProd(double* v, double* w, int n, int np);
double norm_vec(double* v, int n, int np);
void matcsrvecprod(int nn,int* iat,int* ja,double* coef,double* v, double* x, int np);
void daxyps(int nn, double* y, double* x, double alpha, int np);

double scalarProd(double* v, double* w, int n, int np){
    double alpha = 0.; 
    #pragma omp parallel num_threads(np) reduction(+:alpha) 
    {
        #pragma omp for schedule(static)
        for (int i=0; i<n; i++){
            alpha+=v[i]*w[i];
        } 
    }
    return alpha;
}

double norm_vec(double* rhs, int n, int np){
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

/**
 * QR Modified Gram‑Schmidt per H ∈ ℝ^{(k+1)×k}
 * - m = k+1 righe, k colonne
 * - H: input [0..m-1][0..k-1]
 * - Q: output [0..m-1][0..k-1]
 * - R: output [0..k-1][0..k-1]
 */
void qr_mgs(int m, int k,
            double** H, double** Q, double** R,
            int np)
{
    // azzero R
    for(int i=0;i<k;i++)
        for(int j=0;j<k;j++)
            R[i][j] = 0.0;

    for(int j=0;j<k;j++){
        // copia H[:,j] in Q[:,j]
        #pragma omp parallel for num_threads(np)
        for(int i=0;i<m;i++)
            Q[i][j] = H[i][j];

        // ortonormalizza rispetto alle colonne precedenti
        for(int i=0;i<j;i++){
            // R[i][j] = Q[:,i]^T * Q[:,j]
            R[i][j] = scalarProd(Q[i], Q[j], m, np);
            // Q[:,j] -= R[i][j]*Q[:,i]
            #pragma omp parallel for num_threads(np)
            for(int l=0;l<m;l++)
                Q[l][j] -= R[i][j] * Q[l][i];
        }

        // diagonale R
        R[j][j] = norm_vec(Q[j], m, np);
        // normalizza Q[:,j]
        #pragma omp parallel for num_threads(np)
        for(int i=0;i<m;i++)
            Q[i][j] /= R[j][j];
    }
}


void gmres(int n, int* iat, int* ja, double* coef,
           double* rhs, double tol, int maxit,
           int np, double* x)
{
    double beta = norm_vec(rhs, n, np);
    double exit_cond = tol * beta;

    double* resvec = (double*)malloc((maxit+1)*sizeof(double));
    resvec[0] = beta;
    for(int i=1;i<=maxit;i++) resvec[i]=0.0;

    // controlla qui
    if (beta < exit_cond){
        free(resvec);
        return;
    }

    double** V    = (double**) malloc((maxit+1)*sizeof(double*));
    double*  Vbuf =(double*)  malloc(n*(maxit+1)*sizeof(double));
    for(int i=0;i<=maxit;i++)
        V[i] = &Vbuf[n*i];

    // H: (maxit+1)×maxit
    double** H    = (double**) malloc((maxit+1)*sizeof(double*));
    double*  Hbuf = (double*) malloc((maxit+1)*maxit*sizeof(double));
    for(int i=0;i<=maxit;i++)
        H[i] = &Hbuf[i*n];

    // Q: (maxit+1)×maxit, R: maxit×maxit
    double** Q    = (double**)malloc((maxit+1)*sizeof(double*));
    double*  Qbuf = (double*)malloc((maxit+1)*maxit*sizeof(double));
    double** R    = (double**)malloc(maxit*sizeof(double*));
    double*  Rbuf = (double*)malloc(maxit*maxit*sizeof(double));
    for(int i=0;i<=maxit;i++) Q[i] = &Qbuf[i*maxit];
    for(int i=0;i< maxit;i++) R[i] = &Rbuf[i*maxit];

    // vnew
    double* vnew = (double*) malloc(n*sizeof(double));

    // V[:,0] = rhs/beta
    #pragma omp parallel for num_threads(np)
    for(int i=0;i<n;i++)
        V[0][i] = rhs[i]/beta;

    int it = 0;
    while(it < maxit && resvec[it] > exit_cond){
        matcsrvecprod(n, iat, ja, coef, V[it], vnew, np);

        for(int j=0;j<=it;j++){
            H[j][it] = scalarProd(vnew, V[j], n, np);
            daxyps(n, vnew, V[j], -H[j][it], np);
        }
        double h_new = norm_vec(vnew, n, np);
        H[it+1][it] = h_new;

        if (h_new < 1e-10){
            printf("happy breakdown at iter %d\n", it);
            qr_mgs(it+1, it+1, H, Q, R, np);
            break;
        }

        #pragma omp parallel for num_threads(np)
        for(int i=0;i<n;i++)
            V[it+1][i] = vnew[i]/h_new;

        qr_mgs(it+2, it+1, H, Q, R, np);

        resvec[it+1] = fabs(beta * Q[it+1][0]);

        it++;
    }

    double* y    = (double*) malloc(it*sizeof(double));
    double* rhs2 = (double*) malloc(it*sizeof(double));
    for(int i=0;i<it;i++)
        rhs2[i] = beta * Q[0][i];
    for(int i=it-1;i>=0;i--){
        double s = 0.;
        for(int j=i+1;j<it;j++)
            s += R[i][j]*y[j];
        y[i] = (rhs2[i] - s) / R[i][i];
    }

    #pragma omp parallel for num_threads(np)
    for(int i=0;i<n;i++){
        double sum = x[i];
        for(int j=0;j<it;j++)
            sum += V[j][i] * y[j];
        x[i] = sum;
    }

    free(resvec);
    free(Vbuf);  free(V);
    free(Hbuf);  free(H);
    free(Qbuf);  free(Q);
    free(Rbuf);  free(R);
    free(vnew);
    free(y);
    free(rhs2);
}
