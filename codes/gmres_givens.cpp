#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <omp.h>

using namespace std;

// x = A *v
void matcsrvecprod(int n, int *iat, int *ja, double *coef, double *v, double *x, int num_threads)
{
    #pragma omp parallel for num_threads(num_threads)
    for(int i = 0; i < n; i++){
        double sum = 0.0;
        for(int k = iat[i]; k < iat[i+1]; k++)
            sum += coef[k] * v[ja[k]];
        x[i] = sum;
    }
}

// flat H: H_ij = H_at(H, ncolH, i, j ) - reference so that we can read and write
inline double &H_at(double *H_flat, int ncolH, int i, int j){
    return H_flat[i * ncolH + j];
}

// GMRES
// sol is preallocated length n
void gmres(int n, int *iat, int *ja, double *coef, double *rhs, double tol, int maxit, int np, double *sol)
{
    // Mb = diag(A)^{-1} rhs, beta = ||Mb||
    double *diag = (double*)malloc(n * sizeof(double));
    double *Mb   = (double*)malloc(n * sizeof(double));
    double  beta = 0.0;
    #pragma omp parallel for reduction(+:beta) num_threads(np)
    for(int i = 0; i < n; i++){
        // diagonal entry
        double d = 1.0; // if not found set to 1 to avoid problem
        for(int k = iat[i]; k < iat[i+1]; k++){
            if(ja[k] == i){
                d = coef[k];
                break;
            }
        }
        diag[i] = d;
        Mb[i]   = rhs[i] / d;
        beta   += Mb[i] * Mb[i];
    }
    beta = sqrt(beta);

    // zero solution
    if(beta < tol){
        memset(sol, 0, n * sizeof(double));
        free(diag);
        free(Mb);
        return;
    }

    // V for krylov bases
    int mmax = maxit + 1;
    double **V    = (double**) malloc(mmax * sizeof(double*));
    double *Vbuff = (double*) malloc(n * mmax * sizeof(double));
    for(int j = 0; j < mmax; j++)
        V[j] = &Vbuff[j * n];

    // H flat 
    int nrowH = maxit + 1, ncolH = maxit;
    double *H   = (double*) calloc(nrowH * ncolH, sizeof(double));
    // Givens vectors
    double *cs     = (double* )malloc(maxit * sizeof(double));
    double *sn     = (double*) malloc(maxit * sizeof(double));
    double *s      = (double*) calloc(maxit+1, sizeof(double));
    double *resid  = (double*)calloc(maxit+1, sizeof(double));

    double *vnew = (double*)malloc(n * sizeof(double));

    // start initialization
    s[0]       = beta;
    resid[0]   = beta;
    double tb  = tol * beta;
    // V0 = Mb / beta
    #pragma omp parallel for num_threads(np)
    for(int i = 0; i < n; i++) 
        V[0][i] = Mb[i] / beta;

    int it_exit = 0; // number of necessary iterations

    // main loop
    for(int it = 1; it <= maxit; it++){
        int j = it - 1;

        // vnew = A V[j], -> M^-1 vnew 
        #pragma omp parallel for num_threads(np)
        for(int i = 0; i < n; i++){
            double sum = 0.0;
            for(int k = iat[i]; k < iat[i+1]; k++){
                sum += coef[k] * V[j][ja[k]];
            }
            vnew[i] = sum / diag[i]; // M^-1 vnew 
        }

        // arnoldi orthonormalize vnew against previous vectors
        for(int k = 0; k <= j; k++){
            double dot = 0.0;
            #pragma omp parallel for reduction(+:dot)
            for(int i = 0; i < n; i++)
                dot += V[k][i] * vnew[i];
            H_at(H,ncolH, k,   j) = dot;      // H(k,j)
            #pragma omp parallel for num_threads(np)
            for(int i = 0; i < n; i++)
                vnew[i] -= dot * V[k][i];
        }
        // norm vnew
        double vnorm = 0.0;
        #pragma omp parallel for reduction(+:vnorm)
        for(int i = 0; i < n; i++){
            vnorm += vnew[i]*vnew[i];
        }
        vnorm = sqrt(vnorm);
        H_at(H,ncolH, j+1, j) = vnorm;       // H(j+1, j)

        if(vnorm < 1e-15){
            it_exit = it; //happy breakdown
            break;
        }
        #pragma omp parallel for num_threads(np)
        for(int i = 0; i < n; i++)
            V[it][i] = vnew[i] / vnorm;

        // apply previous givens rotations to H(0:j+1, j)
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

        resid[it] = fabs(s[j+1]);

        // exit check
        if(resid[it] <= tb){
            it_exit = it;
            break;
        }
    }
    if(it_exit == 0) it_exit = maxit; // no convergence

    printf("Number of iterations: %u\n", it_exit);
    // R y = s[:it_exit-1]  (R is [saved in] H(0:it_exit-1,0:it_exit-1))
    double *y = (double*) malloc(it_exit * sizeof(double));
    for(int i = it_exit-1; i >= 0; i--){
        double sum = 0.0;
        for(int k = i+1; k < it_exit; k++)
            sum += H_at(H,ncolH, i, k) * y[k];
        y[i] = (s[i] - sum) / H_at(H,ncolH, i, i);
    }

    // sol += V[:,0:it_exit-1] * y
    memset(sol, 0, n * sizeof(double));
    #pragma omp parallel for num_threads(np)
    for(int i = 0; i < n; i++){
        double xi = 0.0;
        for(int k = 0; k < it_exit; k++)
            xi += V[k][i] * y[k];
        sol[i] = xi;
    }

    free(diag);
    free(Mb);
    free(Vbuff);
    free(V);
    free(H);
    free(cs);
    free(sn);
    free(s);
    free(resid);
    free(vnew);
    free(y);
}
