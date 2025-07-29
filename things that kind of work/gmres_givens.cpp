// gmres_givens.cpp
// GMRES with Jacobi preconditioning & Givens‐rotation QR updates

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <omp.h>

// Sparse CSR mat‑vec:  x = A * v
void matcsrvecprod(int n,
                   int *iat, int *ja,
                   double *coef,
                   double *v,
                   double *x,
                   int num_threads)
{
    #pragma omp parallel for num_threads(num_threads)
    for(int i = 0; i < n; i++){
        double sum = 0.0;
        for(int k = iat[i]; k < iat[i+1]; k++)
            sum += coef[k] * v[ja[k]];
        x[i] = sum;
    }
}

// index into row‑major H of size (nrowH)×(ncolH):
//    H_ij = H_flat[ i * ncolH + j ]
inline double &H_at(double *H_flat, int ncolH, int i, int j){
    return H_flat[(size_t)i * ncolH + j];
}

// GMRES solver
// A in CSR (iat,ja,coef), solve A x = rhs
// tol,    maxit,
// np = # threads,
// x_out must be allocated length n
void gmres(int n,
           int *iat, int *ja,
           double *coef,
           double *rhs,
           double tol, int maxit,
           int np,
           double *x_out)
{
    // 1) Jacobi precondition: Mb = D^{-1} rhs, beta = ||Mb||
    double *diag = (double*)std::malloc(n * sizeof(double));
    double *Mb   = (double*)std::malloc(n * sizeof(double));
    double  beta = 0.0;
    #pragma omp parallel for reduction(+:beta) num_threads(np)
    for(int i = 0; i < n; i++){
        // find diagonal entry A[i,i]
        double d = 1.0;
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
    beta = std::sqrt(beta);

    // trivial case
    if(beta < tol){
        std::memset(x_out, 0, n * sizeof(double));
        std::free(diag);
        std::free(Mb);
        return;
    }

    // 2) allocate Krylov basis V[:,0..maxit], each length n
    int mmax = maxit + 1;
    double **V    = (double**)std::malloc(mmax * sizeof(double*));
    double *Vbuff = (double*)std::malloc((size_t)n * mmax * sizeof(double));
    for(int j = 0; j < mmax; j++)
        V[j] = &Vbuff[(size_t)j * n];

    // 3) allocate H (nrowH = maxit+1, ncolH = maxit), row‑major
    int nrowH = maxit + 1, ncolH = maxit;
    double *H   = (double*)std::calloc((size_t)nrowH * ncolH, sizeof(double));

    // 4) Givens data
    double *cs     = (double*)std::malloc(maxit * sizeof(double));
    double *sn     = (double*)std::malloc(maxit * sizeof(double));
    double *s      = (double*)std::calloc(maxit+1, sizeof(double));
    double *resid  = (double*)std::calloc(maxit+1, sizeof(double));

    // temp storage
    double *vnew = (double*)std::malloc(n * sizeof(double));

    // initialize
    s[0]       = beta;
    resid[0]   = beta;
    double tb  = tol * beta;
    // V0 = Mb / beta
    #pragma omp parallel for num_threads(np)
    for(int i = 0; i < n; i++) 
        V[0][i] = Mb[i] / beta;

    int m_actual = 0; // number of necessary iterations

    // GMRES main loop
    for(int m = 1; m <= maxit; m++){
        int j = m - 1;

        // 4.1) vnew = A*V[j], then Jacobi
        matcsrvecprod(n, iat, ja, coef, V[j], vnew, np);
        #pragma omp parallel for num_threads(np)
        for(int i = 0; i < n; i++)
            vnew[i] /= diag[i];

        // 4.2) Arnoldi: orthogonalize vnew against V[0..j]
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
        for(int i = 0; i < n; i++)
            vnorm += vnew[i]*vnew[i];
        vnorm = std::sqrt(vnorm);
        H_at(H,ncolH, j+1, j) = vnorm;       // H(j+1, j)

        if(vnorm < 1e-15){
            m_actual = j+1; //happy breakdown
            break;
        }
        #pragma omp parallel for num_threads(np)
        for(int i = 0; i < n; i++)
            V[m][i] = vnew[i] / vnorm;

        // 4.3) Givens update for column j  
        //    apply previous rotations to H(0:j+1, j)
        for(int k = 0; k < j; k++){
            double h_kj   = H_at(H,ncolH, k,   j),
                   h_k1j  = H_at(H,ncolH, k+1, j);
            double t      = cs[k]*h_kj + sn[k]*h_k1j;
            double t1     = -sn[k]*h_kj + cs[k]*h_k1j;
            H_at(H,ncolH, k,   j) = t;
            H_at(H,ncolH, k+1, j) = t1;
        }
        // build new rotation to zero H(j+1,j)
        double hjj  = H_at(H,ncolH, j,   j),
               hj1j = H_at(H,ncolH, j+1, j),
               rho  = std::sqrt(hjj*hjj + hj1j*hj1j); 
            //    rho = std::hypot(hjj, hj1j); 
        cs[j] =  hjj/rho;
        sn[j] =  hj1j/rho;
        H_at(H,ncolH, j,   j) = rho;
        H_at(H,ncolH, j+1, j) = 0.0;

        // rotate s
        {
            double s_j   = s[j],
                   s_j1  = s[j+1];
            s[j]   =  cs[j]*s_j + sn[j]*s_j1;
            s[j+1] = -sn[j]*s_j + cs[j]*s_j1;
        }
        resid[m] = fabs(s[j+1]);

        // check
        if(resid[m] <= tb){
            m_actual = j+1;
            break;
        }
    }
    if(m_actual == 0) m_actual = maxit;

    printf("Number of iterations: %u\n", m_actual);
    // 5) back‑solve R y = s[0..m_actual-1]  (R is [saved in] H(0:m_actual,0:m_actual))
    double *y = (double*)std::malloc(m_actual * sizeof(double));
    for(int i = m_actual-1; i >= 0; i--){
        double sum = 0.0;
        for(int k = i+1; k < m_actual; k++)
            sum += H_at(H,ncolH, i, k) * y[k];
        y[i] = (s[i] - sum) / H_at(H,ncolH, i, i);
    }

    // 6) form x_out += V[:,0:m_actual-1] * y
    std::memset(x_out, 0, n * sizeof(double));
    #pragma omp parallel for num_threads(np)
    for(int i = 0; i < n; i++){
        double xi = 0.0;
        for(int k = 0; k < m_actual; k++)
            xi += V[k][i] * y[k];
        x_out[i] = xi;
    }

    // --- cleanup ---
    std::free(diag);
    std::free(Mb);
    std::free(Vbuff);
    std::free(V);
    std::free(H);
    std::free(cs);
    std::free(sn);
    std::free(s);
    std::free(resid);
    std::free(vnew);
    std::free(y);
}
