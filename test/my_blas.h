#ifndef MY_BLAS_H
#define MY_BLAS_H

// CSR format structure (only used if you need sparse matrix–vector products)
struct csr {
    int* iat;     // Row pointers.
    int* ja;      // Column indices.
    double* coef; // Nonzero coefficients.
    int nr;       // Number of rows.
};

void sumvec(double* v, double* u, double* res, int nn);
void daxpy(double* v, double* u, double a, int nn);
double scalarprod(double* v, double* u, int nn);
void matvecprod(double** A, double* v, double* res, int nc, int nr);
void matmatprod(double** A, double** B, double** C, int m, int n, int p);
void transpose(double** A, double** At, int nr, int nc);
void sparsematvecprod(csr& A, double* v, double* res);

// QR factorization of A = Q * R.
// A: nr x nc, Q: nr x nc, R: nc x nc.
void qr(double** A, double** Q, double** R, int nr, int nc);

#endif // MY_BLAS_H
