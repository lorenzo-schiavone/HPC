void gmres(int nrows, int* iat, int* ja, double* coef, double* rhs, double tol, int maxit, int np, double* x);
void matcsrvecprod(int nn,int* iat,int* ja,double* coef,double* v, double* x, int np);