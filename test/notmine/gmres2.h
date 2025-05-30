void matcsrvecprod(int nn,int* iat,int* ja,double* coef,double* v, double* x, int np);

void gmres(int n, int* iat, int* ja, double* coef,
           double* rhs, double tol, int maxit,
           int np, double* x);