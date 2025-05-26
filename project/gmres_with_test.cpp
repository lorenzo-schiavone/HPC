#include <cmath>
#include <cstdlib>
#include <fstream>

#include "gmres.h"
#include "random_csr.h"

void testGMRES() {
    printf("\nTesting GMRES decomposition...\n");
    int nrows = 100;
    double* coef = nullptr;
    int* ja = nullptr;
    int* iat = nullptr;
    int nnz = 3 * nrows - 2; 

    generate_stiffness_csr(nrows, coef, ja, iat);
    printf("stiffness created\n");
    double* v = (double*) malloc (nrows*sizeof(double));
    double* rhs = (double*) malloc (nrows*sizeof(double));
    double* x = (double*) malloc (nrows*sizeof(double));
    for (int i=0; i<nrows; i++){
        v[i]=1.;
        rhs[i]=0.;
    }

    // printf("iat: ");
    // for (int i=0; i<nrows+1; i++){
    //     printf("%d ", iat[i]);
    // }
    // printf("\n");

    // printf("ja: ");
    // for (int i=0; i<nnz; i++){
    //     printf("%d ", ja[i]);
    // }
    // printf("\n");

    // printf("coef: ");
    // for (int i=0; i<nnz; i++){
    //     printf("%f ", coef[i]);
    // }
    // printf("\n");

    // printf("%f\n", coef[29]);
    matcsrvecprod(nrows, iat, ja, coef, v, rhs, 4);
    printf("product done created\n");
    // printf("v: ");
    // for (int i=0; i<nrows; i++){
    //     printf("%f ", v[i]);
    // }
    // printf("\n");
    printf("rhs: ");
    for (int i=0; i<nrows; i++){
        printf("%f ", rhs[i]);
    }
    printf("\n");
    gmres(nrows,iat,ja,coef, rhs, 1e-14,100, 4, x);
    // print x
    printf("x: ");
    for (int i=0; i<nrows; i++){
        printf("%f ", x[i]);
    }
    printf("\n");
}

int main(){
    testGMRES();
    return 0;
}