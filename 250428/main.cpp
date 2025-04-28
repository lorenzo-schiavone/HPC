#include <cstdio>
#include <stdlib.h>
#include <chrono>

#include <omp.h>
#include "readMat.h"

int main(int argc, const char* argv[]){

    int np;

    if (argc != 3){
        printf("uso programma: nome_matrice, np\n");
        exit(1);
    }
    np = atoi(argv[2]);
    int nrows, ncols, nterm;
    int* ia;
    int* ja;
    double* coef;

    bool BINREAD = false;
    
    readCSRmat(&nrows, &ncols,&nterm, &ia, &ja, &coef, (char *) argv[1], BINREAD);
    // printf("nrows: %u\n", nrows);
    // printf("ncols: %u\n", ncols);
    // printf("nterm: %u\n", nterm);

    
    // CSR A;
    // A.nr = nrows; A.nt = nterm;
    // A.iat = ia; A.ja = ja; A.coef = coef;

    double* v = (double*) malloc(ncols*sizeof(double));
    for (int i = 0; i < ncols; i++) {
        // fscanf(fid, "%lf", &(v[i]));
        v[i]=i+1;
    }
        // fclose(fid);

    double* res = (double*) malloc(nrows*sizeof(double));
    auto startTime = std::chrono::high_resolution_clock::now();
    #pragma omp parallel for num_threads(np)
    for (int i = 0; i < nrows; i++){
        res[i] = 0.0;
        // printf("%u:\n", i);
        for (int j = ia[i]; j < ia[i+1]; j++){
            // printf("%f ", coef[j]);
            // printf("ja: %u ", ja[j]);
            res[i] += (coef[j] * v[ja[j]]);
            // printf("%f ", res[i]);
        }
        // printf("\n");
    }
    auto endTime = std::chrono::high_resolution_clock::now();
    auto timeTaken = std::chrono::duration<double>(endTime - startTime);

    printf("time taken %f\n", timeTaken.count());
    // printf("Av:\n");
    // for (int i = 0; i < nrows; i++){
    //     printf("%e\n",res[i]);

    // }
    return 0;
}
