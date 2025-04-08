#include <stdlib.h>
#include <cstdio>
#include <omp.h>


int main(int argc, const char* argv[]){
    int nn,np;

    FILE *fid = fopen(argv[1], "r");
    if  (fscanf(fid, "%d %d" ,&nn,&np)!= 2){
        printf("Errore lettura\n");
        exit(2);
    } 

    printf("Vettori di dimensione %d\n", nn);
    printf("Numero di thread %d\n", np);

    double *u = (double*) malloc (nn*sizeof(double));
    double *v = (double*) malloc (nn*sizeof(double));

    for (int i=0;i<nn; i++){
        // u[i]= double(rand())/ RAND_MAX;
        // v[i]= double(rand())/ RAND_MAX;
        u[i] = 1.;
        v[i] = 1.;
    }

    double *y = (double*) malloc (nn*sizeof(double));
    int bsize = nn / np;

    #pragma omp parallel num_threads(np)
    {
        int myid = omp_get_thread_num();
        int istart = myid*bsize;
        int iend = (myid +1)*bsize; // bilanciare il numero di iterazioni in modo diverso custom
        // prendere anche quelle in fondo
        if (myid == np-1){
            iend = nn;
        }
        // #pragma omp for schedule(static)
        for (int i=istart; i<iend; i++){
            y[i]= v[i] + u[i] + myid;
        } 
    }

    printf("Risultato:\n");
    for (int i = 0; i < nn; i++){
        printf("%2.2f + %2.2f = %2.2f \n",u[i], v[i], y[i]);

    }
    // printf("...\n");

    // for (int i = nn-10; i < nn; i++){
    //     printf("%2.2f + %2.2f = %2.2f \n",u[i], v[i], y[i]);

    // }

    
}