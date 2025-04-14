#include <cstdio>
#include <stdlib.h>
#include <omp.h>

struct CSR{
    int nr;
    int nt;
    int* iat;
    int* ja;
    double* coef;
};

int main(int argc, const char* argv[]) {  
    CSR A;
    int nr, nt, nv;
    double* v;
    int np=4;
    if (argc != 3){
        printf("Uso programma: ./MatCSRVet Matrice Vettore\n");
        exit(1);

    } else {
        //Leggo la matrice
        FILE *fid = fopen(argv[1], "r");
        if (fid == NULL) {
            printf("Errore: impossibile aprire il file %s\n", argv[1]);
            exit(1);
        }

        if  (fscanf(fid, "%d %d" ,&nr,&nt)!= 2){
            printf("Errore lettura\n");
            exit(2);
        }

        printf("nr nt %d %d", nr, nt);
        //Leggo il vettore
        A.nr = nr;
        A.nt = nt;
        A.iat = (int*) malloc((nr+1)*sizeof(int));
        A.ja = (int*) malloc(nt*sizeof(int));
        A.coef = (double*) malloc(nt*sizeof(double));
        printf("here\n");
                 
        for (int i = 0; i < nr+1; i++) {
            fscanf(fid, "%d",&(A.iat[i]));
        }

        for (int i = 0; i < nt; i++){
            fscanf(fid, "%d",&(A.ja[i]));
        }

        for (int i = 0; i < nt; i++) {
            fscanf(fid, "%lf",&(A.coef[i]));
        }


        fclose(fid);
        //Leggo il vettore
        fid = fopen(argv[2], "r");
        if (fscanf(fid, "%d" , &nv) != 1){
            printf("Errore lettura\n");
            exit(3);
        }
        v = (double*) malloc(nv*sizeof(double));
        for (int i = 0; i < nv; i++) {
            fscanf(fid, "%lf", &(v[i]));
        }
        fclose(fid);
    }
    

    printf("Matrice in formato coordinate:\n");
    printf("%d nr" , nr);
    for (int i = 0; i < nr; i++){
        printf("--> %d %d\n",A.iat[i],A.iat[i+1]);
            for (int j = A.iat[i]; j < A.iat[i+1]; j++){
                printf("%4d %4d %10.4f\n" ,i,A.ja[j],A.coef[j]);
            }
    }

    printf("Vettore:\n");
    for (int i = 0; i < nv; i++){
        printf(" %10.2f\n",v[i]);

    }

//Effettuiamo il prodotto MxV
    int *iat = A.iat;
    int*ja = A.ja;
    double *coef = A.coef;
    double *x = (double*) malloc(nr*sizeof(double));

    // se tanti processi ogni tanto controllare che vettori siano sincronizzato
    #pragma omp parallel for num_threads(np)
    {
    for (int i = 0; i < nr; i++){
        x[i] = 0.0;
        for (int j = iat[i]; j < iat[i+1]; j++){
            x[i] += coef[j] * v[ja[j]];
        }
    }
    }

    printf("MxV:\n");
    for (int i = 0; i < nr; i++){
        printf(" %10.2f\n",x[i]);

    }

}

// g++ MatCSRVet.cpp -o MatCSRVet.x
// ./a.out MatCSRVer.x Matrix_CSR Vector_CSR
