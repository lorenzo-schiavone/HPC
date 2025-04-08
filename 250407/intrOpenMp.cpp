#include <cstdio>    //impostazioni di input e output sono incluse in queste librerie
#include <stdlib.h>
// #include "somma.h"

int main(int argc, const char* argv[]) {    //gli argomenti servono per passare gli argomenti al codice nel momento del lancio
    int nr, nc, nv;
    double** A;
    double* Abuf;
    double* v;
    if (argc != 3){
        printf("Uso programma: ./MatVet Matrice Vettore\n");
        exit(1);

    } else {
        //Leggo la matrice
        FILE *fid = fopen(argv[1], "r");
        if  (fscanf(fid, "%d %d" ,&nr,&nc)!= 2){
            printf("Errore lettura\n");
            exit(2);
        } 
        A = (double**) malloc(nr*sizeof(double*));
        Abuf = (double*) malloc (nr*nc*sizeof(double));
        for (int i = 0; i < nr; i++) {
            A[i] = &(Abuf[i*nc]);
        }
        for (int i= 0; i < nr; i++){
            for (int j = 0; j < nc; j++)
            fscanf(fid, "%lf",&(A[i][j]));
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
    

    printf("Matrice:\n");
    for (int i = 0; i < nr; i++){
        for (int j = 0; j < nc; j++) printf(" %2.2f", A[i][j]);
        printf("\n");
            
    }

    printf("Vettore:\n");
    for (int i = 0; i < nv; i++){
        printf(" %2.2f\n",v[i]);

    }

//Effettuiamo il prodotto MxV

double *x = (double*) malloc(nr*sizeof(double));
for (int i = 0; i < nr; i++){
    x[i] = 0.0;
    for (int j = 0; j < nc; j++){
        x[i] += A[i][j]*v[j];

    }
}

printf("MxV:\n");
for (int i = 0; i < nr; i++){
     printf(" %10.2f\n",x[i]);

 }

double *y = (double*) malloc(nr*sizeof(double));
int np = 4;
#pragma omp parallel num_threads(np)
{
    #pragma omp for
    for (int i=0; i<nr; i++){
        y[i]= v[i] + x[i];
    } 
}

printf("Vettore somma:\n");
for (int i = 0; i < nr; i++){
    printf(" %2.2f\n",y[i]);

}
}

// g++ MatVet.cpp -o MatVet.x
//