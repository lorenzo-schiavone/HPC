#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

int main(int argc, const char* argv[]){

    int nr, nc; //less local scope
    double **A;
    double *Abuffer; 

    int nv; // dimensione vettore
    double *v;
    double *w;

    double *res;

    if (argc!=4){
        printf("Uso programma: ./pblas matrice vettore vettore2\n");
        exit(1);
        // return 1;
    } else{
        // leggo matrice
        FILE *fid = fopen(argv[1], "r");
        // fscanf ritorna il numero di oggetti che ha letto correttamente
        if (fscanf(fid, "%d %d", &nr, &nc) != 2){
            printf("errore nella lettura della matrice\n");
            exit(2);
        }
        //alloco Abuffer
        A = (double**) malloc(nr*sizeof(double*));
        Abuffer = (double*) malloc(nr*nc*sizeof(double));

        for (int i=0; i<nr;i++){
            A[i] = &Abuffer[(i+1)*nc]; // A[i] è un puntatore a double quindi bisogna assegnargli un riferimento
        }

        for (int i=0; i<nr;i++){
            for (int j=0; j<nc; j++){
                fscanf(fid, "%lf", &A[i][j]); //riferimento dell'elemento da leggere, f sta per float
            }
        }
        fclose(fid);

        // leggo vettore
        fid = fopen(argv[2], "r");
        // fscanf ritorna il numero di oggetti che ha letto correttamente
        if (fscanf(fid, "%d", &nv) != 1){
            printf("errore nella lettura del vettore\n");
            exit(3);
        } else {
        v = (double*) malloc (nv *sizeof(double));

        for (int i=0; i<nv;i++){
            fscanf(fid, "%lf", &v[i]); 
        }

        fclose(fid);

        // leggo vettore2
        fid = fopen(argv[3], "r");
        // fscanf ritorna il numero di oggetti che ha letto correttamente
        if (fscanf(fid, "%d", &nv) != 1){
            printf("errore nella lettura del vettore\n");
            exit(3);
        } else {
        w = (double*) malloc (nv *sizeof(double));

        for (int i=0; i<nv;i++){
            fscanf(fid, "%lf", &w[i]); 
        }

        fclose(fid);
        }
    }
    printf("A: \n");
    // stampa a video della matrice
    for (int i=0; i<nr;i++){
        for (int j=0; j<nc; j++){
            printf("%2.2lf ", A[i][j]); // lf è il formato in doppia precione. f in singola
        }
        printf("\n");
    }

    printf("v: \n");
    // stampa a video della vettore
    for (int i=0; i<nv;i++){
        printf("%2.2lf ", v[i]); // lf è il formato in doppia precione. f in singola
    }
    printf("\n");

    printf("w: \n");
    // stampa a video della vettore
    for (int i=0; i<nv;i++){
        printf("%2.2lf ", w[i]); // lf è il formato in doppia precione. f in singola
    }
    printf("\n");

    // somma vettori v+w
    omp_set_num_threads(2);
    res = (double*) malloc (nv *sizeof(double));
    #pragma omp parallel for
    {
        for (int i = 0; i < nv; i++) {
        res[i] = v[i] + w[i];
        }
    }
    printf("v+w: \n");
    // stampa a video della somma
    for (int i=0; i<nv;i++){
        printf("%2.2lf ", res[i]); // lf è il formato in doppia precione. f in singola
    }
    printf("\n");

    // scalar product
    double result=0;
    #pragma omp parallel for reduction(+:result) 
    {
        for (int i=0; i<nv; i++){
            result += v[i]*w[i];
        }
    }
    printf("prodotto scalare:\n%2.2lf\n", result);
}
}