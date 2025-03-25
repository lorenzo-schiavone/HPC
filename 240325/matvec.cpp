#include <stdlib.h>
#include <stdio.h>

int main(int argc, const char* argv[]){

    int nr, nc; //less local scope
    double **A;
    double *Abuffer; 

    int nv; // dimensione vettore
    double *v;
    if (argc!=3){
        printf("Uso programma: ./matvet matrice vettore\n");
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
        A = (double**) malloc(nr*sizeof(double *));
        Abuffer = (double*) malloc(nr*nc*sizeof(double));

        for (int i=0; i<nr;i++){
            A[i] = &Abuffer[(i+1)*nc]; // A[i] è un puntatore a double quindi bisogna asssegnargli un riferimento
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
        }
    }
    // stampa a video della matrice
    for (int i=0; i<nr;i++){
        for (int j=0; j<nc; j++){
            printf("%2.2lf ", A[i][j]); // lf è il formato in doppia precione. f in singola
        }
        printf("\n");
    }

    // stampa a video della vettore
    for (int i=0; i<nv;i++){
        printf("%2.2lf ", v[i]); // lf è il formato in doppia precione. f in singola
    }
    printf("\n");
}