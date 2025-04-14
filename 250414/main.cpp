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
        u[i]= double(rand())/ RAND_MAX;
        v[i]= double(rand())/ RAND_MAX;
        // u[i] = 1.;
        // v[i] = 2.;
    }

    printf("u         v\n");
    for (int i=0; i<10; i++){
        printf("%10.6f %10.6f\n", u[i],v[i]);
    }

    double *y = (double*) malloc (nn*sizeof(double));
    int bsize = nn / np;

    double alpha = 0.; // alpha inizzializzato con l'elemento neutro dell'operatore
    #pragma omp parallel num_threads(np) reduction(+:alpha) // riduzione scelta da lui. possiamo essere un po' più precisi
    {
        int myid = omp_get_thread_num();
        int istart = myid*bsize;
        int iend = (myid +1)*bsize; 
        if (myid == np-1){
            iend = nn;
        }
        // #pragma omp for schedule(static)
        for (int i=istart; i<iend; i++){
            alpha+=v[i]*u[i];
        } 
    }

    printf("Prodotto scalare1:\n");
    printf("%10.15e\n", alpha);

    alpha = 0.;
    double* ridv = (double*) malloc(np*sizeof(double)); // uno per processore usato
    printf("Prodotto scalare2:\n");

    #pragma omp parallel num_threads(np) reduction(+:alpha) // riduzione scelta da lui. possiamo essere un po' più precisi
    {
        int myid = omp_get_thread_num();
        int istart = myid*bsize;
        int iend = (myid +1)*bsize; 
        if (myid == np-1){
            iend = nn;
        }
        // #pragma omp for schedule(static)
        double alpha_loc;
        for (int i=istart; i<iend; i++){
            alpha_loc+=v[i]*u[i];
        } 
        ridv[myid] = alpha_loc;
        #pragma omp barrier
        // #pragma omp single nowait // non ha nessuna importanza avere una barriera dopo; ha senso metterlo dentro la regione parallela perché se non
        // {
        //     for (int j=0; j<np;j++){
        //         alpha+=ridv[j];
        //     }
        // }
        // fintanto che alpha è privata ogni processo può fare la sua riduzione
        alpha_loc=0; // utile per il gradiente coniugato. tutti i trhead devono avere il prodotto scalare al loro interno in una variabile private
        for (int j=0; j<np;j++){ // pecca: riduzione lineare - a imbuto sarebbe più veloce ma con pochi processori cambia poco
            alpha_loc+=ridv[j];
        } // riduzione comandata
        printf("id %u: %10.15e\n", myid, alpha_loc);
    }


    



    return 0;
}