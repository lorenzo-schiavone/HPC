// include library
#include <cstdio>
#include <stdlib.h>

// include my function
#include "somma.h"  // signature o firma necessaria per compilare

int main(int argc, const char* argv[]){
    printf("Hello World\n"); 
    /*argc variabile intera che contiene il numero di parametri inserito. 
    IL primo parametro è sempre il nome del programma. arg è un vettore di stringe che contiene i parametri */
    printf("il numero di parametri è %d\n", argc);
    for (int i=0; i<argc; i++){
        // printf("%s\n", argv[i]); //%s perchè sono stringe
        if (i != 0){
        int par = atoi(argv[i]); // se non è in grado di convertire da 0.
        printf("in formato intero: %d\n", par);
        }
    }

    int sum = somma(10,20);
    int a = 10;
    int b = 20;

    printf("10 + 20 = %d\n", sum);

    printf("%d + %d = %d\n", a, b, somma_ref_const(a,b));
    somma_ref(a,b);

    printf( "a= %d\n", a );
    printf("%d + %d = %d\n", a, b, somma_ref_const(a,b));
}