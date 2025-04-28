#include <cstdio>
#include <stdlib.h>

#include "my_blas.h"

int main(int argc, const char* argv[]){

	// read file nodes
	int nnodes, dim;
	FILE *fid = fopen(argv[1], "r");
    if (!fid) {
        perror("Error opening node file");
        exit(1);
    }
	if  (fscanf(fid, "%d %d" ,&nnodes, &dim)!= 2){
        	printf("Errore lettura\n");
        	exit(2);
	} 
	int u1,u2;
	if  (fscanf(fid, "%d %d" ,&u1, &u2)!= 2){
                printf("Errore lettura\n");
                exit(2); 
        }
	double** nodes = (double **) malloc(nnodes *sizeof(double*)); // one per row
	double* nodesbuf;
	nodesbuf = (double*) malloc ( nnodes*dim * sizeof(double));
	for (int i=0;i<nnodes;i++){
        	nodes[i] = &nodesbuf[i*dim];
	}
	double* node;
	for (int i=0; i<nnodes; i++){
		node = nodes[i];
        fscanf(fid, "%d", &u1);
		for (int j=0; j<dim; j++){
			fscanf(fid, "%lf", &node[j]);
		}
	}
    fclose(fid);

    // read file elements
	int nelem, nodeperelem;
	fid = fopen(argv[2], "r");
    if (!fid) {
        perror("Error opening node file");
        exit(1);
    }
	if  (fscanf(fid, "%d %d" ,&nelem, &nodeperelem)!= 2){
        printf("Errore lettura\n");
        exit(2);
	} 
	int w1;
	if  (fscanf(fid, "%d" ,&w1)!= 1){
            printf("Errore lettura\n");
            exit(2); 
        }
	int** elements = (int **) malloc(nelem *sizeof(int*)); // one per row
	int* elementsbuf;
	elementsbuf = (int*) malloc ( nelem*nodeperelem * sizeof(int));
	for (int i=0;i<nelem;i++){
        elements[i] = &elementsbuf[i*nodeperelem];
	}
	int* element;
	for (int i=0; i<nelem; i++){
		element = elements[i];
        fscanf(fid, "%d", &w1);
		for (int j=0; j<nodeperelem; j++){
			fscanf(fid, "%d", &element[j]);
		}
        fscanf(fid, "%d", &w1);
	}
    fclose(fid);


    /// first scan: compute non zero per row: number of connectivity
    int* connectivities = (int *) malloc(nnodes *sizeof(int)); // one per node
    for (int i = 0; i < nnodes; i++) connectivities[i] = 0;
    bool** link = (bool **) malloc(nnodes * sizeof(bool));
    bool* linkbuf;
	linkbuf = (bool*) malloc ( nnodes*nnodes* sizeof(bool));
	for (int i=0;i<nnodes;i++){
        link[i] = &linkbuf[i*nnodes];
        for (int j=0; j<nnodes; j++ ){
            link[i][j] = false; 
        }
	}

    for (int i=0; i<nelem; i++){
        element = elements[i];
        for (int j=0; j<nodeperelem; j++){
            int nj = (int)element[j] -1;
            if (nj < 0 || nj >= nnodes) {
                printf("Invalid node index %d in element %d\n", element[j], i+1);
                exit(3);
            }
            for (int k=0; k<nodeperelem; k++){
                
                int nk = (int)element[k] -1;
                if (nk >= 0 && nk < nnodes) { 
                    if (!link[nj][nk]){
                        link[nj][nk] = true;
                        connectivities[nj]++;
                    }
                } else {
                    printf("Invalid node index %d in element %d\n", nk + 1, i);
                    exit(3);
                }
            }
        }
    }
    // printf("Connectivities:\n");
    // for (int i=0;i<nnodes;i++){
    //     printf("nodes %d: %d\n", i, connectivities[i]);
    // }
    // IAT
    int* iat;
    iat = (int*) malloc (nnodes*sizeof(int));
    int cumsum=0;
    printf("iat: \n");
    for (int i=0;i<nnodes; i++){
        iat[i]=cumsum;
        cumsum+=connectivities[i];
        printf("%d ", iat[i]);
    }
    printf("\n");
    // MAP for easy access position: one per row, with the length given by the connectivity
    int** map = (int**) malloc(nnodes*sizeof(int*));
    int* mapbuffer = (int*) malloc(cumsum * sizeof(int));
    // i have to compute for each node the global index ordering of the nodes linked to it:
    // eg node 33 is linked to 17 61 27 45 31 33, i should first order them 17 27 31 33 45 61 and then a map 17 to 0, 61 to 2 27 to 3...
    for (int i=0; i< nnodes;i++){
        map[i] = (int*) malloc(connectivities[i]*sizeof(int));
    }


    // ASSEMBLY MATRIX!
    csr A;
    A.nr = nnodes;
    A.iat = iat;


}
