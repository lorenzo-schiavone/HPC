#include <iostream> 
#include <stdlib.h>
#include <fstream>
#include <cstdlib>
using namespace std;

#include "topol.h"

// MAIN PROGRAM
int main(int argc, const char* argv[]){

   // Check arguments
   if (argc < 2) {
      printf("Too few arguments.\n Usage: [tetra file]\n");
      exit(1);
   }

   // Read tetrahedrons
   ifstream tetra_file(argv[1]);
   // Read header
   int ntet;
   tetra_file >> ntet;
   int *tetra_buf = (int*) malloc(4*ntet*sizeof(int));
   int **tetra = (int**) malloc(ntet*sizeof(int*));
   int k = 0;
   int junk;
   for (int i = 0; i < ntet; i++){
      tetra[i] = &(tetra_buf[k]);
      k += 4;
      tetra_file >> junk;
      for (int j = 0; j < 4; j++){
	 tetra_file >> tetra[i][j];
      }
      tetra_file >> junk;
   }
   // Close the input file
   tetra_file.close();

   // Set C-style of tet
   for (int i = 0; i < ntet; i++)
      for (int j = 0; j < 4; j++) tetra[i][j]--;

   FILE *of = fopen("out_tet","w");
   for (int i = 0; i < ntet; i++){
      for (int j = 0; j < 4; j++) fprintf(of," %10d",tetra[i][j]+1);
      fprintf(of,"\n");
   }
   fclose(of);

   // Find the number of equations
   int nn = 0;
   for (int i = 0; i < ntet; i++)
      for (int j = 0; j < 4; j++) 
         if (tetra[i][j] > nn) nn = tetra[i][j];
   nn++;

   printf("Number of equations: %d\n",nn);

   // Create topology
   int nterm;
   int *iat = nullptr;
   int *ja = nullptr;
   topol(nn,ntet,30,tetra,nterm,iat,ja);

   // Dump topology
   of = fopen("out_pat","w");
   for (int i = 0; i < nn; i++){
      // Increase by 1 row and column indices for MATLAB visualization
      for (int j = iat[i]; j < iat[i+1]; j++) fprintf(of,"%10d %10d 1.0\n",i+1,ja[j]+1);
   }
   fclose(of);

   // Free memory
   free(tetra);
   free(tetra_buf);
   free(iat);
   free(ja);

}
