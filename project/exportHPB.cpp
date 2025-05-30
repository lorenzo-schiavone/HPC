#include <iostream> 
#include <stdlib.h>
#include <fstream>
#include <cstdlib>
using namespace std;
#include <omp.h>
#include <chrono>

#include "topol.h"
#include "gmres.h"

double det3(double* r0,double* r1,double* r2 ){
   return r0[0]*(r1[1]*r2[2]-r2[1]*r1[2]) 
        - r1[0]*(r0[1]*r2[2]-r2[1]*r0[2])
        + r2[0]*(r0[1]*r1[2]-r1[1]*r0[2]);
}
double sign(double x){
   return (x>0) - (x<0);
}

// MAIN PROGRAM
int main(int argc, const char* argv[]){

   // Check arguments
   if (argc < 4) {
      printf("Too few arguments.\n Usage: [tetra file] [coord file] [np]\n");
      exit(1);
   }
   int np = atoi(argv[3]);
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
   
   // // save tet again (?)
   // FILE *of = fopen("out_tet","w");
   // for (int i = 0; i < ntet; i++){
   //    for (int j = 0; j < 4; j++) fprintf(of," %10d",tetra[i][j]+1);
   //    fprintf(of,"\n");
   // }
   // fclose(of);

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
   topol(nn,ntet,50,tetra,nterm,iat,ja); // is 30 enough? connectivity per node is ok, it would be very irregular
   printf("Topology created!\n");
   // here we have iat and ja ready to be used
   // -------------------------------------------------------------------------------------------------------------------------------

   // READ COORD
   ifstream coord_file(argv[2]);
   // Read header (i have added it )
   int nnodes;
   coord_file >> nnodes;
   double *coord_buff = (double*) malloc(3*nnodes*sizeof(double)); // we are in 3d
   double **coord = (double**) malloc(nnodes*sizeof(double*));
   k = 0;
   for (int i = 0; i < nnodes; i++){
      coord[i] = &(coord_buff[k]);
      k += 3;
      coord_file >> junk; // this is the index of the node we are throwing away
      for (int j = 0; j < 3; j++){
	      coord_file >> coord[i][j];
      }
   }
   // Close the input file
   coord_file.close();
   
   // -------------------------------------------------------------------------------------------------------------------------------
   std::chrono::time_point<std::chrono::high_resolution_clock> startTime = std::chrono::high_resolution_clock::now();
   // diffusion and flow velocity
   double* D = (double*) malloc(3*sizeof(double));
   D[0]=.1;D[1]=.1;D[2]=1.;
   double* v = (double*) malloc(3*sizeof(double));
   v[0]=1.;v[1]=1.;v[2]=1.;

   // we have to build H, B, P csr matrix.
   // so just coefH, coefB, coefP
   double* coefH = (double*) malloc (nterm * sizeof(double));
   double* coefB = (double*) malloc (nterm * sizeof(double));
   double* coefP = (double*) malloc (nterm * sizeof(double));
   for (int i=0; i< nterm; i++){
      coefH[i]=0;
      coefB[i]=0;
      coefP[i]=0;
   }

   # pragma omp parallel for num_threads(np)
   for (int jj=0;jj<1;jj++){
      int* tet = tetra[jj];
      double* n0=coord[tet[0]];
      double* n1=coord[tet[1]];
      double* n2=coord[tet[2]];
      double* n3=coord[tet[3]];       
      double vol = (det3(n1,n2,n3)-
                   det3(n0,n2,n3)+
                   det3(n0,n1,n3)-
                   det3(n0,n1,n2))/6;
      printf("volume: %f\n", vol);
      double* a = (double*) malloc(4*sizeof(double));
      double* b = (double*) malloc(4*sizeof(double));
      double* c = (double*) malloc(4*sizeof(double));
      double* d = (double*) malloc(4*sizeof(double));

      for (int ii=0;ii<4;ii++){
         uint j = tet[(ii+1)%4];
         uint k = tet[(ii+2)%4];
         uint m = tet[(ii+3)%4];
         printf("%u %u %u\n",j+1,k+1,m+1);
         double* nj=coord[j];
         double* nk=coord[k];
         double* nm=coord[m];

         // print matrix [nj,nk,nm]
         printf("M: ");
         for (int i=0;i<3;i++){
            printf("%f ", nj[i]);
         } printf("\n");
         for (int i=0;i<3;i++){
            printf("%f ", nk[i]);
         } printf("\n");
         for (int i=0;i<3;i++){
            printf("%f ", nm[i]);
         } printf("\n");



         a[ii] = det3(nj,nk,nm);
         b[ii] = - (nk[1]*nm[2]-nm[1]*nk[2] - (nj[1]*nm[2]-nj[2]*nm[1])+nj[1]*nk[2]-nj[2]*nk[1]);
         c[ii] = (nk[0]*nm[2]-nm[0]*nk[2] - (nj[0]*nm[2]-nj[2]*nm[0])+nj[0]*nk[2]-nj[2]*nk[0]);
         d[ii] = (nk[1]*nm[0]-nm[1]*nk[0] - (nj[1]*nm[0]-nj[0]*nm[1])+nj[1]*nk[0]-nj[0]*nk[1]); // I swap col 1 with col 2 so the det change sign
      }

      // print a, b, c, d
      printf("a: \n");
      for (int i=0;i<4;i++){
         printf("%f ", a[i]);
      }
      printf("\n");
      printf("b: \n");
      for (int i=0;i<4;i++){
         printf("%f ", b[i]);
      }
      printf("\n");
      printf("c: \n");
      for (int i=0;i<4;i++){
         printf("%f ", c[i]);
      }
      printf("\n");
      printf("d: \n");
      for (int i=0;i<4;i++){
         printf("%f ", d[i]);
      }
      printf("\n");

      
      double *Hloc_buf = (double*) malloc(4*4*sizeof(double));
      double **Hloc = (double**) malloc(4*sizeof(double*));
      double *Ploc_buf = (double*) malloc(4*4*sizeof(double));
      double **Ploc = (double**) malloc(4*sizeof(double*));
      double *Bloc_buf = (double*) malloc(4*4*sizeof(double));
      double **Bloc = (double**) malloc(4*sizeof(double*));
      int kk=0;
      for (int i=0;i<4;i++){
         Hloc[i] = &(Hloc_buf[kk]);
         Ploc[i] = &(Ploc_buf[kk]);
         Bloc[i] = &(Bloc_buf[kk]);
         kk+=4;
      }

      for (int i=0;i<4; i++){
         double* Hi = Hloc[i];
         double* Pi = Ploc[i];
         double* Bi = Bloc[i];
         for (int j=0;j<4;j++){
            Hi[j] = (D[0]*b[i]*b[j] + D[1]*c[i]*c[j] + D[2]*d[i]*d[j])/ (36 * abs(vol));
            Pi[j] = abs(vol)/20;
            Bi[j] = sign(vol)/24 * (v[0]*b[j]+v[1]*c[j]+v[2]*d[j]);
         }
         Pi[i]*=2;
      }
      
      // put them inside coefH, coefP, coefB
      // omp atomic or whatever here like assembly only row for nodes deputed to the process
      // maybe later instead of first building the local matrix and then put it inside coef, put it directly there

      for (int i=0;i<4; i++){
         // tet[i] gives the row index
         double* Hi = Hloc[i];
         double* Pi = Ploc[i];
         double* Bi = Bloc[i];
         for (int j=0;j<4;j++){
            // in ja from iat[tet[i]] look for tet[j]
            int curr_node = tet[i];
            for (int ii=iat[curr_node]; ii<iat[curr_node+1]; ii++){
               int jjj = tet[j];
               if (ja[ii] == jjj){
                  #pragma omp atomic
                  coefH[ii] += Hi[j];
                  #pragma omp atomic
                  coefP[ii] += Pi[j];
                  #pragma omp atomic
                  coefB[ii] += Bi[j];
               }
            }
         }
      }  
   }
   std::chrono::time_point<std::chrono::high_resolution_clock> endTime = std::chrono::high_resolution_clock::now();
   std::chrono::duration<double> timeTaken = std::chrono::duration<double>(endTime - startTime);

   printf("assembly time taken %f\n", timeTaken.count());

   // Export P,B, H
   FILE *f;

    // Export iat
    f = fopen("iat.txt", "w");
    for (int i = 0; i <= nn; ++i)
        fprintf(f, "%d\n", iat[i]);
    fclose(f);

    // Export ja
    f = fopen("ja.txt", "w");
    for (int i = 0; i < nterm; ++i)
        fprintf(f, "%d\n", ja[i]);
    fclose(f);

    // Export coef
    f = fopen("coefH.txt", "w");
    for (int i = 0; i < nterm; ++i)
        fprintf(f, "%.16g\n", coefH[i]);
    fclose(f);

    // Export coef
    f = fopen("coefB.txt", "w");
    for (int i = 0; i < nterm; ++i)
        fprintf(f, "%.16g\n", coefB[i]);
    fclose(f);

    // Export coef
    f = fopen("coefP.txt", "w");
    for (int i = 0; i < nterm; ++i)
        fprintf(f, "%.16g\n", coefP[i]);
    fclose(f);

    // Export nnz
    f = fopen("nnz.txt", "w");
    fprintf(f, "%d\n", nterm);
    fclose(f);
   
   
   // Free memory
   free(coord);
   free(coord_buff);
   free(tetra);
   free(tetra_buf);
   free(iat);
   free(ja);
   free(coefH);
   free(coefP);
   free(coefB);

}
