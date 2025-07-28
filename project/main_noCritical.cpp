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
   
   // Find the number of equations
   int nn = 0;
   for (int i = 0; i < ntet; i++)
      for (int j = 0; j < 4; j++) 
         if (tetra[i][j] > nn) nn = tetra[i][j];
   nn++; // C style

   printf("Number of equations: %d\n",nn);

   // Create topology
   int nterm;
   int *iat = nullptr;
   int *ja = nullptr;
   topol(nn,ntet,50,tetra,nterm,iat,ja); 
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
   D[0]=0.02;D[1]=0.01;D[2]=0.01;
   double* v = (double*) malloc(3*sizeof(double));
   v[0]=1.;v[1]=1.;v[2]=1.5;

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

   int bsize = ntet / np;
   # pragma omp parallel num_threads(np)
   {  
      int myid = omp_get_thread_num();
      int istart = myid*bsize;
      int iend = (myid +1)*bsize; 
      if (myid == np-1){
         iend = ntet;
      }

      // thread private alloction
      double* a = (double*) malloc(4*sizeof(double));
      double* b = (double*) malloc(4*sizeof(double));
      double* c = (double*) malloc(4*sizeof(double));
      double* d = (double*) malloc(4*sizeof(double));
      double *Hloc_buf = (double*) malloc(4*4*sizeof(double));
      double **Hloc = (double**) malloc(4*sizeof(double*));
      double *Ploc_buf = (double*) malloc(4*4*sizeof(double));
      double **Ploc = (double**) malloc(4*sizeof(double*));
      double *Bloc_buf = (double*) malloc(4*4*sizeof(double));
      double **Bloc = (double**) malloc(4*sizeof(double*));
      for (int i=0;i<4;i++){
         Hloc[i] = &(Hloc_buf[i*4]);
         Ploc[i] = &(Ploc_buf[i*4]);
         Bloc[i] = &(Bloc_buf[i*4]);
      }

      for (int jj=0;jj<ntet;jj++){
         int* tet = tetra[jj];
         double* n0=coord[tet[0]];
         double* n1=coord[tet[1]];
         double* n2=coord[tet[2]];
         double* n3=coord[tet[3]];  
         bool test[4];
         for (int i=0;i<4;i++){
            test[i] = (tet[i]<iend)&&(tet[i]>=istart);
         }
         // test[0] = (n0<iend)&&(n0>=istart);  
         // test[1] = (n1<iend)&&(n1>=istart);  
         // test[2] = (n2<iend)&&(n2>=istart);  
         // test[3] = (n3<iend)&&(n3>=istart);
         if (!(test[0]||test[1]||test[2]||test[3])) {
            continue; // no node in the range
         }
         double vol = (det3(n1,n2,n3)-
                     det3(n0,n2,n3)+
                     det3(n0,n1,n3)-
                     det3(n0,n1,n2))/6;
         

         for (int ii=0;ii<4;ii++){
            int idx[3]; int cnt=0;
            for (int jj=0;jj<4;j++){
               if (ii==jj){
                  continue;
               }
               idx[cnt]=jj;
               cnt+=1;
            }
            double segno = (-1)^(ii+1);
            double* nj=coord[tet[(ii+1)%4]];
            double* nk=coord[tet[(ii+2)%4]];
            double* nm=coord[tet[(ii+3)%4]];
            a[ii] = segno * det3(nj,nk,nm);
            b[ii] = segno * (-1) * (nk[1]*nm[2]-nm[1]*nk[2] - (nj[1]*nm[2]-nj[2]*nm[1])+nj[1]*nk[2]-nj[2]*nk[1]);
            c[ii] = segno * (nk[0]*nm[2]-nm[0]*nk[2] - (nj[0]*nm[2]-nj[2]*nm[0])+nj[0]*nk[2]-nj[2]*nk[0]);
            d[ii] = segno * (nk[1]*nm[0]-nm[1]*nk[0] - (nj[1]*nm[0]-nj[0]*nm[1])+nj[1]*nk[0]-nj[0]*nk[1]); // I swap col 1 with col 2 so the det change sign
         }

         for (int i=0;i<4; i++){
            if (!(test[i])){
               continue;
            }
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
            int curr_node = tet[i];
            if (test[i]){
               double* Hi = Hloc[i];
               double* Pi = Ploc[i];
               double* Bi = Bloc[i];
               for (int j=0;j<4;j++){
               // in ja from iat[tet[i]] look for tet[j]
               
                  for (int ii=iat[curr_node]; ii<iat[curr_node+1]; ii++){
                     int jjj = tet[j];
                     if (ja[ii] == jjj){
                        coefH[ii] += Hi[j];
                        coefP[ii] += Pi[j];
                        coefB[ii] += Bi[j];
                     }
                  }
               }
            }
         }  
      }
   }
   std::chrono::time_point<std::chrono::high_resolution_clock> endTime = std::chrono::high_resolution_clock::now();
   std::chrono::duration<double> timeTaken = std::chrono::duration<double>(endTime - startTime);

   printf("assembly time taken %f\n", timeTaken.count());

   //// HERE WE TRY TO MAKE A GMRES
   startTime = std::chrono::high_resolution_clock::now();
   double* coefA = (double*) malloc (nterm * sizeof(double));
   double dt = 0.1;
   for (int i=0;i<nterm; i++){
      coefA[i] = coefB[i]+coefH[i] + coefP[i]/dt;
   }
   // printf("coefA build\n");
   double* x_true = (double*) malloc(nnodes*sizeof(double));
   for(int i=0;i<nnodes;i++){
      x_true[i] = 1.0;
   }

   double* q = (double*) malloc(nnodes*sizeof(double));
   matcsrvecprod(nnodes, iat, ja, coefA, x_true, q, np);
   // printf("rhs build\n");


   double* x = (double *) malloc(nnodes*sizeof(double));
   for (int i=0;i<nnodes;i++){
      x[i] = 0.;
   }
   double tol = 1e-13;
   int maxit = 100;
   gmres(nnodes, iat, ja, coefA, q, tol, maxit, np, x);
   endTime = std::chrono::high_resolution_clock::now();
   timeTaken = std::chrono::duration<double>(endTime - startTime);

   printf("gmres time taken %f\n", timeTaken.count());
   printf("q: \n");
   for (int i=0; i<20; i++){
      printf("%f ", q[i]);
   }
   printf("\n");
   printf("x: \n");
   for (int i=0; i<20; i++){
      printf("%f ", x[i]);
   }
   printf("\n");
   
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
