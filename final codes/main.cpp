#include <iostream> 
#include <stdlib.h>
#include <fstream>
#include <cstdlib>
using namespace std;
#include <math.h>
#include <omp.h>
#include <chrono>

#include "topol.h"
// #include "gmres.h"
#include "gmres_restart.h"

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
   if (argc < 6) {
      printf("Too few arguments.\n Usage: [tetra file] [coord file] [np] [maxit] [restart]\n");
      exit(1);
   }
   int np = atoi(argv[3]);
   int maxit = atoi(argv[4]);
   int restart = atoi(argv[5]);
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
   nn++;

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
   int nnodes; // potrei usare nn. ntet invece va imposto a mano
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

   // diffusion and flow velocity
   double* D = (double*) malloc(3*sizeof(double));
   D[0]=.4;D[1]=.1;D[2]=.1;
   double* v = (double*) malloc(3*sizeof(double));
   v[0]=1.;v[1]=1.;v[2]=2;

   printf("D:\n");
   printf("%f %f %f\n", D[0],D[1],D[2]);
   printf("v:\n");
   printf("%f %f %f\n", v[0],v[1],v[2]);

   // ASSEMBLY
   double start;
   double end;
   start = omp_get_wtime();  
   // we have to build H, B, P csr matrix. -> so just coefH, coefB, coefP
   double* coefH = (double*) malloc (nterm * sizeof(double));
   double* coefB = (double*) malloc (nterm * sizeof(double));
   double* coefP = (double*) malloc (nterm * sizeof(double));
   for (int i=0; i< nterm; i++){
      coefH[i]=0;
      coefB[i]=0;
      coefP[i]=0;
   }

   # pragma omp parallel num_threads(np)
   {
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
      int kk=0;
      for (int i=0;i<4;i++){
         Hloc[i] = &(Hloc_buf[kk]);
         Ploc[i] = &(Ploc_buf[kk]);
         Bloc[i] = &(Bloc_buf[kk]);
         kk+=4;
      }
      # pragma omp for 
      for (int jj=0;jj<ntet;jj++){
         int* tet = tetra[jj];
         double* n0=coord[tet[0]];
         double* n1=coord[tet[1]];
         double* n2=coord[tet[2]];
         double* n3=coord[tet[3]];       
         double vol = (det3(n1,n2,n3)-
                     det3(n0,n2,n3)+
                     det3(n0,n1,n3)-
                     det3(n0,n1,n2))/6;
         
         for (int ii=0;ii<4;ii++){
            int idx[3]; int cnt=0;
            for (int jj=0;jj<4;jj++){
               if (ii==jj){
                  continue;
               }
               idx[cnt]=jj;
               cnt+=1;
            }
            double segno = pow(-1, ii); 
            double* nj=coord[tet[idx[0]]];
            double* nk=coord[tet[idx[1]]];
            double* nm=coord[tet[idx[2]]];
            a[ii] = segno * det3(nj,nk,nm);
            b[ii] = segno * (-1) * (nk[1]*nm[2]-nm[1]*nk[2] - (nj[1]*nm[2]-nj[2]*nm[1])+nj[1]*nk[2]-nj[2]*nk[1]);
            c[ii] = segno * (nk[0]*nm[2]-nm[0]*nk[2] - (nj[0]*nm[2]-nj[2]*nm[0])+nj[0]*nk[2]-nj[2]*nk[0]);
            d[ii] = segno * (nk[1]*nm[0]-nm[1]*nk[0] - (nj[1]*nm[0]-nj[0]*nm[1])+nj[1]*nk[0]-nj[0]*nk[1]); // I swap col 1 with col 2 so the det change sign
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
   }
   end = omp_get_wtime();
   printf("assembly time taken %f\n", end - start);

   //// GMRES
   start = omp_get_wtime();
   double* coefA = (double*) malloc (nterm * sizeof(double));
   double* coefAmod = (double*) malloc (nterm * sizeof(double)); // to be modified for imposing dirichlet bc

   double dt = 1;
   for (int i=0;i<nterm; i++){
      coefA[i] = coefB[i]+coefH[i] + coefP[i]/dt; // if multiple steps i would need it
      coefAmod[i] = coefA[i]; 
   }


   double* bdval = (double*) malloc(nnodes*sizeof(double));
   for (int i=0;i<nnodes;i++){
      // bdval[i]=0.;
      if ((coord[i][0]< 1e-10)&&(coord[i][1]< 1e-10)&&(coord[i][2]< .3 + 1e-10)){
         bdval[i]=1.;
      }
      else {
         bdval[i]=0.;
      }
   }

   double* Abdval = (double *) malloc(nnodes*sizeof(double));
   for (int i=0;i<nnodes;i++){
      Abdval[i] = 0.;
   }
   matcsrvecprod(nn,iat,ja,coefA, bdval, Abdval, np); // xx is the boundary adjust


   // boundary impose
   for (int iii = 0; iii < nn; iii++){
      if (bdval[iii]){ // bc to impose, set the row equal to zero except for diagonal element
         for (int j=iat[iii]; j<iat[iii+1];j++){
            if (ja[j]==iii){
               coefAmod[j]=1.;
            }
            else{
               coefAmod[j]=0.;
            }
         }
      }
      else{ // set columns where q[j]=1 equal to 0;
         for (int j=iat[iii]; j<iat[iii+1];j++){
            if (bdval[ja[j]]){
               coefAmod[j]=0.;
            }
         }
      }
   }

   // building rhs
   double* rhs = (double *) malloc(nnodes*sizeof(double));
   for (int i=0;i<nnodes;i++){
      rhs[i] = 0.;
   }

   matcsrvecprod(nn,iat,ja, coefP, bdval, rhs, np);
   for (int i=0;i<nnodes;i++){
      if (bdval[i]){
         rhs[i] = bdval[i];
      }
      else{
         rhs[i] = rhs[i]/dt - Abdval[i]; // (P/dt)u_prev - Abdval
      }
   } 

   // vector for the solution
   double* sol = (double *) malloc(nnodes*sizeof(double));
   for (int i=0;i<nnodes;i++){
      sol[i] = 0.;
   }

   double tol = 1e-9;

   
   gmres(nnodes, iat, ja, coefAmod, rhs, tol, maxit, restart, np, sol);
   end = omp_get_wtime();
   printf("gmres time taken %f\n", end - start);

   printf("sol: \n");
   for (int i=0; i<20; i++){
      printf("%1.2e\n", sol[i]);
   }
   printf("\n");

   FILE *f;

   // Export sol
   f = fopen("sol.txt", "w");
   for (int i = 0; i <= nn; ++i)
      fprintf(f, "%f\n", sol[i]);
   fclose(f);
   
   // free memory
   free(coord);
   free(coord_buff);
   free(tetra);
   free(tetra_buf);
   free(iat);
   free(ja);
   free(coefH);
   free(coefP);
   free(coefB);
   free(coefA);
   free(coefAmod);
   free(bdval);
   free(Abdval);
   free(sol);
   free(rhs);

}
