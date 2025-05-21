#include <iostream> 
#include <stdlib.h>
#include <fstream>
#include <cstdlib>
using namespace std;

#include "topol.h"

double det3(double* r0,double* r1,double* r2 ){
   return r0[0]*(r1[1]*r2[2]-r2[1]*r1[2]) 
        - r1[0]*(r0[1]*r2[2]-r2[1]*r0[2])
        + r2[0]*(r0[1]*r1[2]-r1[1]*r0[2]);
}
int sign(double x){
   return (x>0) - (x<0);
}

// MAIN PROGRAM
int main(int argc, const char* argv[]){

   // Check arguments
   if (argc < 2) {
      printf("Too few arguments.\n Usage: [tetra file] [coord file]\n");
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
   topol(nn,ntet,30,tetra,nterm,iat,ja); // is 30 enough? connectivity per node is ok, it would be very irregular
   printf("Topology created!\n");
   // here we have iat and ja ready to be used
   // -------------------------------------------------------------------------------------------------------------------------------

   // READ COORD
   ifstream coord_file(argv[2]);
   // Read header (i have added)
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

   // Set C-style of coord
   for (int i = 0; i < nnodes; i++)
      for (int j = 0; j < 3; j++) coord[i][j]--;

   // -------------------------------------------------------------------------------------------------------------------------------

   // diffusion and flow velocity
   double* D = (double*) malloc(3*sizeof(double));
   D[0]=1.;D[1]=1.;D[2]=1.;
   double* v = (double*) malloc(3*sizeof(double));
   v[0]=1.;v[1]=1.;v[2]=1.;

   double* volumes = (double*) malloc(ntet*sizeof(double));
   
   for (int i=0;i<ntet; i++){
      // volumes[i] = 1/6 * det[ 1....] det3x3
      int* tet = tetra[i];
      double* n0=coord[tet[0]];
      double* n1=coord[tet[1]];
      double* n2=coord[tet[2]];
      double* n3=coord[tet[3]];

      volumes[i] = (det3(n1,n2,n3)-
                   det3(n0,n2,n3)+
                   det3(n0,n1,n3)-
                   det3(n0,n1,n2))/6;
   }
   printf("Volumes: \n");
   for (int i=0;i<10; i++){
      printf("%f\n", volumes[i]);
   }
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

   for (int jj=0;jj<ntet;jj++){
      int* tet = tetra[jj];
      double vol = volumes[jj];
      double a[4];
      double b[4];
      double c[4];
      double d[4];

      for (int ii=0;ii<4;ii++){
         double* nj=coord[tet[(ii+1)%4]];
         double* nk=coord[tet[(ii+2)%4]];
         double* nm=coord[tet[(ii+3)%4]];
         a[ii] = det3(nj,nk,nm);
         b[ii] = - (nk[1]*nm[2]-nm[1]*nk[2] - (nj[1]*nm[2]-nj[2]*nm[1])+nj[1]*nk[2]-nj[2]*nk[1]);
         c[ii] = (nk[0]*nm[2]-nm[0]*nk[2] - (nj[0]*nm[2]-nj[2]*nm[0])+nj[0]*nk[2]-nj[2]*nk[0]);
         d[ii] = (nk[1]*nm[0]-nm[1]*nk[0] - (nj[1]*nm[0]-nj[0]*nm[1])+nj[0]*nk[2]-nj[0]*nk[1]); // I swap col 1 with col 2 so the det change sign
      }

      double *Hloc_buf = (double*) malloc(4*4*sizeof(double));
      double **Hloc = (double**) malloc(4*sizeof(double*));
      double *Ploc_buf = (double*) malloc(4*4*sizeof(double));
      double **Ploc = (double**) malloc(4*sizeof(double*));
      double *Bloc_buf = (double*) malloc(4*4*sizeof(double));
      double **Bloc = (double**) malloc(4*sizeof(double*));
      k=0;
      for (int i=0;i<4;i++){
         Hloc[i] = &(Hloc_buf[k]);
         Ploc[i] = &(Ploc_buf[k]);
         Bloc[i] = &(Bloc_buf[k]);
         k+=4;
      }

      for (int i=0;i<4; i++){
         double* Hi = Hloc[i];
         double* Pi = Hloc[i];
         double* Bi = Hloc[i];
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
            for (int ii=iat[tet[i]]; ii<iat[tet[i]+1]; ii++){
               int jjj = tet[j];
               if (ja[ii]== jjj){
                  coefH[ii] += Hi[j];
                  coefP[ii] += Pi[j];
                  coefB[ii] += Bi[j];
               }
            }
         }
      }
      
   }

   

   // before we may need to compute signed volume -> determinant 3x3
   // even before the coor matrix

   // Free memory
   free(coord);
   free(coord_buff);
   free(tetra);
   free(tetra_buf);
   free(iat);
   free(ja);

}
