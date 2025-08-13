#include <stdlib.h> 
#include <stdio.h> 

void topol(const int n, const int nt, const int n1, const int * const *tetra, int &nterm,
           int *&iat, int *&ja){

int I1[4];
        
// Overestimate number of non-zeroes
nterm = n1*n;

// Init ja_tmp
int *ja_tmp = (int*) malloc(nterm*sizeof(int));
for (int i = 0; i< nterm; i++) ja_tmp[i] = -1;

for (int k = 0; k < nterm; k++){
   if ( (k/n1)*n1 == k ) ja_tmp[k]=k/n1;
}

// Analize all tetrahedrons ( nt )
for (int k = 0; k < nt; k++){

   for (int i = 0; i < 4; i++) I1[i] = tetra[k][i];
    
   // Sort tet nodes
   for (int j = 0; j < 3; j++){
      for (int i = 0; i < 3 - j; i++){
         if ( I1[i] > I1[i+1] ){
            int aux=I1[i];
            I1[i]=I1[i+1];
            I1[i+1]=aux;
         }
      }
   }
    
   // Build ja_tmp array
   for (int i = 0; i < 3; i++){
      int j = i+1;
      for (int l = j; l < 4; l++){
         int m = n1*I1[i]+l-i;
         int MCONTR = n1*(I1[i]+1);
         bool continua = true;
         while ( ((I1[l]-ja_tmp[m]) < 0) || ((I1[l]-ja_tmp[m]) > 0) ){
            if ((I1[l]-ja_tmp[m]) > 0){
               if (ja_tmp[m] == -1 ){
                  ja_tmp[m] = I1[l];
               } else {
                  m++;
                  if ( (m-MCONTR)>=0 ){
                     printf("1) Error for tetra %d unknown %d\n",k,I1[l]);
                     printf("Increase n1\n");
                     int ind = n1*I1[i];
                     for (int ii = ind; ii < ind+n1; ii++) printf("%d\n",ja_tmp[ii]);
                     fflush(stdin);
                     return;
                  }
               }
            } else if ( (I1[l]-ja_tmp[m]) < 0 ){
               int mm = m + 1;
               while (ja_tmp[mm] > -1){
                  mm++;
                  if ( m-MCONTR>=0){
                     printf("2) Error for tetra %d unknown %d\n",k,I1[l]);
                     printf("Increase n1\n");
                     int ind = n1*I1[i];
                     for (int ii = ind; ii < ind+n1; ii++) printf("%d\n",ja_tmp[ii]);
                     fflush(stdin);
                     return;
                  }
               }
               while (ja_tmp[mm] == -1 ){
                 ja_tmp[mm] = ja_tmp[mm-1];
                 mm--;
                 while ( mm > m ){
                    ja_tmp[mm] = ja_tmp[mm-1];
                    mm--;
                 }
               }
               ja_tmp[m]=I1[l];
            }
         } // while loop
      } // l loop

   } // i loop

} // k loop

// Build iat array
int *iat_tmp = (int*) malloc((n+1)*sizeof(int));
iat_tmp[0] = 0;
int m = 0;
int j = 0;
for (int k = 0; k < nterm; k += n1){
   for (int i = 0; i < n1; i++){
      if ( ja_tmp[k+i] >= 0){
         m++;
      } 
   }
   j++;
   iat_tmp[j] = m;
}

// Compact ja_tmp
m = 0;
for (int k = 0; k < nterm; k++){
   if (ja_tmp[k] >= 0){
      ja_tmp[m] = ja_tmp[k];
      m++;
   }
}

// Complete the matrix
nterm = 2*m-n;
iat = (int*) malloc((n+1)*sizeof(int));
ja = (int*) malloc(nterm*sizeof(int));
int *ptrow = (int*) malloc(n*sizeof(int));
for (int i = 0; i < n; i++) ptrow[i] = 0;
for (int i = 0; i < n; i++){
   for (int j = iat_tmp[i]; j < iat_tmp[i+1]; j++){
      int jcol = ja_tmp[j];
      ptrow[i]++;
      if (jcol != i) ptrow[jcol]++;
   }
}
iat[0] = 0;
for (int i = 0; i < n; i++){
   iat[i+1] = iat[i] + ptrow[i];
   ptrow[i] = iat[i];
}

for (int i = 0; i < n; i++){
   for (int j = iat_tmp[i]; j < iat_tmp[i+1]; j++){
      int jcol = ja_tmp[j];
      ja[ptrow[i]++] = jcol;
      if (i != jcol) ja[ptrow[jcol]++] = i;
   }
}
free(iat_tmp);
free(ja_tmp);
free(ptrow);

}
