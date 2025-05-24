#include <cmath>
#include <cstdlib>
#include <fstream>

double norm(double* v, int n){
  double acc=0.;
  for (int i=0;i<n;i++){
    acc += v[i]*v[i];
  }
  return sqrt(acc);
}

double scalarprod(double* v, double* w, int n){
  double acc=0.;
  for (int i=0;i<n;i++){
    acc += v[i]*w[i];
  }
  return acc;
}
// QR FACTORIZATION FOR COLUMN BASED MATRICES
void qr(double** A, int nrow, int ncol, double*** Q_out, double*** R_out){
  double** Q = (double**) malloc(nrow * sizeof(double*));
  double* Qbuf = (double*) malloc(nrow*nrow * sizeof(double));
  for (int i=0; i<nrow;i++){
    Q[i] = &Qbuf[i*nrow];
  }
  double** R = (double**) malloc(ncol * sizeof(double*));
  double* Rbuf = (double*) malloc(nrow*ncol * sizeof(double));
  for (int i=0; i<ncol;i++){
    R[i] = &Rbuf[i*nrow];
  }
  // set R, Q ==0
  for (int j=0; j<nrow; j++){
     for (int i=0; i<nrow; i++){
      Q[j][i]= 0.;
     }
  }
  for (int j=0; j<ncol; j++){
     for (int i=0; i<nrow; i++){
      R[j][i]= 0.;
     }
  }
  // copy column of A into Q -- if column left put one in the diagonal
  int min_d;
  if (nrow<ncol){
    min_d = nrow;
  }
  else{
    min_d = ncol;
  }
  printf("nrow: %u\n", nrow);
  printf("ncol: %u\n", ncol);


  printf("min_d = %u\n", min_d);
  for (int j=0; j<min_d; j++){
     for (int i=0; i<nrow; i++){
      Q[j][i]= A[j][i];
     }
  }
  if (min_d< nrow) {
    for (int i=min_d; i< nrow; i++){
      Q[i][i] =1.;
    }
  }
  // here begin mgs
  double alpha;
  for (int i=0;i<min_d;i++){
    // normalizzo colonna Q[i]
    R[i][i] = norm(Q[i], nrow);
    if (R[i][i] < 1e-10){
      break; // è una colonna di zeri non bene
    }
    for (int j=0;j<nrow;j++){
        Q[i][j]/=R[i][i];
    }

    // rimuovo componenti parallele a Q[i] da Q[j] per j> i
    
    for (int j=i+1;j<nrow; j++){
      alpha = scalarprod(Q[i], Q[j], nrow);
      for (int jj=0;jj<nrow;jj++){
              Q[j][jj]-=(alpha*Q[i][jj]);
          }
      if ((j<ncol) && (i<nrow)){
         R[j][i] = alpha;
      }
     
    }
  }
  
  for (int i = min_d; i<nrow; i++){
    printf("Q[%u]: \n", i);
    for (int j=0; j<nrow; j++){
        printf("%f ", Q[i][j]);

    }
    printf("\n");
    
    alpha = norm(Q[i], nrow);
    printf("alpha: %f\n", alpha);
    for (int j=0;j<nrow;j++){
        Q[i][j]/=alpha;
    }

  }
  *Q_out = Q;
  *R_out = R;
}