#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>

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

  // if nrow > ncol: Q is nrow x nrow, R is nrow x ncol 
  // this is the only case i am interested in, or ncol = nrow

  double** Q = (double**) malloc(nrow * sizeof(double*));
  double* Qbuf = (double*) calloc(nrow*nrow, sizeof(double));
  for (int i=0; i<nrow;i++){
    Q[i] = &Qbuf[i*nrow];
  }

  double** R = (double**) malloc(ncol * sizeof(double*));
  double* Rbuf = (double*) calloc(nrow*ncol,sizeof(double));
  for (int i=0; i<ncol;i++){
    R[i] = &Rbuf[i*nrow];
  }

  int min_d = (nrow<ncol) ? nrow : ncol;

  // copy columns of A into Q 
  for (int j=0; j<min_d; j++){
     for (int i=0; i<nrow; i++){
      Q[j][i]= A[j][i];
     }
  }
  // if column left put random
  if (min_d< nrow) {
    for (int i=min_d; i< nrow; i++){
        for (int j=0; j<nrow; j++){
            Q[i][j] = (double)rand() / RAND_MAX;
        }
    }
  }
  // here begin mgs
  double alpha;
  for (int i=0;i<min_d;i++){
    // normalizzo colonna Q[i]
    alpha= norm(Q[i], nrow);
    for (int j=0;j<nrow;j++){
        Q[i][j]/=alpha;
    }
    R[i][i] = alpha;

    // rimuovo componenti parallele a Q[i] da Q[j] per j> i
    for (int j=i+1;j<ncol; j++){
      alpha = scalarprod(Q[i], Q[j], nrow);
      for (int jj=0;jj<nrow;jj++){
        Q[j][jj]-=(alpha*Q[i][jj]);
      }
      R[j][i] = alpha;
    }
    for (int j=ncol;j<nrow; j++){
      alpha = scalarprod(Q[i], Q[j], nrow);
      for (int jj=0;jj<nrow;jj++){
        Q[j][jj]-=(alpha*Q[i][jj]);
      }
    }
  }
  // completo base ortonormale di Q -> non devo scrivere su R 
  for (int i = min_d; i<nrow; i++){
    alpha = norm(Q[i], nrow);
    for (int j=0;j<nrow;j++){
      Q[i][j]/=alpha;
    }
    for (int j=i+1;j<nrow; j++){
      alpha = scalarprod(Q[i], Q[j], nrow);
      for (int jj=0;jj<nrow;jj++){
        Q[j][jj]-=(alpha*Q[i][jj]);
      }
    }

  } 
  *Q_out = Q;
  *R_out = R;
}