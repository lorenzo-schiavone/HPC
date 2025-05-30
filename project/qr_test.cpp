#include <cmath>
#include <cstdlib>
#include <fstream>

#include "qr.h"

// double norm(double* v, int n){
//   double acc=0.;
//   for (int i=0;i<n;i++){
//     acc += v[i]*v[i];
//   }
//   return sqrt(acc);
// }

// double scalarprod(double* v, double* w, int n){
//   double acc=0.;
//   for (int i=0;i<n;i++){
//     acc += v[i]*w[i];
//   }
//   return acc;
// }

void matrixmatrixprod( double** A, double** B, double*** C_out, int nrowA, int ncolA, int ncolB){
    double** C = (double**) malloc(ncolB * sizeof(double*));
    double* Cbuf = (double*) malloc(ncolB*nrowA * sizeof(double));
    for (int i=0; i<ncolB;i++){
        C[i] = &Cbuf[i*nrowA];
        for (int j=0;j<nrowA;j++){
            C[i][j] = 0.;
        }
    }
    double red = 0.;
    for (int j=0;j<ncolB;j++){
        for (int i=0;i<nrowA;i++){
            red=0;
            for (int k=0; k<ncolA;k++){
                red += A[k][i]*B[j][k]; // C_ij = A_ik * B_kj
            }
            C[j][i] = red;
        }
    }
    *C_out = C;
}


// void qr(double** A, int nrow, int ncol, double*** Q_out, double*** R_out){
//   double** Q = (double**) malloc(nrow * sizeof(double*));
//   double* Qbuf = (double*) malloc(nrow*nrow * sizeof(double));
//   for (int i=0; i<nrow;i++){
//     Q[i] = &Qbuf[i*nrow];
//   }
//   double** R = (double**) malloc(ncol * sizeof(double*));
//   double* Rbuf = (double*) malloc(nrow*ncol * sizeof(double));
//   for (int i=0; i<ncol;i++){
//     R[i] = &Rbuf[i*nrow];
//   }
//   // set R, Q ==0
//   for (int j=0; j<nrow; j++){
//      for (int i=0; i<nrow; i++){
//       Q[j][i]= 0.;
//      }
//   }
//   for (int j=0; j<ncol; j++){
//      for (int i=0; i<nrow; i++){
//       R[j][i]= 0.;
//      }
//   }
//   // copy column of A into Q -- if column left put one in the diagonal
//   int min_d = (nrow<ncol ? nrow : ncol);
//   for (int j=0; j<min_d; j++){
//      for (int i=0; i<nrow; i++){
//       Q[j][i]= A[j][i];
//      }
//   }
//   if (min_d< nrow) {
//     for (int i=min_d; i< nrow; i++){
//       Q[i][i] =1.;
//     }
//   }
//   // here begin mgs
//   double alpha;
//   for (int i=0;i<min_d;i++){
//     // normalizzo colonna Q[i]
//     R[i][i] = norm(Q[i], nrow);
//     if (R[i][i] < 1e-10){
//       continue; // è una colonna di zeri non bene
//     }
//     for (int j=0;j<nrow;j++){
//         Q[i][j]/=R[i][i];
//     }

//     // rimuovo componenti parallele a Q[i] da Q[j] per j> i
    
//     for (int j=i+1;j<nrow; j++){
//       alpha = scalarprod(Q[i], Q[j], nrow);
//       for (int jj=0;jj<nrow;jj++){
//               Q[j][jj]-=(alpha*Q[i][jj]);
//           }
//       if ((j<ncol) && (i<nrow)){
//          R[j][i] = alpha;
//       }
     
//     }
//   }
//   for (int i = min_d; i<nrow; i++){
//     alpha = norm(Q[i], nrow);
//     for (int j=0;j<nrow;j++){
//         Q[i][j]/=alpha;
//     }

//   }
//   *Q_out = Q;
//   *R_out = R;
// }


/*
come usarlo
double** Q;
double** R;
qr(H, m+1, m, &Q, &R);

Q_H[m][0] dovrà essere usato per il residuo

*/

void printMatrix(double** mat, int nrows, int ncols, const char* name) {
    printf("%s:\n", name);
    for (int i = 0; i < nrows; ++i) {
        for (int j = 0; j < ncols; ++j) {
            printf("%2.2f ", mat[j][i]);
        }
        printf("\n");
    }
}

void testMatrixMatrixProduct() {
    printf("Testing matrixmatrixprod...\n");
    // A 2x3
    int nrow = 2, ncol = 3;
    double** A = (double**)malloc(ncol * sizeof(double*));
    double* Abuf = (double*)malloc(ncol * nrow * sizeof(double));
    for (int i = 0; i < ncol; ++i) {
        A[i] = &Abuf[i * nrow];
        for (int j = 0; j < nrow; ++j) {
            A[i][j] = 0.0;
        }
    }
    A[0][0]=1.;
    A[2][1]=1;
    A[2][0]=1.;
    A[1][1]=.5;
    // A: 
    // [[1 0 1]
    //  [0 .5 1]

    // B 3x2
    double** B = (double**)malloc(nrow * sizeof(double*));
    double* Bbuf = (double*)malloc(ncol * nrow * sizeof(double));
    for (int i = 0; i < nrow; ++i) {
        B[i] = &Bbuf[i * ncol];
        for (int j = 0; j < ncol; ++j) {
            B[i][j] = 0.;
        }
    }
    B[0][0]=1.; B[1][0]=1; B[1][1]=1.; B[0][2] = .5;
    // B: 
    // [[1  1]
    //  [0  1]
    //  [.5 0]
    printMatrix(A, nrow, ncol, "A");
    printMatrix(B, ncol, nrow, "B");

    double** C;
    matrixmatrixprod(A, B, &C, nrow, ncol, nrow);
    printMatrix(C, nrow, nrow, "C = AB");

    free(A[0]); free(A);
    free(B[0]); free(B);
    free(C[0]); free(C);

}

void testQRDecomposition() {
    printf("\nTesting QR decomposition...\n");
    // Test1
    int nrow = 24, ncol = 20;
    double** A = (double**)malloc(ncol * sizeof(double*));
    double* Abuf = (double*)malloc(ncol * nrow * sizeof(double));
    for (int i = 0; i < ncol; ++i) {
        A[i] = &Abuf[i * nrow];
        for (int j = 0; j < nrow; ++j) {
            A[i][j] = (double)rand() / RAND_MAX;
        }
    }

    printMatrix(A, nrow, ncol, "A");
    double** Q, ** R;
    qr(A, nrow, ncol, &Q, &R);
    printMatrix(Q, nrow, nrow, "Q");
    printMatrix(R, nrow, ncol, "R");

     // Check QR = A
    double** CRes;
    matrixmatrixprod(Q, R, &CRes, nrow, nrow, ncol);
    for (int i = 0; i < ncol; ++i) {
        for (int j = 0; j < nrow; ++j) {
            CRes[i][j] -= A[i][j];
        }
    }
    printMatrix(CRes, nrow, ncol, "QR-A");

    // Check Q^TQ = I
    double** QT = (double**)malloc(nrow * sizeof(double*));
    double* QTbuf = (double*)malloc(nrow * nrow * sizeof(double));
    for (int i = 0; i < nrow; ++i) {
        QT[i] = &QTbuf[i * nrow];
        for (int j = 0; j < nrow; ++j) {
            QT[i][j] =Q[j][i];
        }
    }
    double** QTQ;
    matrixmatrixprod(QT, Q, &QTQ, nrow, nrow, nrow);
    for (int i = 0; i < nrow; ++i) {
      QTQ[i][i] -= 1;
    }
    printMatrix(QTQ, nrow, nrow, "Q^T * Q-I");

    free(A[0]); free(A);
    free(Q[0]); free(Q);
    free(R[0]); free(R);
    free(QTQ[0]); free(QTQ);
}

int main() {
    // testMatrixMatrixProduct();
    testQRDecomposition();
    return 0;
}