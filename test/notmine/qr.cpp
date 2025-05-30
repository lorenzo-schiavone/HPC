#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "my_blas.h"

// Compute the QR factorization of A = Q * R using classical Gram–Schmidt.
// A is of size (nr x nc), Q will be (nr x nc) with orthonormal columns,
// and R will be (nc x nc) upper triangular.
void qr(double** A, double** Q, double** R, int nr, int nc) {
    for (int j = 0; j < nc; j++) {
        // Copy the j-th column of A into a temporary vector v.
        double* v = new double[nr];
        for (int i = 0; i < nr; i++) {
            v[i] = A[i][j];
        }
        // For each previous column, subtract its component.
        for (int i = 0; i < j; i++) {
            double dot = 0.0;
            // Compute the dot product between Q(:,i) and A(:,j).
            for (int k = 0; k < nr; k++) {
                dot += Q[k][i] * A[k][j];
            }
            R[i][j] = dot;
            // Subtract the projection: v = v - dot * Q(:,i)
            for (int k = 0; k < nr; k++) {
                v[k] -= dot * Q[k][i];
            }
        }
        // Compute the norm of the modified vector v.
        double norm = 0.0;
        for (int i = 0; i < nr; i++) {
            norm += v[i] * v[i];
        }
        norm = sqrt(norm);
        R[j][j] = norm;
        // Normalize v to form the j-th column of Q.
        for (int i = 0; i < nr; i++) {
            Q[i][j] = v[i] / norm;
        }
        delete[] v;
    }
}

int main() {
    // Define the dimensions of A.
    // For example, try a rectangular matrix with more rows than columns.
    int nr = 11; // number of rows
    int nc = 10;  // number of columns

    // Allocate A (nr x nc)
    double** A = new double*[nr];
    A[0] = new double[nr * nc];
    for (int i = 1; i < nr; i++) {
        A[i] = A[i - 1] + nc;
    }

    // Initialize A with random values.
    for (int i = 0; i < nr; i++) {
        for (int j = 0; j < nc; j++) {
            A[i][j] = (double)rand() / RAND_MAX;
        }
    }

    // Allocate Q (nr x nc)
    double** Q = new double*[nr];
    Q[0] = new double[nr * nc];
    for (int i = 1; i < nr; i++) {
        Q[i] = Q[i - 1] + nc;
    }

    // Allocate R (nc x nc) (upper triangular)
    double** R = new double*[nc];
    R[0] = new double[nc * nc];
    for (int i = 1; i < nc; i++) {
        R[i] = R[i - 1] + nc;
    }

    // Compute the QR factorization of A.
    qr(A, Q, R, nr, nc);

    // Verify the factorization: compute C = Q * R, which should approximate A.
    // Allocate C (nr x nc)
    double** C = new double*[nr];
    C[0] = new double[nr * nc];
    for (int i = 1; i < nr; i++) {
        C[i] = C[i - 1] + nc;
    }
    // Multiply Q (nr x nc) by R (nc x nc) to obtain C (nr x nc).
    matmatprod(Q, R, C, nr, nc, nc);

    // Print the verification: (C - A) should be (nearly) zero.
    printf("Verification (C - A):\n");
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++) {
            printf("% .2e ", C[i][j] - A[i][j]);
        }
        printf("\n");
    }

    // Cleanup memory.
    delete[] A[0]; delete[] A;
    delete[] Q[0]; delete[] Q;
    delete[] R[0]; delete[] R;
    delete[] C[0]; delete[] C;

    return 0;
}
