#include <cmath>
#include <cstdlib>

double scalprod(double* v, double* w, int k) {
    double acc = 0.0;
    for (int i = 0; i < k; ++i) {
        acc += v[i] * w[i];
    }
    return acc;
}

double norm(double* v, int k) {
    return std::sqrt(scalprod(v, v, k));
}

// A,Q,R are pointer to columns here i don't know if that may be a problem
void qr(double** A, double** Q, double** R, int m) {
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            Q[j][i] = A[j][i]; 
        }
    }
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            R[i][j] = 0.0;
        }
    }

    for (int j = 0; j < m; ++j) {
        R[j][j] = norm(Q[j], m);
        for (int i = 0; i < m; ++i) {
            Q[j][i] /= R[j][j];
        }
        
        for (int k = j + 1; k < m; ++k) {
            R[j][k] = scalprod(Q[j], Q[k], m);
            for (int i = 0; i < m; ++i) {
                Q[k][i] -= R[j][k] * Q[j][i];
            }
        }
    }
}