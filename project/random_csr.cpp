#include <cstdlib>
#include <cassert>
#include <cstdio> 

void generate_stiffness_csr(int nrows, double*& coef, int*& ja, int*& iat) {
    assert(nrows >= 2);  // Ensure the matrix has at least 2 rows

    int nnz = 3 * nrows - 2;

    coef = (double*) malloc(nnz * sizeof(double));
    ja = (int*) malloc(nnz * sizeof(int));
    iat = (int*) malloc((nrows + 1) * sizeof(int));

    int idx = 0;
    for (int i = 0; i < nrows; ++i) {
        iat[i] = idx;

        if (i - 1 >= 0) {
            if (idx >= nnz) { fprintf(stderr, "Overflow at i=%d (left)\n", i); exit(1); }
            coef[idx] = -1.0;
            ja[idx] = i - 1;
            ++idx;
        }

        if (idx >= nnz) { fprintf(stderr, "Overflow at i=%d (diag)\n", i); exit(1); }
        coef[idx] = 4.0;
        ja[idx] = i;
        ++idx;

        if (i + 1 < nrows) {
            if (idx >= nnz) { fprintf(stderr, "Overflow at i=%d (right)\n", i); exit(1); }
            coef[idx] = -1.0;
            ja[idx] = i + 1;
            ++idx;
        }
    }

    iat[nrows] = idx;
}
