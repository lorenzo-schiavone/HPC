#include <iostream>

struct vec {
    int N;
    double* v;
    void print() {
        for (int j = 0; j < N; j++) {
            std::cout << v[j] << " ";
        }
        std::cout << std::endl;
    }   
};

int alloc_vec(vec &v) {
    v.v = (double *)malloc(v.N * sizeof(double));
    if (v.v == nullptr) return 1;
    return 0;
}

int main() {
    vec v;
    v.N = 3;
    int alloc_result = alloc_vec(v);
    if (alloc_result != 0) {
        std::cerr << "Memory allocation failed with code " << alloc_result << std::endl;
        return alloc_result;
    }
    v.v[0] = 1.0;
    v.v[1] = 1.0;
    v.v[2] = 1.0;
    v.print();
    free(v.v); // Don't forget to free the allocated memory
    return 0;
}