#include <iostream>
// #include <omp.h> // 
#include "/usr/local/opt/libomp/include/omp.h"

int main(){
    omp_set_num_threads(4);
    #pragma omp parallel
    {
        int myid = omp_get_thread_num();
        std::cout << "Hello world from processor " << myid << std::endl;
    }
    return 0;
}