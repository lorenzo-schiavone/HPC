#include <iostream>
#include <cstdio>
#include <omp.h> // 

int main(){
    omp_set_num_threads(4);
    int myid;
    #pragma omp parallel 
    {
        myid = omp_get_thread_num();
        // std::cout << "Hello World from processor " << myid << std::endl;
        fprintf(stdout, "Hello world from processor %d\n ", myid);
    }
    return 0;
}