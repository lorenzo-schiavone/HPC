g++ -Xpreprocessor -fopenmp topol.cpp main.cpp -L/usr/local/opt/libomp/lib -lomp -I/usr/local/opt/libomp/include -o vol

g++ -Xpreprocessor -fopenmp qr.cpp gmres.cpp random_csr.cpp gmres_with_test.cpp -L/usr/local/opt/libomp/lib -lomp -I/usr/local/opt/libomp/include -o gmrestest