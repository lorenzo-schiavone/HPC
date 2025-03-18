#include <iostream>

struct csr{
    int N, NT;
    int* IAT;
    int* JA;
    double* coef;
};

int alloc_csr(csr &A, const int N, const int NT){
    A.N = N; A.NT = NT;
    A.IAT = (int*) malloc( (N+1) * sizeof(int));
    A.JA = (int*) malloc( (NT) * sizeof(int));
    A.coef = (double*) malloc( (NT) * sizeof(double));
    if (A.IAT == nullptr) return 1;
    if (A.JA == nullptr) return 2;
    if (A.coef == nullptr) return 3;
    return 0;
}

int dealloc_csr( csr &A){
    free(A.coef);
    free(A.IAT);
    free(A.JA);
    A.N = 0;
    A.NT=0;
    return 0;
}

csr sum_csr(const csr &A,const csr &B){
    // somma matrici sparse con lo stesso pattern
    csr C;
    int ierr = alloc_csr(C, A.N, A.NT);
    // if (ierr!=0){}
    double* coef_a = A.coef; // riferimenti diretti o handle
    double* coef_b = B.coef;
    double* coef_c = C.coef;

    for (int i = 0; i< A.NT; i++){
        coef_c[i] = coef_a[i]+ coef_b[i];
    }
    return C;
}

struct vec {
    int N;
    double* v;
    void print(){
        for (int j=0; j< N; j++){
            std::cout << v[j] << " ";
        }
        std::cout << std::endl;
    }   
};

int alloc_vec(vec &v){
    v.v = (double *) malloc(v.N*sizeof(double));
    if (v.v==nullptr) return 1;
    return 0;
}

double my_sum(vec v){
    double s = 0.0;
    for (int i=0; i< v.N; i++){
        s+=v.v[i];
    }
    return s;
}

struct mat {
    int N;
    int M;
    int** A;
    void print(){
        for (int i=0; i<N; i++){
            for (int j=0; j<M; j++){
            std::cout << A[i][j] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    } 
};

int alloc_mat(mat &A){
    A.A = (int **) malloc(A.M *sizeof(int*));
    int * Abuf;
    Abuf = (int*) malloc ( A.M*A.N * sizeof(int));
    if (Abuf==nullptr) return 1;
    for (int i=0;i<A.M;i++){
        A.A[i] = &Abuf[i*A.N];
    }
    return 0;
}

double scal_prod( const vec& v, const vec& w){
    double tot=0.;
    for (int i=0;i<v.N; i++){
        tot+= v.v[i]*w.v[i];
    }
    return tot;
}

vec matvecprod(const mat &A, const vec &v){
    vec w;
    w.N = A.M;
    alloc_vec(w);
    for (int i =0; i < A.N; i++){
        // crea handle
        w.v[i]=0.;
        for (int j =0; j< A.M; j++){
            w.v[i] += A.A[i][j] * v.v[j];
        }
    }
    return w;
}

mat matmatprod(const mat &A, const mat &B){
    mat C;
    C.M = A.M;
    C.N = B.N;
    alloc_mat(C);
    for (int i=0;i < C.M; i++){
        for (int j=0;j< C.N; j++){
            C.A[i][j]=0;
            for (int k=0; k<A.N; k++){
                C.A[i][j] += A.A[i][k]*B.A[k][j];
            }
        }
    }
    return C;
}

vec sparsemat_vec_prod(const csr &A,vec v){
    vec res;
    res.N = v.N;
    alloc_vec(res);
    for (int i=0;i<A.N; i++){
        res.v[i]=0.;
        int istart = A.IAT[i];
        int iend = A.IAT[i+1];
        for (int j=istart;j<iend;j++){
            res.v[i] += v.v[A.JA[j]]*A.coef[j];
        }
    }
    return res;
}

int main(){
    std::cout << "Hello World\n";
    vec v;
    v.N = 3;
    alloc_vec(v);
    v.v[0]=1.0;v.v[1]=1.0; v.v[2]=1.0;
    v.print();
    mat A;
    A.M = 3; A.N=3;
    alloc_mat(A);
    for (int i = 0; i < A.N; i++){
        for (int j = 0; j < A.M; j++)
            A.A[i][j] = 0;}
    A.A[0][0] = 1;
    A.A[0][2] = 1;
    A.A[1][0] = -1;
    A.A[2][1] = -1;

    vec w = matvecprod(A, v);

    A.print();

    w.print();

    csr B;
    int N =3;
    int NT = 4;
    alloc_csr(B, N, NT);
    B.coef[0] = 1;
    B.coef[1]= 1;
    B.coef[2]=-1;
    B.coef[3]=-1;
    B.JA[0]=0;
    B.JA[1]=1;
    B.JA[2]=0;
    B.JA[3]=1;
    B.IAT[0]=0;
    B.IAT[1]=2;
    B.IAT[2]=3;
    B.IAT[3]=4;
    
    vec z = sparsemat_vec_prod(B,v);
    z.print();

    return 0;
 }