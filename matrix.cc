#include<iostream>
#include<sys/time.h>
#include<cblas.h>
#include<vector>
#include<unistd.h>
#include<cstdlib>

using namespace std;

template<class T>
void MatrixPrint(vector< vector<T> > a, size_t cols, size_t rows){
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++){
            cout << a[i][j] << " ";
        }
        cout << endl;
    }
}

long diff(struct timeval* b, struct timeval* e){
    return (e->tv_sec - b->tv_sec) * 1000000 + e->tv_usec - b->tv_usec;
}


void RandomMatrix(double* m, const int mSize)
{
    srand(time(NULL));
    for (int i = 0; i < mSize; i++)
        m[i] = 1.0 / (double)(i);
}


int main(){

    const int N = 1024;//size of matrix, to simplify, we use square matrix
    const int n = 100;//number  of test time
    struct timeval tv;
    struct timeval etv;
    // const double *a = (double*) malloc(sizeof(double) * N * N);
    // const double *b = (double*) malloc(sizeof(double) * N * N);
    // const double *c = (double*) malloc(sizeof(double) * N * N);

    double* a = new double[N*N];
    double* b = new double[N*N];
    double* res = new double[N*N];
    RandomMatrix(a, N * N);
    RandomMatrix(b, N * N);
    //memset(res, 0, sizeof(double) * N * N);

    for(int i = 0; i < 10; i++){
        double tmp = 0;
        for(int h = 0; h < N; h++){
            for(int k = 0; k < N; k++){
                for(int j  = 0; j < N; j++){
                    tmp +=  *(a + h * N + j) *  *(b + j * N + k);
                }
                *(res + N * h + k) = tmp;
            }
        }
        // Change the matrix to avoid from being "optimized"
        a[0] += 1.0;
        b[0] += 2.0;
    }

    gettimeofday(&tv, NULL);
    for(int i = 0; i < n; i++){
        double tmp = 0;
        for(int h = 0; h < N; h++){
            for(int j  = 0; j < N; j++){
                for(int k = 0; k < N; k++){
                
                    tmp +=  *(a + h * N + j) *  *(b + j * N + k);
                }
                *(res + N * h + k) = tmp;
            }
        }
        // Change the matrix to avoid from being "optimized"
        a[0] += 1.0;
        b[0] += 2.0;
    }
    gettimeofday(&etv, NULL);
    std::cout << "matrix multiply my implementation time : " << diff(&tv, &etv) <<endl;


    gettimeofday(&tv, NULL);
    for(int i = 0; i < n; i++){
        
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1, a, N, b, N, 0, res, N);
        // Change the matrix to avoid from being "optimized"
        a[0] += 1.0;
        b[0] += 2.0;
    }
    gettimeofday(&etv, NULL);
    std::cout << "matrix multiply openblas time : " << diff(&tv, &etv)<<endl ;
    free(a);
    free(b);
    free(res);

    return 0;
    
}