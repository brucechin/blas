#include<iostream>
#include<sys/time.h>
#include<cblas.h>
#include<vector>
#include<unistd.h>
#include<cstdlib>

using namespace std;
#define MAX 100

long diff(struct timeval* b, struct timeval* e){
    return (e->tv_sec - b->tv_sec) * 1000000 + e->tv_usec - b->tv_usec;
}


template<class T>
void RandomFill(std::vector<T>& numbers, size_t size){
    numbers.resize(size);
    for(size_t i = 0; i != size; i++){
        numbers[i] = static_cast<T>(rand() % 5000);
        numbers[i] *= 0.1;
    }
}

template<class T>
void MatrixRandomFill(std::vector<std::vector<T> >& numbers, size_t cols, size_t rows){
    vector<T> tmp(cols);
    numbers.resize(rows, tmp);
    for(size_t i = 0; i != rows; i++){
        for(size_t j = 0; j != cols; j++){
            numbers[i][j] = static_cast<T>(rand() % 5000);
            numbers[i][j] *= 0.1;
        }
    }
}


template<class T>
T Dot(T* a, T* b, size_t size){
    //vector dot vector
    T sum = 0.0;
    for(int i = 0; i < size; i++){
        sum += (a[i] * b[i]);
    }
    return sum;
}



void Level1Test(int tests, int size)
{
    //vector multiply vector test

    int number_of_tests = tests;
    int size_of_vector = size;

    double res = 0.0;
    struct timeval tv;
    struct timeval etv;

//double precision input
    while(size_of_vector < 100001){

        vector<double> da, db;
        vector<float> fa, fb;
        RandomFill(da, size_of_vector);
        RandomFill(db, size_of_vector);
        RandomFill(fa, size_of_vector);
        RandomFill(fb, size_of_vector);

        for(int i = 0; i < 10; i++){
            Dot(da.data(), db.data(), da.size());
        }    

        gettimeofday(&tv, NULL);

        for(int i = 0; i < number_of_tests; i++){
            res = Dot(da.data(), db.data(), da.size());
        }    
        gettimeofday(&etv, NULL);

        std::cout << "my double precision dot time : " << diff(&tv, &etv);
        printf(" res : %f\n", res);


        res = 0.0;
        gettimeofday(&tv, NULL);
        for(int i = 0; i < number_of_tests; i++){
            res = cblas_ddot(da.size(), da.data(), 1, db.data(), 1);
        }    
        gettimeofday(&etv, NULL);
        std::cout << "openblas ddot time: " << diff(&tv, &etv) ;
        printf(" res : %f\n", res);

    // single precision input


        for(int i = 0; i < 10; i++){
            Dot(fa.data(), fb.data(), fa.size());
        }    

        gettimeofday(&tv, NULL);

        for(int i = 0; i < number_of_tests; i++){
            res = Dot(fa.data(), fb.data(), fa.size());
        }    
        gettimeofday(&etv, NULL);

        std::cout << "my single precision dot time : " << diff(&tv, &etv) ;
        printf(" res : %f\n", res);

        res = 0.0;

        gettimeofday(&tv, NULL);
        for(int i = 0; i < number_of_tests; i++){
            res = cblas_sdot(fa.size(), fa.data(), 1, fb.data(), 1);
        }    
        gettimeofday(&etv, NULL);

        std::cout << "openblas sdot time : " << diff(&tv, &etv) ;
        printf(" res : %f\n", res);



        res = 0.0;

        gettimeofday(&tv, NULL);
        for(int i = 0; i < number_of_tests; i++){
            res = cblas_dsdot(fa.size(), fa.data(), 1, fb.data(), 1);
        }    
        gettimeofday(&etv, NULL);

        std::cout << "openblas dsdot time : " << diff(&tv, &etv) ;
        printf(" res : %f\n", res);




        res = 0.0;

        gettimeofday(&tv, NULL);
        for(int i = 0; i < number_of_tests; i++){
            res = cblas_sdsdot(fa.size(), 0.0, fa.data(), 1, fb.data(), 1);
        }    
        gettimeofday(&etv, NULL);

        std::cout << "openblas sdsdot: " << diff(&tv, &etv) << std::endl;
        printf(" res : %f\n", res);
        
        size_of_vector *= 2;    
    }    
}

template<class T>
void MatrixPrint(vector< vector<T> > a, size_t cols, size_t rows){
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++){
            cout << a[i][j] << " ";
        }
        cout << endl;
    }
}


template<class T>
T* MatrixDot(T a[MAX][MAX], T b[MAX][MAX], T res[MAX][MAX], const size_t H, const size_t J,const size_t K){
    //matrix multiply matrix
    //matrix a size is h*j, matrix b size is j*k

    T tmp = 0.0;

    for(int h = 0; h < H; h++){

        for(int k = 0; k < K; k++){
            for(int j  = 0; j < J; j++){
                tmp += a[h][j] * b[j][k];
            }
            res[h][k] = tmp;
        }

    }

    return res;
}

void Level2Test(int tests, size_t cols, size_t rows, size_t t){
    //matrix multiply vector test

    struct timeval tv;
    struct timeval etv;
    vector< vector<double> > da, db, res;
    double arr_a[cols][rows];
    double arr_b[t][cols];
    double arr_res[t][rows];
    MatrixRandomFill(da, cols, rows);
    MatrixRandomFill(db, t, cols);
    MatrixRandomFill(res, t, rows);
    

    int number_of_tests = tests;



    
    //MatrixPrint(res, 1, rows);

    // for(int i = 0; i < 10; i++){
    //     //warm up
    //     res = MatrixDot(da, db, res, rows, cols, t);
    // }    

    // gettimeofday(&tv, NULL);
    // for(int i = 0; i < number_of_tests; i++){
    //     res = MatrixDot(da, db, res, rows, cols, t);
    // }    
    // gettimeofday(&etv, NULL);
    // std::cout << "my implementation dot time : " << diff(&tv, &etv);

/*

void cblas_dgemm (const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb, const MKL_INT m, const MKL_INT n, 
const MKL_INT k, const double alpha, const double *a, const MKL_INT lda, const double *b, const MKL_INT ldb, const double beta, double *c, const MKL_INT ldc);
*/


    gettimeofday(&tv, NULL);
    for(int i = 0; i < number_of_tests; i++){
        // double* a = &da[0][0];
        // double* b = &db[0][0];
        // double* r = &res[0][0];
        //cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, rows, cols, t, 1, a, cols, b, t, 0, r, t);
        printf("ok");
    }    
    gettimeofday(&etv, NULL);

    std::cout << "openblas dgemm time : " << diff(&tv, &etv) ;

}

int main(){
    

    Level2Test(10,10,10,10);

}