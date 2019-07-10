/*
 * @Author: lianke.qin@gmail.com 
 * @Date: 2019-07-09 10:18:18 
 * @Last Modified by: lianke.qin@gmail.com
 * @Last Modified time: 2019-07-09 15:39:11
 */


#include"MatrixCalculator.h"
#include"Matrix.h"

//cblas_daxpy(const MKL_INT n, const double a, const double *x, const MKL_INT incx, double *y, const MKL_INT incy);
static Matrix* MatrixCalculator::add(Matrix mat1, Matrix mat2){
    int nrow = mat1.getNRow();
    int ncol = mat1.getNCol();
    Matrix* res = new Matrix(mat2);
    
    cblas_daxpy(nrow * ncol, 1.0, mat1.value, 1, res->value, 1);
    
    return res;
}

static Matrix* MatrixCalculator::add(Matrix mat1, Matrix mat2, int num){
    int nrow = mat1.getNRow();
    int ncol = mat1.getNCol();
    Matrix* res = new Matrix(mat2);
    int colnum = min(num, ncol);
    int colStart = max(0, ncol - num);
    
    //TODO 

    return res;
}


static Matrix MatrixCalculator::sub(Matrix mat1, Matrix mat2){
    int nrow = mat1.getNRow();
    int ncol = mat1.getNCol();
    Matrix* res = new Matrix(mat1);
    
    cblas_daxpy(nrow * ncol, -1.0, mat2.value, 1, res->value, 1);
    
    return res;
}

static Matrix* MatrixCalculator::sub(Matrix mat1, Matrix mat2, int num){
    int nrow = mat1.getNRow();
    int ncol = mat1.getNCol();
    Matrix* res = new Matrix(mat1);
    
    //TODO    
    return res;
}


static Matrix* MatrixCalculator::div(Matrix mat1, Matrix mat2){
    int nrow = mat1.getNRow();
    int ncol = mat1.getNCol();
    Matrix* res = new Matrix(mat1);
    
    //TODO

    return res;
}

static Matrix* MatrixCalculator::div(Matrix mat1, Matrix mat2, int num){
    int nrow = mat1.getNRow();
    int ncol = mat1.getNCol();
    Matrix* res = new Matrix(mat1);
    
    //TODO

    return res;
}


static Matrix* MatrixCalculator::mul(Matrix mat1, Matrix mat2){
    int nrow = mat2.getNRow();
    int ncol = mat2.getNCol();
    
}

static Matrix* MatrixCalculator::mul(Matrix mat1, Matrix mat2, int num){
    int nrow = mat1.getNRow();
    int ncol = mat1.getNCol();

}

//void cblas_domatcopy (char ordering, char trans, size_t rows, size_t cols, const double alpha, const double * A, size_t lda, double * B, size_t ldb);
static Matrix* MatrixCalculator::mul(double val1, Matrix mat2){
    int nrow = mat2.getNRow();
    int ncol = mat2.getNCol();
    Matrix res = Matrix(nrow, ncol);
    cblas_domatcopy(CblasRowMajor, CblasNoTrans, nrow, ncol, val1, mat2.value, ncol, res->value, ncol);
}

static Matrix* MatrixCalculator::mul(double val1, Matrix mat2, int num){

}
static double* MatrixCalculator::mul(double val1, double* vec2){
    int len = sizeof(vec2) / sizeof(double);
    double* res = new double[len];    
    cblas_daxpy(len, val1, vec2, 1, res, 1);
    return res;
}

//void cblas_dgemm (const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb, const MKL_INT m, 
//                  const MKL_INT n, const MKL_INT k, const double alpha, const double *a, const MKL_INT lda, const double *b, 
//                  const MKL_INT ldb, const double beta, double *c, const MKL_INT ldc);

static Matrix* matrixMul(double* vec1, Matrix mat2){
    int nrow = mat2.getNRow();
    int ncol = mat2.getNCol();
    Matrix* res = new Matrix(1, ncol);
    Matrix mat1 = Matrix(1. nrow, vec1);
    
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 1, nrow, ncol, 1.0, mat1.value, nrow, mat2.value, ncol, 0.0, res->value, ncol);
    return res;
}

static Matrix* matrixMul(Matrix mat1, double* vec2){
    int nrow = mat1.getNRow();
    int ncol = mat1.getNCol();
    Matrix* res = new Matrix(nrow, 1);
    Matrix mat2 = Matrix(ncol, 1, vec2);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nrow, ncol, 1, 1.0, mat1.value, ncol, mat2.value, 1, 0.0, res->value, 1);
    return res;
}


static Matrix* matrixMul(Matrix mat1, Matrix mat2){
    int nrow1 = mat1.getNRow();
    int ncol1 = mat1.getNCol();
    int nrow2 = mat2.getNRow();
    int ncol2 = mat2.getNCol();
    Matrix* res = new Matrix(nrow1, ncol2);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nrow1, ncol1, ncol2, 1.0, mat1.value, ncol1, mat2.value, ncol2, 0.0, res->value, ncol2);
    return res;
}


static Matrix* MatrixCalculator::max(Matrix mat1, Matrix mat2){
    int nrow = mat1.getNRow();
    int ncol = mat1.getNCol();
    Matrix* res = new Matrix(nrow, ncol);
    vdFmax(nrow * ncol, mat1.value, mat2.value, res->value);
    return res;
}

static Matrix* MatrixCalculator::max(Matrix mat1, Matrix mat2, int num){

}

static Matrix* MatrixCalculator::min(Matrix mat1, Matrix mat2){
    int nrow = mat1.getNRow();
    int ncol = mat1.getNCol();
    Matrix* res = new Matrix(nrow, ncol);
    vdFmin(nrow * ncol, mat1.value, mat2.value, res->value);
    return res;
}

static Matrix* MatrixCalculator::min(Matrix mat1, Matrix mat2, int num){

}

static LogicMatrix* MatrixCalculator::bigger(Matrix mat1, Matrix mat2){

}
static LogicMatrix* MatrixCalculator::bigger(Matrix mat1, Matrix mat2, int num){

}
static LogicMatrix* MatrixCalculator::bigger(Matrix mat1, double val){

}
static LogicMatrix* MatrixCalculator::bigger(Matrix mat1, double val, int num){

}

static LogicMatrix* MatrixCalculator::smaller(Matrix mat1, Matrix mat2){

}
static LogicMatrix* MatrixCalculator::smaller(Matrix mat1, Matrix mat2, int num){

}
static LogicMatrix* MatrixCalculator::smaller(Matrix mat1, double val){

}
static LogicMatrix* MatrixCalculator::smaller(Matrix mat1, double val, int num){

}

static LogicMatrix* MatrixCalculator::equal(Matrix mat1, Matrix mat2){

}
static LogicMatrix* MatrixCalculator::equal(Matrix mat1, Matrix mat2, int num){

}
static LogicMatrix* MatrixCalculator::equal(Matrix mat1, double val){

}
static LogicMatrix* MatrixCalculator::equal(Matrix mat1, double val, int num){

}