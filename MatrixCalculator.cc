/*
 * @Author: lianke.qin@gmail.com 
 * @Date: 2019-07-09 10:18:18 
 * @Last Modified by: lianke.qin@gmail.com
 * @Last Modified time: 2019-07-09 15:39:11
 */


#include"MatrixCalculator.h"
#include"Matrix.h"
#include"MatrixFactory.h"
#include<algorithm>
#include<list>

static double MatrixCalculator::int2Double(int n){
    return doubleIntArr[n];
}

static double MatrixCalculator::intDoubleDivide(int a, int b){
    if(a < 100  && b < 100 && a >= 0 && b >= 0){
        return doubleIntDivideArr[b][a];
    }else{
        return int2Double[a] / int2Double[b];
    }
}

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
    int nrow = mat1.getNRow();
    int ncol = mat1.getNCol();
    int colnuum = min(num, ncol);
    int colStart = max(0, ncol - num);
    Matrix* res = new Matrix(nrow, colnum);

    for(int i = 0; i < nrow; i++){
        int colid = 0;
        for(int j =  colStart; j < ncol; j++){
            res->value[ i * ncol + colid] = max(mat1.value[i * ncol + j], mat2.value[i * ncol + j]);
            colid++;
        }
    }
    return res;
}

static Matrix* MatrixCalculator::min(Matrix mat1, Matrix mat2){
    int nrow = mat1.getNRow();
    int ncol = mat1.getNCol();
    Matrix* res = new Matrix(nrow, ncol);
    vdFmin(nrow * ncol, mat1.value, mat2.value, res->value);
    return res;
}

static Matrix* MatrixCalculator::min(Matrix mat1, Matrix mat2, int num){
    int nrow = mat1.getNRow();
    int ncol = mat1.getNCol();
    int colnuum = min(num, ncol);
    int colStart = max(0, ncol - num);
    Matrix* res = new Matrix(nrow, colnum);

    for(int i = 0; i < nrow; i++){
        int colid = 0;
        for(int j =  colStart; j < ncol; j++){
            res->value[ i * ncol + colid] = min(mat1.value[i * ncol + j], mat2.value[i * ncol + j]);
            colid++;
        }
    }
    return res;
}

static LogicMatrix* MatrixCalculator::bigger(Matrix mat1, Matrix mat2){
    int nrow = mat1.getNRow();
    int ncol = mat1.getNCol();
    int len = nrow * ncol;
    LogicMatrix* res = new Matrix(nrow, ncol);
    bool* r = &res->value[0];
    double* v1 = &mat1.value[0];
    double* v2 = &mat2.value[0];
    for(int i = 0; i < len; i++){
        *r++ = (*v1++) > (*v2++) ? true : false;
    }
    return res;
}
static LogicMatrix* MatrixCalculator::bigger(Matrix mat1, Matrix mat2, int num){

}
static LogicMatrix* MatrixCalculator::bigger(Matrix mat1, double val){
    int nrow = mat1.getNRow();
    int ncol = mat1.getNCol();
    int len = nrow * ncol;
    LogicMatrix* res = new Matrix(nrow, ncol);
    bool* r = &res->value[0];
    double* v1 = &mat1.value[0];
    for(int i = 0; i < len; i++){
        *r++ = (*v1++) > val ? true : false;
    }
    return res;
}
static LogicMatrix* MatrixCalculator::bigger(Matrix mat1, double val, int num){

}

static LogicMatrix* MatrixCalculator::smaller(Matrix mat1, Matrix mat2){
    int nrow = mat1.getNRow();
    int ncol = mat1.getNCol();
    int len = nrow * ncol;
    LogicMatrix* res = new Matrix(nrow, ncol);
    bool* r = res->value[0];
    double* v1 = &mat1.value[0];
    double* v2 = &mat2.value[0];
    for(int i = 0; i < len; i++){
        *r++ = (*v1++) < (*v2++) ? true : false;
    }
    return res;
}
static LogicMatrix* MatrixCalculator::smaller(Matrix mat1, Matrix mat2, int num){

}
static LogicMatrix* MatrixCalculator::smaller(Matrix mat1, double val){
    int nrow = mat1.getNRow();
    int ncol = mat1.getNCol();
    int len = nrow * ncol;
    LogicMatrix* res = new Matrix(nrow, ncol);
    bool* r = &res->value[0];
    double* v1 = &mat1.value[0];
    for(int i = 0; i < len; i++){
        *r++ = (*v1++) < val;
    }
    return res;
}
static LogicMatrix* MatrixCalculator::smaller(Matrix mat1, double val, int num){

}

static LogicMatrix* MatrixCalculator::equal(Matrix mat1, Matrix mat2){
    int nrow = mat1.getNRow();
    int ncol = mat1.getNCol();
    int len = nrow * ncol;
    LogicMatrix* res = new Matrix(nrow, ncol);
    bool* r = res->value[0];
    double* v1 = &mat1.value[0];
    double* v2 = &mat2.value[0];
    for(int i = 0; i < len; i++){
        *r++ = (*v1++) == (*v2++);
    }
    return res;
}
static LogicMatrix* MatrixCalculator::equal(Matrix mat1, Matrix mat2, int num){

}
static LogicMatrix* MatrixCalculator::equal(Matrix mat1, double val){
    int nrow = mat1.getNRow();
    int ncol = mat1.getNCol();
    int len = nrow * ncol;
    LogicMatrix* res = new Matrix(nrow, ncol);
    bool* r = &res->value[0];
    double* v1 = &mat1.value[0];
    for(int i = 0; i < len; i++){
        *r++ = (*v1++) == val;
    }
    return res;
}
static LogicMatrix* MatrixCalculator::equal(Matrix mat1, double val, int num){

}

static LogicMatrix* MatrixCalculator::between(Matrix mat1, double lowerbound, double upperbound){
    int nrow = mat1.getNRow();
    int ncol = mat1.getNCol();
    int len = nrow * ncol;
    LogicMatrix* res = new Matrix(nrow, ncol);
    bool* r = &res->value[0];
    double* v1 = &mat1.value[0];
    for(int i = 0; i < len; i++){
        *r++ = ((*v1) > lowerbound) && ((*v1++) <= upperbound);
    }
    return res;
}
static Matrix* MatrixCalculator::betweenValue(Matrix mat1, double lowerbound, double upperbound){

}
static LogicMatrix* MatrixCalculator::and(LogicMatrix mat1, LogicMatrix mat2){
    int nrow = mat1.getNRow();
    int ncol = mat1.getNCol();
    int len = nrow * ncol;
    LogicMatrix* res = new Matrix(nrow, ncol);
    bool* r = &res->value[0];
    bool* v1 = &mat1.value[0];
    bool* v2 = &mat2.value[0];
    for(int i = 0; i < len; i++){
        *r++ = (*v1++) && (*v2++);
    }
    return res;
}
static LogicMatrix* MatrixCalculator::and(LogicMatrix mat1, LogicMatrix mat2, int num){

}
static LogicMatrix* MatrixCalculator::or(LogicMatrix mat1, LogicMatrix mat2){
    int nrow = mat1.getNRow();
    int ncol = mat1.getNCol();
    int len = nrow * ncol;
    LogicMatrix* res = new Matrix(nrow, ncol);
    bool* r = &res->value[0];
    bool* v1 = &mat1.value[0];
    bool* v2 = &mat2.value[0];
    for(int i = 0; i < len; i++){
        *r++ = (*v1++) || (*v2++);
    }
    return res;
}
static LogicMatrix* MatrixCalculator::or(LogicMatrix mat1, LogicMatrix mat2, int num){

}
static LogicMatrix* MatrixCalculator::not(LogicMatrix mat1){
    int nrow = mat1.getNRow();
    int ncol = mat1.getNCol();
    int len = nrow * ncol;
    LogicMatrix* res = new Matrix(nrow, ncol);
    bool* r = &res->value[0];
    bool* v1 = &mat1.value[0];
    for(int i = 0; i < len; i++){
        *r++ = !(*v1++);
    }
    return res;
}
static LogicMatrix* MatrixCalculator::not(LogicMatrix mat1, int num){

}

static Matrix* MatrixCalculator::condition(LogicMatrix mat1, Matrix mat2, Matrix mat3){
    int nrow = mat1.getNRow();
    int ncol = mat1.getNCol();
    int len = nrow * ncol;
    LogicMatrix* res = new Matrix(nrow, ncol);
    bool* r = &res->value[0];
    double* v1 = &mat1.value[0];
    double* v2 = &mat2.value[0];
    double* v3 = &mat3.value[0];
    for(int i = 0; i < len; i++)
        double left = *v2++;
        double right = *v3++;
        *r++ = !(*v1++) ? left : right;
    }
    return res;
}
static Matrix* MatrixCalculator::condition(LogicMatrix mat1, Matrix mat2, Matrix mat3. int num){

}

static double MatrixCalculator::rankFirst(double* vec, int highIndex, int lowIndex){
    int len = highIndex - lowIndex;
    double first = vec[0];
    int biggerCount = 0;
    int count = 0;
    double res = NAN;
    if(!isnan(first)){
        for(int i = highIndex; i >= lowIndex; i--){
            double val = vec[i];
            if(!isnan(val)){
                if(val > first){
                    biggerCount++;
                }
                count++;
            }
        }
        if(intDoubleDivide(count, n) > VALIDITY_PERCENTAGE_REQUIREMENT){
            res = intDoubleDivide(biggerCount, count);
        }

    }

    return res;
}
static double* MatrixCalculator::rank(double* vec){
    int n = sizeof(vec) / sizeof(double);
    double* res = MatrixFactory.getInstanceOfNaNVector(n);
    //TODO
}
static Matrix* MatrixCalculator::rank(Matrix mat){
    int nrow = mat.getNRow();
    int ncol = mat.getNCol();
    Matrix* res = new Matrix(nrow, ncol);

    for(int j = 0 ; j < ncol; j++){
        double* vec = new double[nrow];
        for(int i = 0; i < nrow; i++){
            vec[i] = mat.value[i * ncol + j];
        }
        double* curResult = rank(vec);
        for(int i = 0; i < nrow; i++){
            res[i * col + j] = curResult[i];
        }
    }
    return res;
}
static Matrix* MatrixCalculator::rank(Matrix mat, int num){
    int nrow = mat.getNRow();
    int ncol = mat.getNCol();
    int colnum = min(num, ncol);
    int colStart = max(0, ncol - num);
    Matrix* res = new Matrix(nrow, colnum);
    int colid = 0;
    for(int j = colStart; j < ncol; j++){
        double* vec = new double[nrow];
        for(int i = 0; i < nrow; i++){
            vec[i] = mat.value[i * ncol + j];
        }
        double* curResult = rank(vec);
        for(int i = 0; i < nrow; i++){
            res[i * col + colid] = curResult[i];
        }
        colid++;
    }
    return res;
}

static Matrix MatrixCalculator::round(Matrix mat){
    int nrow = mat.getNRow();
    int ncol = mat.getNCol();
    Matrix* res = new Matrix(nrow, ncol);
    vdRound(nrow * ncol, mat.value, res->value);
    return res;
}
static Matrix MatrixCalculator::round(Matrix mat, int num){

}
static Matrix MatrixCalculator::floor(Matrix mat){
    int nrow = mat.getNRow();
    int ncol = mat.getNCol();
    Matrix* res = new Matrix(nrow, ncol);
    vdFloor(nrow * ncol, mat.value, res->value);
    return res;
}
static Matrix MatrixCalculator::floor(Matrix mat, int num){

}
static Matrix MatrixCalculator::abs(Matrix mat){
    int nrow = mat.getNRow();
    int ncol = mat.getNCol();
    Matrix* res = new Matrix(nrow, ncol);
    vdAbs(nrow * ncol, mat.value, res->value);
    return res;
}
static Matrix MatrixCalculator::abs(Matrix mat, int num){
    
}
static Matrix MatrixCalculator::minus(Matrix mat){
    int nrow = mat.getNRow();
    int ncol = mat.getNCol();
    Matrix* res = new Matrix(nrow, ncol);
    Matrix* zero = MatrixFactory.getInstanceOfZeroMatrix(mat);
    vdSub(nrow * ncol, zero->value, mat.value, res->value);
    return res;
}
static Matrix MatrixCalculator::minus(Matrix mat, int num){

}
static double* MatrixCalculator::minus(double* vec){
    int len = sizeof(vec) / sizeof(double);
    double* res = new double[len];
    double* zero = MatrixFactory.getInstanceOfZeroVector(len);
    vdSub(nrow * ncol, zero, vec, res);
    return res;
}
static Matrix MatrixCalculator::sqrt(Matrix mat){
    int nrow = mat.getNRow();
    int ncol = mat.getNCol();
    Matrix* res = new Matrix(nrow, ncol);
    vdSqrt(nrow * ncol, mat.value, res->value);
    return res;
}
static Matrix MatrixCalculator::sqrt(Matrix mat, int num){

}
static Matrix MatrixCalculator::log(Matrix mat){
    int nrow = mat.getNRow();
    int ncol = mat.getNCol();
    Matrix* res = new Matrix(nrow, ncol);
    vdLn(nrow * ncol, mat.value, res->value);//base is e
    return res;
}
static Matrix MatrixCalculator::log(Matrix mat, int num){

}
static Matrix MatrixCalculator::exp(Matrix mat){
    int nrow = mat.getNRow();
    int ncol = mat.getNCol();
    Matrix* res = new Matrix(nrow, ncol);
    vdExp(nrow * ncol, mat.value, res->value);
    return res;
}
static Matrix MatrixCalculator::exp(Matrix mat, int num){

}
static Matrix MatrixCalculator::sign(Matrix mat){

}
static Matrix MatrixCalculator::sign(Matrix mat, int num){

}
static Matrix MatrixCalculator::inverse(Matrix mat){
    int nrow = mat.getNRow();
    int ncol = mat.getNCol();
    Matrix* res = new Matrix(nrow, ncol);
    vdInv(nrow * ncol, mat.value, res->value);
    return res;
}
static Matrix MatrixCalculator::inverse(Matrix mat, int num){

}
static Matrix MatrixCalculator::signedpow(Matrix mat, double index){

}
static Matrix MatrixCalculator::signedpow(Matrix mat, double index, int num){

}

//move matrix value n units left
static Matrix MatrixCalculator::shift(Matrix mat, int n){
    int nrow = mat.getNRow();
    int ncol = mat.getNCol();
    Matrix* res = new Matrix(nrow, ncol);
    for(int i = 0; i < nrow; i++){
        for(int j = 0; j < ncol; j++){
            if(j >= ncol - n){
                res->value[i * ncol + j] = NAN;
            }else{
                res->value[i * ncol + j] = mat.value[i * ncol + j + n];
            }
        }
    }
    return res;
}
static Matrix MatrixCalculator::delay(Matrix mat, int n){
    int nrow = mat.getNRow();
    int ncol = mat.getNCol();
    Matrix* res = new Matrix(nrow, ncol);
    for(int i = 0; i < nrow; i++){
        for(int j = 0; j < ncol; j++){
            if(j < n){
                res->value[i * ncol + j] = NAN;
            }else{
                res->value[i * ncol + j] = mat.value[i * ncol + j - n];
            }
        }
    }
    return res;
}
static Matrix MatrixCalculator::delay(Matrix mat, int n, int num){

}
static Matrix MatrixCalculator::delta(Matrix mat, int n){
    int nrow = mat.getNRow();
    int ncol = mat.getNCol();
    Matrix* res = new Matrix(nrow, ncol);
    for(int i = 0; i < nrow; i++){
        for(int j = 0; j < ncol; j++){
            if(j < n){
                res->value[i * ncol + j] = NAN;
            }else{
                res->value[i * ncol + j] = mat.value[i * ncol + j] - mat.value[i * ncol + j - n];
            }
        }
    }
    return res;
}
static Matrix MatrixCalculator::delta(Matrix mat, int n, int num){

}
static Matrix MatrixCalculator::ratio(Matrix mat, int n){
    int nrow = mat.getNRow();
    int ncol = mat.getNCol();
    Matrix* res = new Matrix(nrow, ncol);
    for(int i = 0; i < nrow; i++){
        for(int j = 0; j < ncol; j++){
            if(j < n){
                res->value[i * ncol + j] = NAN;
            }else{
                res->value[i * ncol + j] = mat.value[i * ncol + j] / mat.value[i * ncol + j - n];
            }
        }
    }
    return res;
}
static Matrix MatrixCalculator::ratio(Matrix mat, int n, int num){

}
static Matrix MatrixCalculator::sum(Matrix mat, int n){
    int nrow = mat.getNRow();
    int ncol = mat.getNCol();
    Matrix* res = new Matrix(nrow, ncol);

    for(int i = 0; i < nrow; i++){
        for(int j = 0; j < ncol; j++){
            double sum = 0.0;
            int count = 0;
            for(int k = 0; k < n && k <= j; k++){
                double val = mat.value[i * ncol + j - k];
                if(!isnan(val)){
                    sum += val;
                    count++;
                }
            }

            if(intDoubleDivide(count, n) > VALIDITY_PERCENTAGE_REQUIREMENT){
                res->value[i * ccol + j] = sum;
            }else{
                res->value[i * ncol + j] = NAN;
            }
        }
    }

    return res;
}
static Matrix MatrixCalculator::sum(Matrix mat, int n, int num){

}
static Matrix MatrixCalculator::product(Matrix mat, int n){
    int nrow = mat.getNRow();
    int ncol = mat.getNCol();
    Matrix* res = new Matrix(nrow, ncol);

    for(int i = 0; i < nrow; i++){
        for(int j = 0; j < ncol; j++){
            double prod = 1.0;
            int count = 0;
            for(int k = 0; k < n && k <= j; k++){
                double val = mat.value[i * ncol + j - k];
                if(!isnan(val)){
                    prod *= val;
                    count++;
                }
                
            }

            if(intDoubleDivide(count, n) > VALIDITY_PERCENTAGE_REQUIREMENT){
                res->value[i * ncol + j] = prod;
            }else{
                res->value[i * ncol + j] = NAN;
            }
        }
    }
    return res;

}
static Matrix MatrixCalculator::product(Matrix mat, int n, int num){
  
}

static Matrix* MatrixCalculator::tsMax(Matrix mat, int n){
    int nrow = mat.getNRow();
    int ncol = mat.getNCol();
    Matrix* res = new Matrix(nrow, ncol);

    for(int i = 0; i < nrow; i++){
        for(int j = 0; j < ncol; j++){
            double max = NAN;
            bool maxIsNaN = true;
            int count = 0;
            for(int k = 0; k < n && k <= j; k++){
                double val = mat.value[i * ncol + j - k];
                if(!isnan(val)){
                    max = val;
                    maxIsNaN = false;
                }else{
                    max = max(max, val);
                }
                count++;
            }

            if(intDoubleDivide(count, n) > VALIDITY_PERCENTAGE_REQUIREMENT){
                res->value[i * ncol + j] = max;
            }else{
                res->value[i * ncol + j] = NAN;
            }
        }
    }
    return res;
}
static Matrix* MatrixCalculator::tsMax(Matrix mat, int n, int num){

}
static Matrix* MatrixCalculator::tsMin(Matrix mat, int n){
    int nrow = mat.getNRow();
    int ncol = mat.getNCol();
    Matrix* res = new Matrix(nrow, ncol);

    for(int i = 0; i < nrow; i++){
        for(int j = 0; j < ncol; j++){
            double min = NAN;
            bool minIsNaN = true;
            int count = 0;
            for(int k = 0; k < n && k <= j; k++){
                double val = mat.value[i * ncol + j - k];
                if(!isnan(val)){
                    min = val;
                    minIsNaN = false;
                }else{
                    min = min(min, val);
                }
                count++;
            }

            if(intDoubleDivide(count, n) > VALIDITY_PERCENTAGE_REQUIREMENT){
                res->value[i * ncol + j] = min;
            }else{
                res->value[i * ncol + j] = NAN;
            }
        }
    }
    return res;
}
static Matrix* MatrixCalculator::tsMin(Matrix mat, int n, int num){

}
static Matrix* MatrixCalculator::tsArgmax(Matrix mat, int n){
    int nrow = mat.getNRow();
    int ncol = mat.getNCol();
    Matrix* res = new Matrix(nrow, ncol);

    for(int i = 0; i < nrow; i++){
        for(int j = 0; j < ncol; j++){
            double idmax = NAN;
            double max = NAN;
            bool maxIsNaN = true;
            int count = 0;
            for(int k = 0; k < n && k <= j; k++){
                double val = mat.value[i * ncol + j - k];
                if(!isnan(val)){
                    if(maxIsNaN){
                        max = val;
                        maxIsNaN = false;
                        idmax = k;
                    }
                    
                }else{
                    if(val > max){
                        max = val;
                        idmax = k;
                    }
                }
                count++;
            }

            if(intDoubleDivide(count, n) > VALIDITY_PERCENTAGE_REQUIREMENT){
                res->value[i * ncol + j] = idmax;
            }else{
                res->value[i * ncol + j] = NAN;
            }
        }
    }
    return res;
}
static Matrix* MatrixCalculator::tsArgmax(Matrix mat, int n, int num){

}
static Matrix* MatrixCalculator::tsArgmin(Matrix mat, int n){
    int nrow = mat.getNRow();
    int ncol = mat.getNCol();
    Matrix* res = new Matrix(nrow, ncol);

    for(int i = 0; i < nrow; i++){
        for(int j = 0; j < ncol; j++){
            double idmin = NAN;
            double min = NAN;
            bool minIsNaN = true;
            int count = 0;
            for(int k = 0; k < n && k <= j; k++){
                double val = mat.value[i * ncol + j - k];
                if(!isnan(val)){
                    if(minIsNaN){
                        min = val;
                        minIsNaN = false;
                        idmin = k;
                    }
                    
                }else{
                    if(val < min){
                        min = val;
                        idmin = k;
                    }
                }
                count++;
            }

            if(intDoubleDivide(count, n) > VALIDITY_PERCENTAGE_REQUIREMENT){
                res->value[i * ncol + j] = idmin;
            }else{
                res->value[i * ncol + j] = NAN;
            }
        }
    }
    return res;
}
static Matrix* MatrixCalculator::tsArgmin(Matrix mat, int n, int num){

}
static Matrix* MatrixCalculator::tsRank(Matrix mat, int n){
    int nrow = mat.getNRow();
    int ncol = mat.getNCol();
    Matrix* res = new Matrix(nrow, ncol);
    bool* notNanArr = new bool[ncol];
    for(int i = 0; i < nrow; i++){
        Matrix* rowRefa = mat.getRowVector(i);
        for(int j = 0; j < ncol; j++){
            notNanArr[j] = !isnan(mat.value[i * ncol + j]);
        }
        for(int j = 0; j < ncol; j++){
            int count = 0;
            int upperBound = min(j + 1, n);
            int lowestIndex = j + 1 - upperBound;
            for(int k = 0; k < upperBound; k++){
                if(notNanArr[j - k]){
                    count++;
                }
            }
            if(intDoubleDivide(count, n) > VALIDITY_PERCENTAGE_REQUIREMENT){
                double val = rankfirst(rowRefa->value, j, lowestIndex);
                res->value[i * ncol + j] = val;
            }else{
                res->value[i * ncol + j] = NAN;
            }

        }
    }
    return res;
}
static Matrix* MatrixCalculator::tsRank(Matrix mat, int n, int num){

}
static Matrix* MatrixCalculator::tsMean(Matrix mat, int n){
    int nrow = mat.getNRow();
    int ncol = mat.getNCol();
    Matrix* res = new Matrix(nrow, ncol);
    bool* notNanArr = new bool[ncol];
    for(int i = 0; i < nrow; i++){
        Matrix* rowRefa = mat.getRowVector(i);
        for(int j = 0; j < ncol; j++){
            notNanArr[j] = !isnan(mat.value[i * ncol + j]);
        }
        for(int j = 0; j < ncol; j++){
            int count = 0;
            int upperBound = min(j + 1, n);
            int lowestIndex = j + 1 - upperBound;
            for(int k = 0; k < upperBound; k++){
                if(notNanArr[j - k]){
                    count++;
                }
            }
            if(intDoubleDivide(count, n) > VALIDITY_PERCENTAGE_REQUIREMENT){
                double val = summaryMean(rowRefa, j, lowestIndex);
                res->value[i * ncol + j] = val;
            }else{
                res->value[i * ncol + j] = NAN;
            }

        }
    }
    return res;
}
static Matrix* MatrixCalculator::tsMean(Matrix mat, int n, int num){

}
static Matrix* MatrixCalculator::tsStd(Matrix mat, int n){
    int nrow = mat.getNRow();
    int ncol = mat.getNCol();
    Matrix* res = new Matrix(nrow, ncol);
    bool* notNanArr = new bool[ncol];
    for(int i = 0; i < nrow; i++){
        Matrix* rowRefa = mat.getRowVector(i);
        for(int j = 0; j < ncol; j++){
            notNanArr[j] = !isnan(mat.value[i * ncol + j]);
        }
        for(int j = 0; j < ncol; j++){
            int count = 0;
            int upperBound = min(j + 1, n);
            int lowestIndex = j + 1 - upperBound;
            for(int k = 0; k < upperBound; k++){
                if(notNanArr[j - k]){
                    count++;
                }
            }
            if(intDoubleDivide(count, n) > VALIDITY_PERCENTAGE_REQUIREMENT){
                double val = sqrt(summaryVariance(rowRefa, j, lowestIndex));
                res->value[i * ncol + j] = val;
            }else{
                res->value[i * ncol + j] = NAN;
            }

        }
    }
    return res;
}
static Matrix* MatrixCalculator::tsStd(Matrix mat, int n, int num){

}
static Matrix* MatrixCalculator::tsSkewness(Matrix mat, int n){
    int nrow = mat.getNRow();
    int ncol = mat.getNCol();
    Matrix* res = new Matrix(nrow, ncol);
    bool* notNanArr = new bool[ncol];
    for(int i = 0; i < nrow; i++){
        Matrix* rowRefa = mat.getRowVector(i);
        for(int j = 0; j < ncol; j++){
            notNanArr[j] = !isnan(mat.value[i * ncol + j]);
        }
        for(int j = 0; j < ncol; j++){
            int count = 0;
            int upperBound = min(j + 1, n);
            int lowestIndex = j + 1 - upperBound;
            for(int k = 0; k < upperBound; k++){
                if(notNanArr[j - k]){
                    count++;
                }
            }
            if(intDoubleDivide(count, n) > VALIDITY_PERCENTAGE_REQUIREMENT){
                double val = summarySkewness(rowRefa, j, lowestIndex);
                res->value[i * ncol + j] = val;
            }else{
                res->value[i * ncol + j] = NAN;
            }

        }
    }
    return res;
}
static Matrix* MatrixCalculator::tsSkewness(Matrix mat, int n, int num){

}
static Matrix* MatrixCalculator::tsKurtosis(Matrix mat, int n){
    int nrow = mat.getNRow();
    int ncol = mat.getNCol();
    Matrix* res = new Matrix(nrow, ncol);
    bool* notNanArr = new bool[ncol];
    for(int i = 0; i < nrow; i++){
        Matrix* rowRefa = mat.getRowVector(i);
        for(int j = 0; j < ncol; j++){
            notNanArr[j] = !isnan(mat.value[i * ncol + j]);
        }
        for(int j = 0; j < ncol; j++){
            int count = 0;
            int upperBound = min(j + 1, n);
            int lowestIndex = j + 1 - upperBound;
            for(int k = 0; k < upperBound; k++){
                if(notNanArr[j - k]){
                    count++;
                }
            }
            if(intDoubleDivide(count, n) > VALIDITY_PERCENTAGE_REQUIREMENT){
                double val = summaryKurtosis(rowRefa, j, lowestIndex);
                res->value[i * ncol + j] = val;
            }else{
                res->value[i * ncol + j] = NAN;
            }

        }
    }
    return res;
}
static Matrix* MatrixCalculator::tsKurtosis(Matrix mat, int n, int num){

}
static Matrix* MatrixCalculator::tsCov(Matrix mat1, Matrix mat2, int n){
    int nrow = mat.getNRow();
    int ncol = mat.getNCol();
    Matrix* res = new Matrix(nrow, ncol);
    bool* notNanArr = new bool[ncol];
    for(int i = 0; i < nrow; i++){
        Matrix* rowRefa = mat1.getRowVector(i);
        Matrix* rowRefb = mat2.getRowVector(i);
        for(int j = 0; j < ncol; j++){
            notNanArr[j] = !isnan(mat.value[i * ncol + j]);
        }
        for(int j = 0; j < ncol; j++){
            int count = 0;
            int upperBound = min(j + 1, n);
            int lowestIndex = j + 1 - upperBound;
            for(int k = 0; k < upperBound; k++){
                if(notNanArr[j - k]){
                    count++;
                }
            }
            if(intDoubleDivide(count, n) > VALIDITY_PERCENTAGE_REQUIREMENT){
                double val = summaryCovariance(rowRefa, rowRefb, j, lowestIndex);
                res->value[i * ncol + j] = val;
            }else{
                res->value[i * ncol + j] = NAN;
            }

        }
    }
    return res;
}
static Matrix* MatrixCalculator::tsCov(Matrix mat1, Matrix mat2, int n, int num){

}
static Matrix* MatrixCalculator::tsCorr(Matrix mat1, Matrix mat2, int n){
    int nrow = mat.getNRow();
    int ncol = mat.getNCol();
    Matrix* res = new Matrix(nrow, ncol);
    bool* notNanArr = new bool[ncol];
    for(int i = 0; i < nrow; i++){
        Matrix* rowRefa = mat1.getRowVector(i);
        Matrix* rowRefb = mat2.getRowVector(i);
        for(int j = 0; j < ncol; j++){
            notNanArr[j] = !isnan(mat.value[i * ncol + j]);
        }
        for(int j = 0; j < ncol; j++){
            int count = 0;
            int upperBound = min(j + 1, n);
            int lowestIndex = j + 1 - upperBound;
            for(int k = 0; k < upperBound; k++){
                if(notNanArr[j - k]){
                    count++;
                }
            }
            if(intDoubleDivide(count, n) > VALIDITY_PERCENTAGE_REQUIREMENT){
                double val = summaryCorrelation(rowRefa, rowRefb, j, lowestIndex);
                res->value[i * ncol + j] = val;
            }else{
                res->value[i * ncol + j] = NAN;
            }

        }
    }
    return res;
}
static Matrix* MatrixCalculator::tsCorr(Matrix mat1, Matrix mat2, int n, int num){

}
static Matrix* MatrixCalculator::tsCountNaN(Matrix mat, int n){

}
static Matrix* MatrixCalculator::tsCountNaN(Matrix mat, int n, int num){

}
static Matrix* MatrixCalculator::tsCountTrue(LogicMatrix mat, int n){

}
static Matrix* MatrixCalculator::tsCountTrue(LogicMatrix mat, int n, int num){

}
static Matrix* MatrixCalculator::tsCountConsecutiveTrue(LogicMatrix mat, int n){

}
static Matrix* MatrixCalculator::tsCountConsecutiveTrue(LogicMatrix mat, int n, int num){

}

static Matrix* MatrixCalculator::decayLinear(Matrix mat, int n){}
static Matrix* MatrixCalculator::decayLinear(Matrix mat, int n, int num){}
static Matrix* MatrixCalculator::decayExponential(Matrix mat, int n){}
static Matrix* MatrixCalculator::decayExponential(Matrix mat, int n, int num){}
static void MatrixCalculator::smoothByDecayLinear(Matrix mat, int n){}
static void MatrixCalculator::inputNaN(Matrix mat, double val){}


static void MatrixCalculator::activate(Matrix mat, double threshold){}
static Matrix* MatrixCalculator::normalize(Matrix mat, double scale, double mean, double bound){}
static Matrix* MatrixCalculator::normalize(Matrix mat, double scale, double mean, double bound, int num){}
static void MatrixCalculator::normalizeBySpec(Matrix mat, double scale, double mean, double bound){}
static Matrix* MatrixCalculator::neutralize(Matrix mat){}
static Matrix* MatrixCalculator::neutralize(Matrix mat, int num){}
static Matrix* MatrixCalculator::mean(Matrix mat){}
static Matrix* MatrixCalculator::unify(Matrix mat){}
static Matrix* MatrixCalculator::unify(Matrix mat, int num){}
static Matrix* MatrixCalculator::unifyByL2(Matrix mat){}
static double* MatrixCalculator::evalValidPct(Matrix alpha){}
static double* MatrixCalculator::evalAbsSum(Matrix alpha){}
static double* MatrixCalculator::evalMean(Matrix alpha){}
static double* MatrixCalculator::evalVariance(Matrix alpha){}
static double* MatrixCalculator::evalInnerProduction(Matrix alpha, Matrix target){}
static double* MatrixCalculator::evalCovariance(Matrix alpha, Matrix target){}
static double* MatrixCalculator::evalCorrelation(Matrix alpha, Matrix target){}




static double  MatrixCalculator::Det(Matrix mat, int N){}
static double  MatrixCalculator::Inverse(Matrix* mat1, int N, Matrix* mat3){}
static Matrix* MatrixCalculator::inv(Matrix mat){}
static double  MatrixCalculator::treat(Matrix mat){}
static double* MatrixCalculator::diag(Matrix mat){
    int n = mat.getNCol();
    double* res = new double[n];
    for(int i = 0; i < n; i++){
        res[i] = mat.value[i * n + i];
    }
    return res;
}
static double* MatrixCalculator::inverseDiag(Matrix mat){
    int n = mat.getNCol();
    double* res = new double[n];
    for(int i = 0; i < n; i++){
        res[i] = 1.0 / mat.value[i * n + i];
    }
    return res;
}
static double  MatrixCalculator::evalInnerProductionByLongVector(Matrix alpha, Matrix target){}
static double  MatrixCalculator::evalCorrelationByLongVector(Matrix alpha, Matrix target){}
static double* MatrixCalculator::evalBeta(Matrix alpha, Matrix target){}
static double  MatrixCalculator::evalBetaByLongVector(Matrix alpha, Matrix target){}



static double* MatrixCalculator::cumSum(Matrix ts){
    int len = ts.getNCol();
    double* res = MatrixFactory.getInstanceOfZeroVector(len);
    double cumsum = 0;
    for(int i = 0; i < len; i++){
        double val = ts.value[i];
        cumsum += isnan(val) ? 0.0 : val;
        res[i] = cumsum;
    }
    return res;
}
static double  MatrixCalculator::summaryMean(Matrix ts, int highIndex, int lowIndex){
    //ts is 1D matrix
    int len = highIndex - lowIndex + 1;
    double x, sumx;
    sumx = 0;
    int n = 0;
    double mean = NAN;
    for(int k = highIndex; k >= lowIndex; k--){
        x = ts.value[k];
        if(!isnan(x) && !isinf(x)){
            n++;
            sumx += x;
        }
    }
    double f_n = int2Double(n);
    if(f_n > int2Double(len) * VALIDITY_PERCENTAGE_REQUIREMENT){
        double meanx = sumx / f_n;
        mean = meanx;
        mean = (isnan(mean) || isinf(mean)) ? 0 : mean ;
    }
    return mean;
}
static double  MatrixCalculator::summaryVariance(Matrix ts, int highIndex, int lowIndex){
    //ts is 1D matrix
    int len = highIndex - lowIndex + 1;
    double x, sumx, sumxx;
    sumx = 0;
    sumxx = 0;
    int n = 0;
    double var = NAN;
    for(int k = highIndex; k >= lowIndex; k--){
        x = ts.value[k];
        if(!isnan(x) && !isinf(x)){
            n++;
            sumx += x;
            sumxx += x * x;
        }
    }
    double f_n = int2Double(n);
    if(f_n > int2Double(len) * VALIDITY_PERCENTAGE_REQUIREMENT){
        double varx = (sumxx - sumx * sumx /f_n) / f_n;
        var = varx;
        var = (isnan(var) || isinf(var)) ? 0 : var;
    }
    return var;

}
static double  MatrixCalculator::summarySkewness(Matrix ts, int highIndex, int lowIndex){
    //ts is 1D matrix
    int len = highIndex - lowIndex + 1;
    double x, sumx, sumxx, sumxxx;
    sumx = 0;
    sumxx = 0;
    sumxxx = 0;
    int n = 0;
    double skewness = NAN;
    for(int k = highIndex; k >= lowIndex; k--){
        x = ts.value[k];
        if(!isnan(x) && !isinf(x)){
            n++;
            sumx += x;
            sumxx += x * x;
            sumxxx += x * x * x;
        }
    }
    double f_n = int2Double(n);
    if(f_n > int2Double(len) * VALIDITY_PERCENTAGE_REQUIREMENT){
        double u = sumx / f_n;
        double sigma - sqrt(sumxx / f_n - u * u);
        double ex3 = sumxxx / f_n;
        skewness = (ex3 = 3 * u * sigma * sigma - u * u * u) / (sigma * sigma * sigma);
        skewness = (isnan(skewness) || isinf(skewness)) ? 0 : skewness;
    }
    return skewness;
}
static double  MatrixCalculator::summaryKurtosis(Matrix ts, int highIndex, int lowIndex){
    //ts is 1D matrix
    int len = highIndex - lowIndex + 1;
    double x, sumx, sumxx, sumxxx;
    sumx = 0;
    sumxx = 0;
    sumxxx = 0;
    sumxxxx = 0;
    int n = 0;
    double kurtosis = NAN;
    for(int k = highIndex; k >= lowIndex; k--){
        x = ts.value[k];
        if(!isnan(x) && !isinf(x)){
            n++;
            sumx += x;
            sumxx += x * x;
            sumxxx += x * x * x;
            sumxxxx += x * x * x * x;
        }
    }
    double f_n = int2Double(n);
    if(f_n > int2Double(len) * VALIDITY_PERCENTAGE_REQUIREMENT){
        double ex1 = sumx / f_n;
        double ex2 = sumxx / f_n;
        double ex3 = sumxxx / f_n;
        double ex4 = sumxxxx / f_n;
        double sigm2 = ex2 - ex1 * ex1;
        kurtosis = (ex4 - 4 * ex1 * ex3 + 6 * ex1 * ex1 * ex2 - 3 * ex1 * ex1 *ex1 * ex1) / (sigma2 * sigma2) - 3.0;
        kurtosis = (isnan(kurtosis) || isinf(kurtosis)) ? 0 : kurtosis;
    }
    return kurtosis;
}
static double  MatrixCalculator::summaryCovariance(Matrix ts1, Matrix ts2, int highIndex, int lowIndex){
    //ts is 1D matrix
    int len = highIndex - lowIndex + 1;
    double x, sumx, sumy, sumxy;
    sumx = 0;
    sumy = 0;
    sumxy = 0;
    int n = 0;
    double cov = NAN;
    for(int k = 0; k < len; k++){
        x = ts1.value[k];
        y = ts2.value[k];
        if(!isnan(x) && !isinf(x)){
            y = (isnan(y) || isinf(y)) ? 0 : y;
            n++;
            sumx += x;
            sumy += y;
            sumxy += x * y;
        }
    }
    double f_n = int2Double(n);
    if(n > (double)len * VALIDITY_PERCENTAGE_REQUIREMENT){
        cov = (sumxy - sumx * sumy / f_n) / f_n;
        cov = (isnan(cov) || isinf(cov)) ? 0 : cov;
    }
    return cov;

}
static double  MatrixCalculator::summaryCorrelation(Matrix ts1, Matrix ts2, int highIndex, int lowIndex){
    double x, y, sumx, sumy, sumxx, sumyy, sumxy, n;
    sumx = 0;
    sumxx = 0;
    sumy = 0;
    sumyy = 0;
    sumxy = 0;
    n = 0;
    double corr = NAN;
    int len = highIndex - lowIndex + 1;
    for(int k = 0; k < len; k++){
        x = ts1.value[k];
        y = ts2.value[k];
        if(!isnan(x) && !isinf(x)){
            y = (isnan(y) || isinf(y)) ? 0 : y;
            n++;
            sumx += x;
            sumx += x * x;
            sumy += y;
            sumyy += y * y;
            sumxy += x * y;
        }
    }

    if(n > (double)len * VALIDITY_PERCENTAGE_REQUIREMENT){
        double varx = (sumxx - sumx * sumx / n) / n;
        double vary = (sumyy - sumy * sumy / n) / n;
        double cov = (sumxy - sumx * sumy / n) / n;
        corr = cov / sqrt(varx * vary);
        corr = isnan(corr) ? 0 : corr;
        corr = (corr > 1) ? 1 : corr;
        corr = (corr < -1) ? -1 : corr;
    }
    return corr;
}
static double  MatrixCalculator::summarySum(Matrix ts){
    int len = ts.getNCol();
    double x, sumx, n;
    sumx = 0;
    n = 0;
    double sum = NAN;
    for(int k = 0; k < len; k++){
        x = ts.value[k];
        if(!isnan(x) && !isinf(x)){
            n++;
            sumx += x;
        }
    }

    if(n > (double)len * VALIDITY_PERCENTAGE_REQUIREMENT){
        sum = isnan(sumx) ? 0 : sumx;
        sum = isinf(sum) ? 0 : sum;
    }
    return sum;
}
static double  MatrixCalculator::summaryMean(Matrix ts){
    int len = ts.getNCol();
    double x, sumx;
    sumx = 0;
    int n = 0;
    double mean = NAN;
    for(int k = highIndex; k >= lowIndex; k--){
        x = ts.value[k];
        if(!isnan(x) && !isinf(x)){
            n++;
            sumx += x;
        }
    }
    double f_n = int2Double(n);
    if(f_n > int2Double(len) * VALIDITY_PERCENTAGE_REQUIREMENT){
        double meanx = sumx / f_n;
        mean = meanx;
        mean = (isnan(mean) || isinf(mean)) ? 0 : mean ;
    }
    return mean;
}
static double  MatrixCalculator::summaryVariance(Matrix ts){
    int len = ts.getNCol();
    double x, sumx, sumxx;
    sumx = 0;
    sumxx = 0;
    int n = 0;
    double var = NAN;
    for(int k = highIndex; k >= lowIndex; k--){
        x = ts.value[k];
        if(!isnan(x) && !isinf(x)){
            n++;
            sumx += x;
            sumxx += x * x;
        }
    }
    double f_n = int2Double(n);
    if(f_n > int2Double(len) * VALIDITY_PERCENTAGE_REQUIREMENT){
        double varx = (sumxx - sumx * sumx /f_n) / f_n;
        var = varx;
        var = (isnan(var) || isinf(var)) ? 0 : var;
    }
    return var;
}
static double  MatrixCalculator::summarySkewness(Matrix ts){
    int len = ts.getNCol();
    double x, sumx, sumxx, sumxxx;
    sumx = 0;
    sumxx = 0;
    sumxxx = 0;
    int n = 0;
    double skewness = NAN;
    for(int k = highIndex; k >= lowIndex; k--){
        x = ts.value[k];
        if(!isnan(x) && !isinf(x)){
            n++;
            sumx += x;
            sumxx += x * x;
            sumxxx += x * x * x;
        }
    }
    double f_n = int2Double(n);
    if(f_n > int2Double(len) * VALIDITY_PERCENTAGE_REQUIREMENT){
        double u = sumx / f_n;
        double sigma - sqrt(sumxx / f_n - u * u);
        double ex3 = sumxxx / f_n;
        skewness = (ex3 = 3 * u * sigma * sigma - u * u * u) / (sigma * sigma * sigma);
        skewness = (isnan(skewness) || isinf(skewness)) ? 0 : skewness;
    }
    return skewness;
}
static double  MatrixCalculator::summaryKurtosis(Matrix ts){
    int len = ts.getNCol();
    double x, sumx, sumxx, sumxxx;
    sumx = 0;
    sumxx = 0;
    sumxxx = 0;
    sumxxxx = 0;
    int n = 0;
    double kurtosis = NAN;
    for(int k = highIndex; k >= lowIndex; k--){
        x = ts.value[k];
        if(!isnan(x) && !isinf(x)){
            n++;
            sumx += x;
            sumxx += x * x;
            sumxxx += x * x * x;
            sumxxxx += x * x * x * x;
        }
    }
    double f_n = int2Double(n);
    if(f_n > int2Double(len) * VALIDITY_PERCENTAGE_REQUIREMENT){
        double ex1 = sumx / f_n;
        double ex2 = sumxx / f_n;
        double ex3 = sumxxx / f_n;
        double ex4 = sumxxxx / f_n;
        double sigm2 = ex2 - ex1 * ex1;
        kurtosis = (ex4 - 4 * ex1 * ex3 + 6 * ex1 * ex1 * ex2 - 3 * ex1 * ex1 *ex1 * ex1) / (sigma2 * sigma2) - 3.0;
        kurtosis = (isnan(kurtosis) || isinf(kurtosis)) ? 0 : kurtosis;
    }
    return kurtosis;
}
static double  MatrixCalculator::summaryCorrelation(Matrix ts1, Matrix ts2){
    double x, y, sumx, sumy, sumxx, sumyy, sumxy, n;
    sumx = 0;
    sumxx = 0;
    sumy = 0;
    sumyy = 0;
    sumxy = 0;
    n = 0;
    double corr = NAN;
    int len = ts1.getNCol();
    for(int k = 0; k < len; k++){
        x = ts1.value[k];
        y = ts2.value[k];
        if(!isnan(x) && !isinf(x)){
            y = (isnan(y) || isinf(y)) ? 0 : y;
            n++;
            sumx += x;
            sumx += x * x;
            sumy += y;
            sumyy += y * y;
            sumxy += x * y;
        }
    }

    if(n > (double)len * VALIDITY_PERCENTAGE_REQUIREMENT){
        double varx = (sumxx - sumx * sumx / n) / n;
        double vary = (sumyy - sumy * sumy / n) / n;
        double cov = (sumxy - sumx * sumy / n) / n;
        corr = cov / sqrt(varx * vary);
        corr = isnan(corr) ? 0 : corr;
        corr = (corr > 1) ? 1 : corr;
        corr = (corr < -1) ? -1 : corr;
    }
    return corr;
}
static double  MatrixCalculator::summaryMaxDrawDown(Matrix ts){
    int len = ts.getNCol();
    double x, peak, drawdown, n;
    peak = 0;
    drawdown = 0;
    n = 0;
    double drawdownRatio = NAN;
    for(int k = 0; k < len; k++){
        x = ts.value[k];
        if(!isnan(x) && !isinf(x)){
            n++;
            peak = max(x, peak);
            drawdown = max(peak - x, drawdown);
        }
    }

    if(n > (double)len * VALIDITY_PERCENTAGE_REQUIREMENT){
        drawdownRatio = drawdown / max(peak, 0.0);
    }
    return drawdownRatio;
}
static double  MatrixCalculator::summaryGini(Matrix ts, int ngroup){
    int len = ts.getNCol();
    double* sumArray = new double[ngroup];
    for(int i = 0; i < len; i++){
        double val = ts.value[i];
        int idgroup = (i * ngroup) / len;
        sumArray[idgroup] += val;
    }
    list<double> sumList = new list<double>();
    for(int i = 0; i < ngroup; i++){
        sumList.
    }
    //TODO

}