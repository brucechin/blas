/*
 * @Author: lianke.qin@gmail.com 
 * @Date: 2019-07-08 10:52:40 
 * @Last Modified by: lianke.qin@gmail.com
 * @Last Modified time: 2019-07-08 16:56:40
 */

#include<iostream>
#include "include/cblas.h"
#include<vector>
#include<cstdlib>
#include<string>
#include"Matrix.h"
#include"LogicMatrix.h"
#include<cmath>
#include<random>
#include<time.h>

class MatrixFactory{

    MatrixFactory(){}
public:
    static Matrix* getInstanceOfEmptyMatrix(){
        return new Matrix(0, 0);
    }

    static Matrix* getInstanceOfEyeMatrix(int n){
        Matrix* res = new Matrix(n, n);
        for(int i = 0; i < n; i++){
            res->value[i * n + i] = 1.0;
        }
        return res;
    }

    static Matrix* getInstanceOfZeroMatrix(int n){
        Matrix* res = new Matrix(n, n);
        return res;
    }

    static Matrix* getInstanceOfZeroMatrix(Matrix* mat){
        Matrix* res = new Matrix(mat->getNRow(), mat->getNCol());
        return res;
    }

    static Matrix* getInstanceOfNaNMatrix(int n){
        Matrix* res = new Matrix(n, n);
        res->setValue(NAN);
    }

    static Matrix* getInstanceOfNaNMatrix(int n, int m){
        Matrix* res = new Matrix(n, m);
        res->setValue(NAN);
    }

    static Matrix* getExpandingColumnInstanceOfMatrix(Matrix* mat, int n){
        int nrow = mat->getNRow();
        int ncol = mat->getNCol();

        Matrix* res = new Matrix(nrow, ncol + n);
        for(int i = 0; i < nrow; i++){
            for(int j = 0; j < ncol; j++){
                double val = mat->value[i * ncol + j];
                res->value[i * ncol + j] = val;
            }
        }
        return res;
    }

    static Matrix* getInstanceOfDiagMatrix(double* diag){
        int n = sizeof(diag) / sizeof(double);
        Matrix* res = getInstanceOfZeroMatrix(n);

        for(int i = 0; i < n; i++){
            double val = diag[i];
            res->value[i * n + i] = val;
        }
        return res;

    }

    static Matrix* getInstanceOfRowMatrix(double* vec){
        int n = sizeof(vec) / sizeof(double);
        Matrix* res = new Matrix(1, n);
        for(int i = 0; i < n; i++){
            double val = vec[i];
            res->value[0 * n + i] = val;
        }        
        return res;
    }

    //TODO : slow implementation, how to accelarate it using BLAS??? mat.value should be set as public attribute?
    static Matrix* mergeMatrixHorizon(Matrix* mat1, Matrix* mat2){
        int nrow = mat1->getNRow();
        int ncol1 = mat1->getNCol();
        int ncol2 = mat2->getNCol();

        Matrix* res = new Matrix(nrow, ncol1 + ncol2);
        int ncol = ncol1 + ncol2;
        for(int i = 0; i < nrow; i++){
            for(int j = 0; j < ncol1; j++){
                res->value[i * ncol + j] = mat1->value[i * ncol1 + j];  
            }

            for(int j = 0; j < ncol2; j++){
                res->value[i * ncol + j + ncol1] = mat1->value[i * ncol2 + j];  
            }
        }
        return res;
    }

    static Matrix* mergeMatrixVertical(Matrix* mat1, Matrix* mat2){
        int ncol = mat1->getNCol();
        int nrow1 = mat1->getNRow();
        int nrow2 = mat2->getNRow();

        Matrix* res = new Matrix(nrow1 + nrow2, ncol);
        
        //for better cache locality
        for(int i = 0; i < nrow1; i++){
            for(int j = 0; j < ncol; j++){
                res->value[i * ncol + j] = mat1->value[i * ncol + j];
            }
        }

        for(int i = 0; i < nrow2; i++){
            for(int j = 0; j < ncol; j++){
                res->value[(i + nrow1) * ncol + j] = mat1->value[i * ncol + j];
            }
        }
        
    }

    static Matrix* subMatrixHorizonByPeriod(Matrix* mat1, int period){
        int nrow = mat1->getNRow();
        int ncol = mat1->getNCol();

        int colnum = (int) ncol / period;
        Matrix* res = new Matrix(nrow, ncol);
        for(int i = 0; i < nrow; i++){
            int colid = 0;
            for(int j = ncol - 1; j >= 0 && colid < colnum; j -= period){
                int currid = colnum - colid - 1;
                res->value[i * ncol + currid] = mat1->value[i * ncol + j];
                colid++;
            }
        }
        return res;
    }

    static Matrix* subMatrixHorizonByPeriodAndTruncate(Matrix* mat1, int period, int num){
        int nrow = mat1->getNRow();
        int ncol = mat1->getNCol();

        int colnum = min(num, (int) ncol / period);
        Matrix* res = new Matrix(nrow, colnum);

        for(int i = 0; i < nrow; i++){
            int colid = 0;
            for(int j = ncol - 1; j >= 0 && colid < colnum; j -= period){
                int currid = colnum - colid - 1;
                res->value[i * ncol + currid] = mat1->value[i * ncol + j];
                colid++;
            }
        }
        return res;
    }

    static LogicMatrix* replicateMatrixVertical(bool* logic, int nrow){
        int ncol = sizeof(logic) / sizeof(bool);
        LogicMatrix* res = new LogicMatrix(nrow, ncol);
        for(int i = 0; i < nrow; i++){
            for(int j = 0; j < ncol; j++){
                bool val = logic[j];//for better cache locality and less memory access
                res->value[i * ncol + j] = val;
            }
        }
    
        return res;
        
    }

    static LogicMatrix* replicateMatrixHorizon(bool* logic, int ncol){
        int nrow = sizeof(logic) / sizeof(bool);
        LogicMatrix* res = new LogicMatrix(nrow, ncol);

        for(int i = 0; i < nrow; i++){
            bool val = logic[i];
            for(int j = 0; j < ncol; j++){
                res->value[i * ncol + j] = val;
            }
        }

        return res;
    }

    static LogicMatrix* subMatrixHorizonByPeriod(LogicMatrix* mat1, int period){
        int nrow = mat1->getNRow();
        int ncol = mat1->getNCol();
        int colnum = (int) ncol / period;
        LogicMatrix* res = new LogicMatrix(nrow, colnum);

        for(int i = 0; i < nrow; i++){
            int colid = 0;
            for(int j = ncol - 1; j >= 0 && colid < colnum; j -= period){
                int currid = colnum - colid - 1;
                res->value[i * ncol + currid] = mat1->value[i * ncol + j];
                colid++;
            }
        }

        return res;
    }

    static LogicMatrix* subMatrixHorizonByPeriodAndTruncate(LogicMatrix* mat1, int period, int num){
        int nrow = mat1->getNRow();
        int ncol = mat1->getNCol();
        int colnum = min(num, (int) ncol / period);
        LogicMatrix* res = new LogicMatrix(nrow, colnum);

        for(int i = 0; i < nrow; i++){
            int colid = 0;
            for(int j = ncol - 1; j >= 0 && colid < colnum; j -= period){
                int currid = colnum - colid - 1;
                res->value[i * ncol + currid] = mat1->value[i * ncol + j];
                colid++;
            }
        }

        return res;
    }

    static Matrix* subMatrixHorizon(Matrix* mat1, int colStart, int colEnd){
        int nrow = mat1->getNRow();
        int ncol = mat1->getNCol();
        Matrix* res = new Matrix(nrow, colEnd - colStart + 1);
        for(int i = 0; i < nrow; i++){
            for(int j = colStart; j <= colEnd; j++){
                double val = mat1->value[i * ncol + j];
                res->setElement(i, j - colStart, val);
            }
        }
        return res;
    }

    //utilize memcpy function here to copy memory contents in batch???
    static double* getInstanceOfVectorCopy(double* vec){
        int len = sizeof(vec) / sizeof(double);
        double* res = new double[len];
        return res;
    }

    static double* getInstanceOfEmptyVector(){
        return new double[0];
    }

    static double* getInstanceOfZeroVector(int n){
        return new double[n];
    }

    static double* getInstanceOfNaNVector(int n){
        double * res = new double[n];
        for(int i = 0; i < n; i++){
            res[i] = NAN;
        }
        return res;
    }

    static double* getInstanceOfUnitVector(int n){
        double * res = new double[n];
        for(int i = 0; i < n; i++){
            res[i] = 1.0;
        }
        return res;
    }

    double getRandDouble(int min, int max){
        double m1 = (double)(rand()%101)/101.0;
        min++;
        double m2 = (double)(rand()%(max - min + 1) + min);
        m2 -= 1;
        return m1 + m2;
    }

    static Matrix* getInstanceOfRandomMatrix(int n, int m, int min, int max){
        srand((unsigned)time(NULL));
        Matrix* res = new Matrix(n, m);
        double* p = &res->value[0];
        for(int i = 0; i < n; i++){
            for(int j = 0; j < m; j++){
                double m1 = (double)(rand()%101)/101.0;
                min++;
                double m2 = (double)(rand()%(max - min + 1) + min);
                m2 -= 1;
                (*p++) = m1 + m2;
            }
        }
        return res;
    }

    static double* mergeVector(double* vec1, double* vec2){
        int ncol1 = sizeof(vec1) / sizeof(double);
        int ncol2 = sizeof(vec2) / sizeof(double);
        double* res = getInstanceOfZeroVector(ncol1 + ncol2);

        for(int j = 0; j < ncol1; j++){
            //or val = *vec++;
            double val = vec1[j];
            res[j] = val;
        }
        for(int j = 0; j < ncol2; j++){
            double val = vec2[j];
            res[j + ncol1] = val;
        }
        return res;
    }

    //TODO : how to express a vector of string using arrays in C++
    //static string* mergeStringVector(string[] vec1, string[] vec2){}
    //static subStringVector(string[] vec1, int colStart, int colEnd){}

    static bool* getInstanceOfTrueBoolVector(int n){
        bool* res = new bool[n];
        for(int i = 0; i < n; i++){
            res[i] = true;
        }
        return res;
    }

    static bool getInstanceOfFalseBoolVector(int n){
        bool* res = new bool[n];
        for(int i = 0; i < n; i++){
            res[i] = false;
        }
        return res;
    }

    static LogicMatrix* getInstanceOfRowLogicMatrix(bool* vec){
        int n = sizeof(vec) / sizeof(bool);
        LogicMatrix* res = new LogicMatrix(1, n);
        for(int i = 0; i < n; i++){
            bool val = vec[1];
            res->value[0 * n + i] = val;
        }
        return res;
    }

    static void copyVectorRight(double* vec, double* host){
        int id = (sizeof(host) - sizeof(vec)) / sizeof(double);
        int len = sizeof(vec) / sizeof(double);
        for(int i = 0; i < len; i++){
            host[i + id] = vec[i];
        }
    }

    static void copyVectorLeft(double* vec, double* host){
        int len = sizeof(vec) / sizeof(double);
        for(int i = 0; i < len; i++){
            host[i] = vec[i];
        }
    }




};