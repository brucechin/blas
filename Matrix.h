/*
 * @Author: lianke.qin@gmail.com 
 * @Date: 2019-07-08 10:53:02 
 * @Last Modified by: lianke.qin@gmail.com
 * @Last Modified time: 2019-07-08 16:29:50
 */


#include<iostream>
//#include<sys/time.h>
//#include<cblas.h>
#include<vector>
#include<cstdlib>
#include<string>
#include<ctime>

using namespace std;

class Matrix{

private:
    int nrow;
    int ncol;
    double* value;

public:

    Matrix(){
        nrow = 0;
        ncol = 0;
        value = new double[0];
    }

    Matrix(int n, int m){
        nrow = n;
        ncol = m;
        value = new double[n * m];
    }

    //TODO : how to input a 2-D array as parameter
    Matrix(double* mat){
        
    }

    //copy constructor
    Matrix(Matrix& mat){
        nrow = mat.nrow;
        ncol = mat.ncol;
        value = mat.value;
    }

    void clear(){
        nrow = 0;
        ncol = 0;
        delete[] value;
    }

    int getNRow(){
        return nrow;
    }

    int getNCol(){
        return ncol;
    }

    //return the start pointer of the matrix????
    double* getMatrix(){
        return value;
    }

    double getElement(int i, int j){
        return value[i * ncol + j];
    }

    void setElement(int i, int j, double x){
        value[i * ncol + j] = x;
    }


    //loop unrolling trick is not used here.
    void setValue(double x){
        //can be optimized using OpenBLAS
        for(int i = 0; i < nrow; i++){
            double *cp = &value[i * ncol];
            for(int j = 0; j < ncol; j++){
                *cp++ = x;
            }
        }
    }

    double* getRowVector(int i){
        double* res = new double[ncol];
        for(int j = 0; j < ncol; j++){
            res[j] = getElement(i,j);
        }
        return res;
    }

    double* getColVector(int j){
        double* res = new double[nrow];
        for(int i = 0; i < nrow; i++){
            res[i] = getElement(i,j);
        }        
        return res;
    }


    void print(){
        printf("print called\n");
    }

    //add a row/column copy function
    
    void copyTo(Matrix mat){
        return ;
    }

    void copyToUpperLeft(Matrix mat){
        return ;
    }

    void copyToLowerRight(Matrix mat){
        return ;
    }

    void saveMatrix(string filename){
        return ;
    }


    void readMatrix(string filename){
        return ;
    }

};