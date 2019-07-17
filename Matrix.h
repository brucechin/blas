/*
 * @Author: lianke.qin@gmail.com 
 * @Date: 2019-07-08 10:53:02 
 * @Last Modified by: lianke.qin@gmail.com
 * @Last Modified time: 2019-07-09 15:25:58
 */
#ifndef MATRIX_H
#define MATRIX_H

#include<iostream>
//#include<sys/time.h>
#include "include/cblas.h"
#include<vector>
#include<cstdlib>
#include<string>
#include<cstring>
#include<ctime>
#include<memory.h>

class Matrix{

private:
    int nrow;
    int ncol;
    

public:
    double* value;//public attribute for easier access during calculation

    Matrix(){
        nrow = 0;
        ncol = 0;
        value = new double[0];
    }

    Matrix(int n, int m){
        nrow = n;
        ncol = m;
        value = new double[n * m];
        //memset
    }

    //TODO : how to input a 2-D array as parameter
    Matrix(int n, int m, double* vec){
        nrow = n;
        ncol = m;
        value = new double[n * m];
        std::memcpy(value, vec, nrow * ncol * sizeof(double));
    }

    //copy constructor
    Matrix(Matrix& mat){
        nrow = mat.nrow;
        ncol = mat.ncol;
        value = new double[nrow * ncol];
        std::memcpy(value, mat.value, nrow * ncol * sizeof(double));
    }

    //copy some rows or columns usingg memcpy

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

    Matrix* getRowVector(int i){
        Matrix* res = new Matrix(1, ncol);
        for(int j = 0; j < ncol; j++){
            res->value[j] = value[i * ncol + j];
        }
        return res;
    }

    Matrix* getColVector(int j){
        Matrix* res = new Matrix(nrow, 1);
        for(int i = 0; i < nrow; i++){
            res->value[i] = value[i * ncol + j];
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

    void saveMatrix(std::string filename){
        return ;
    }


    void readMatrix(std::string filename){
        return ;
    }

};

#endif