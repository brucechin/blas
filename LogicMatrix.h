/*
 * @Author: lianke.qin@gmail.com 
 * @Date: 2019-07-08 10:52:55 
 * @Last Modified by: lianke.qin@gmail.com
 * @Last Modified time: 2019-07-09 11:24:28
 */

#ifndef LOGICMATRIX_H
#define LOGICMATRIX_H

#include<iostream>
#include "include/cblas.h"
#include<vector>
//#include<unistd.h>
#include<cstdlib>
#include<string>


//we need to validate the input pamameters before using them


class LogicMatrix{
private:

    int nrow;
    int ncol;
    

public:
    bool* value;
    
    LogicMatrix(){
        nrow = 0;
        ncol = 0;
        value = new bool[0];
    }

    LogicMatrix(int n, int m){
        nrow = n;
        ncol = m;
        value = new bool[n * m];
    }

    LogicMatrix(bool* mat){
        //may not apply here
    }

    LogicMatrix(LogicMatrix& logic){
        //TODO
    }

    int getNRow(){
        return nrow;
    }

    int getNCol(){
        return ncol;
    }

    bool* getLogicMatrix(){
        return value;
    }

    bool getElement(int i, int j){
        return value[i * ncol + j];
    }

    bool getRowVector(int i){
        bool* res = new bool[ncol];
        for(int j = 0; j < ncol; j++){
            res[j] = getElement(i, j);
        }

        return res;
    }

    double* getRowVectorAsDouble(int i){
        double* res = new double[ncol];
        for(int j = 0; j < ncol; j++){
            res[j] = value[i * ncol + j]? 1.0 : 0.0;
        }
        return res;
    }

    void setElement(int i, int j, bool x){

    }

    //set all value to x
    //TODO : optimize it using BLAS?
    void setValue(bool x){
        for(int i = 0; i < nrow; i++){
            for(int j = 0; j < ncol; j++){
                value[i * ncol + j] = x;
            }
        }
    }

    void copyTo(LogicMatrix& mat){
        mat.nrow = nrow;
        mat.ncol = ncol;
        mat.value = new bool[nrow * ncol];
        for(int i = 0; i < nrow; i++){
            for(int j = 0; j < ncol; j++){
                mat.value[i * ncol + j] = value[i * ncol + j];
            }
        }
    }


};

#endif