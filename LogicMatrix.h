/*
 * @Author: lianke.qin@gmail.com 
 * @Date: 2019-07-08 10:52:55 
 * @Last Modified by:   lianke.qin@gmail.com 
 * @Last Modified time: 2019-07-08 10:52:55 
 */


#include<iostream>
#include<sys/time.h>
#include<cblas.h>
#include<vector>
#include<unistd.h>
#include<cstdlib>
#include<string>


//we need to validate the input pamameters before using them


class LogicMatrix{
private:

    int nrow;
    int ncol;
    bool* value;

public:

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

    LogicMatrix(LogicMatrix logic){
        //TODO
    }

    int getNRow(){
        return this.nrow;
    }

    int getNCol(){
        return this.ncol;
    }

    bool* getLogicMatrix(){
        return this.value
    }

    bool getElement(int i, int j){
        return this.value[i][j];
    }

    bool getRowVector(int i){
        bool res = new bool[ncol];
        for(int j = 0; j < ncol; j++){
            res[j] = getElement(i, j);
        }

        return res;
    }

    double getRowVectorAsDouble(int i){
        bool res = new bool[ncol];
        for(int j = 0; j < ncol; j++){
            res[j] = getElement(i, j)? 1.0 : 0.0;
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
                value[i][j] = x;
            }
        }
    }

    void copyTo(LogicMatrix& mat){
        mat.nrow = nrow;
        mat.ncol = ncol;
        mat.value = new bool[nrow * ncol];
        for(int i = 0; i < nrow; i++){
            for(int j = 0; j < ncol; j++){
                mat.value[i][j] = value[i][j];
            }
        }
    }


}