/*
 * @Author: lianke.qin@gmail.com 
 * @Date: 2019-07-08 10:52:55 
 * @Last Modified by: mikey.zhaopeng
 * @Last Modified time: 2019-07-17 11:31:05
 */

#ifndef LOGICMATRIX_H
#define LOGICMATRIX_H

#include<iostream>
#include <cblas.h>
#include<vector>
#include<cstdlib>
#include<string>
#include<fstream>

//we need to validate the input pamameters before using them

using namespace std;

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

    bool compareMatrix(LogicMatrix* other) {
		if (nrow != other->getNRow()) return false;
		if (ncol != other->getNCol()) return false;
		bool* a = value;
		bool* b = other->value;
		int len = nrow * ncol;
		for (int i = 0; i < len; i++) {
			if ((*a++) != (*b++)) return false;
		}
		return true;
	}

    void print(){
        bool* p = value;
        for(int i = 0; i < nrow; i++){
            for(int j = 0; j < ncol; j++){
                std::cout << *p++ << " ";
            }
            std::cout << std::endl;
        }
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

    void saveMatrix(std::string filename){
		ofstream file;
		file.open(filename, ios::binary | ios::out);
		if (file.is_open()) {
			file.write((char*)&nrow, sizeof(int));
			file.write((char*)&ncol, sizeof(int));
			file.write((char *)value, sizeof(bool) * nrow * ncol);
			file.close();
		}
		else {
			cout << "file open failure" << endl;
		}
		
        return ;
    }


    void readMatrix(std::string filename){
		ifstream file;
		file.open(filename, ios::binary | ios::in);
		if (file.is_open()) {
			file.read((char*)&nrow, sizeof(int));
			file.read((char*)&ncol, sizeof(int));
			delete[] value;
			value = new bool[nrow * ncol];
			file.read((char*)value, sizeof(bool) * nrow * ncol);
			file.close();
		}
		else {
			cout << "file open failure" << endl;
		}
		
        return ;
    }


};

#endif
