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
#include <cblas.h>
#include<vector>
#include<cstdlib>
#include<string>
#include<cstring>
#include<ctime>
#include<memory.h>
#include<fstream>
#include<cmath>
using namespace std;

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
		~Matrix(){
			delete[] value;
	//		std::cout << "deconstruct"<<std::endl;
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
		void setNCol(int m){
			ncol = m;
		}

		void setNRow(int n){
			nrow = n;
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
		bool compareMatrix(Matrix* other) {
			if (nrow != other->getNRow()) return false;
			if (ncol != other->getNCol()) return false;
			double* a = value;
			double* b = other->value;
			int len = nrow * ncol;
			double error_bound = 0.00001;
			for (int i = 0; i < len; i++) {
				double v1 = std::abs(a[i]);
				double v2 = std::abs(b[i]);
				if(v1 > 1){
					if (std::abs(v1 - v2) > error_bound * std::abs(v1 + v2)){
						std::cout<< i / ncol << " "<< i % ncol << " "<< v1 << " " << v2 << std::endl;	
						return false;
					}
				}else{
					if (std::abs(v1 - v2) > error_bound){
						std::cout<< i / ncol << " "<< i % ncol << " "<< v1 << " " << v2 << std::endl;	
						return false;
					}
				}
			}
			return true;
		}

		void print(){
			double* p = value;
			for(int i = 0; i < nrow; i++){
				for(int j = 0; j < ncol; j++){
					std::cout << *p++ << " ";
				}
				std::cout << std::endl;
			}
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
			ofstream file;
			file.open(filename, ios::binary | ios::out);
			if (file.is_open()) {
				file.write((char*)&nrow, sizeof(int));
				file.write((char*)&ncol, sizeof(int));
				file.write((char *)value, sizeof(double) * nrow * ncol);
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
				value = new double[nrow * ncol];
				file.read((char*)value, sizeof(double) * nrow * ncol);
				file.close();
			}
			else {
				cout << "file open failure" << endl;
			}

			return ;
		}

};

#endif
