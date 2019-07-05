#include<iostream>
#include<sys/time.h>
#include<cblas.h>
#include<vector>
#include<unistd.h>
#include<cstdlib>
#include<string>

class Matrix{

private:
    int nrow;
    int ncol;
    double* value;

public:

    Matrix(){
        this.nrow = 0;
        this.ncol = 0;
        this.value = new double[0][0];
    }

    Matrix(int n, int m){
        this.nrow = n;
        this.ncol = m;
        this.value = new double[n][m];
    }

    //TODO : how to input a 2-D array as parameter
    Matrix(double* mat){
        if(mat)
    }

    //copy constructor
    Matrix(Matrix mat){
        this = mat.getMatrix();
    }

    void clear(){
        this.nrow = 0;
        this.ncol = 0;
        delete[] this.value;
    }

    int getNRow(){
        return this.nrow;
    }

    int getNCol(){
        return this.ncol;
    }

    //return the start pointer of the matrix????
    double* getMatrix(){
        return this.value;
    }

    double getElement(int i, int j){
        return this.value[i][j];
    }

    void setElement(int i, int j, double x){
        this.value[i][j] = x;
    }

    void setValue(double x){
        //can be optimized using OpenBLAS
        for(int i = 0; i < this.nrow; i++){
            for(int j = 0; j < this.ncol; j++){
                this.value[i][j]= x;
            }
        }
    }

    double* getRowVector(int i){
        int ncol = this.getNCol(); 
        double* res = new double[ncol];
        for(int j = 0; j < ncol; j++){
            res[j] = this.getElement(i,j);
        }
        return res;
    }

    double getColVector(int j){
        int nrow = this.getNRow();
        double* res = new double[nrow];
        for(int i = 0; i > nrow; i++){
            res[i] = this.getElement(i,j);
        }        
        return res;
    }


    void print(){
        printf("print called\n");
    }

    void copyTo(Matrix mat){
        return 0;
    }

    void copyToUpperLeft(Matrix mat){
        return 0;
    }

    void copyToLowerRight(Matrix mat){
        return 0;
    }

    void saveMatrix(string filename){
        return 0;
    }


    void readMatrix(string filename){
        return 0;
    }

}