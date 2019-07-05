#include<iostream>
#include<sys/time.h>
#include<cblas.h>
#include<vector>
#include<unistd.h>
#include<cstdlib>
#include<string>
#include"Matrix.h"
#include"LogicMatrix.h"

class MatrixFactory{
    static Matrix emptyMatrix = new Matrix(0,0);
    static double* emptyVector = new double[0];


    MatrixFactory(){}

    static Matrix getInstanceOfEmptyMatrix(){

    }

    static Matrix getInstanceOfZeroMatrix(int n){

    }

    static Matrix getInstanceOfZeroMatrix(Matrix mat){

    }

    static Matrix getInstanceOfNaNMatrix(int n){

    }

    static Matrix getInstanceOfNaNMatrix(int n, int m){

    }

    static Matrix getExpandingColumnInstanceOfMatrix(Matrix mat, int n){

    }

    static Matrix getInstanceOfDiagMatrix(double* diag){

    }

    static Matrix getInstanceOfRowMatrix(double* vec){

    }

    static Matrix mergeMatrixHorizon(Matrix mat1, Matrix mat2){

    }

    static Matrix mergeMatrixVertical(Matrix mat1, Matrix mat2){

    }

    static Matrix subMatrixHorizonByPeriod(Matrix mat1, int period){

    }

    static Matrix subMatrixHorizonByPeriodAndTruncate(Matrix mat1, int period, int num){

    }

    static LogicMatrix replicateMatrixVertical(bool* logic, int nrow){

    }

    static LogicMatrix replicateMatrixHorizon(bool* logic, int ncol){

    }

    static LogicMatrix subMatrixHorizonByPeriod(LogicMatrix mat1, int period){

    }

    static LogicMatrix subMatrixHorizonByPeriodAndTruncate(LogicMatrix mat1, int period, int num){

    }

    static Matrix subMatrixHorizon(Matrix mat1, int colStart, int colEnd){

    }

    static double* getInstanceOfVectorCopy(double* vec){

    }

    static double* getInstanceOfEmptyVector(){

    }

    static double* getInstanceOfZeroVector(int n){

    }

    static double* getInstanceOfNaNVector(int n){

    }

    static double* getInstanceOfUnitVector(int n){

    }

    static double* mergeVector(double* vec1, double* vec2){

    }

    //static string* mergeStringVector(){}

    static bool* getInstanceOfTrueBoolVector(int n){

    }

    static bool getInstanceOfFalseBoolVector(int n){

    }

    static LogicMatrix getInstanceOfRowLogicMatrix(bool* vec){

    }

    static void copyVectorRight(double* vec, double* host){

    }

    static void copyVectorLeft(double* vec, double* host){
        
    }








}