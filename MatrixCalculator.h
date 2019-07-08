/*
 * @Author: lianke.qin@gmail.com 
 * @Date: 2019-07-08 10:53:17 
 * @Last Modified by: lianke.qin@gmail.com
 * @Last Modified time: 2019-07-08 20:55:04
 */


#include<iostream>
#include<sys/time.h>
#include<cblas.h>
#include<vector>
#include<unistd.h>
#include<cstdlib>
#include"LogicMatrix.h"



class MatrixCalculator{
private:

    static double doubleIntAr[];
    static double doubleIntDivideArr[][];

    static double int2Double(int n);
    static double intDoubleDivide(int a, int  b);

public:
    const static int MAXHISTORYLENGTH = 90;

    static Matrix add(Matrix mat1, Matrix mat2);
    static Matrix add(Matrix mat1, Matrix mat2, int num);
    static Matrix sub(Matrix mat1, Matrix mat2);
    static Matrix sub(Matrix mat1, Matrix mat2, int num);
    static Matrix div(Matrix mat1, Matrix mat2);    
    static Matrix div(Matrix mat1, Matrix mat2, int num);   
    static Matrix mul(Matrix mat1, Matrix mat2);    
    static Matrix mul(Matrix mat1, Matrix mat2, int num); 
    static Matrix mul(double val1, Matrix mat2);
    static Matrix mul(double val1, Matrix mat2, int num);
    static double* mul(double val1, double* vec2);
    static Matrix matrixMul(double* vec1, Matrix mat2);
    static Matrix matrixMul(Matrix mat1, double* vec2);
    static Matrix matrixMul(Matrix mat1, Matrix mat2);




    static max(Matrix mat1, Matrix mat2);
    static max(Matrix mat1, Matrix mat2, int num);
    static min(Matrix mat1, Matrix mat2);
    static min(Matrix mat1, Matrix mat2, int num);

    static LogicMatrix bigger(Matrix mat1, Matrix mat2);
    static LogicMatrix bigger(Matrix mat1, Matrix mat2, int num);
    static LogicMatrix bigger(Matrix mat1, double val);
    static LogicMatrix bigger(Matrix mat1, double val, int num);

    static LogicMatrix smaller(Matrix mat1, Matrix mat2);
    static LogicMatrix smaller(Matrix mat1, Matrix mat2, int num);
    static LogicMatrix smaller(Matrix mat1, double val);
    static LogicMatrix smaller(Matrix mat1, double val, int num);

    static LogicMatrix equal(Matrix mat1, Matrix mat2);
    static LogicMatrix equal(Matrix mat1, Matrix mat2, int num);
    static LogicMatrix equal(Matrix mat1, double val);
    static LogicMatrix equal(Matrix mat1, double val, int num);

    static LogicMatrix between(Matrix mat1, double lowerbound, double upperbound);
    static Matrix betweenValue(Matrix mat1, double lowerbound, double upperbound);
    static LogicMatrix and(LogicMatrix mat1, LogicMatrix mat2);
    static LogicMatrix and(LogicMatrix mat1, LogicMatrix mat2, int num);
    static LogicMatrix or(LogicMatrix mat1, LogicMatrix mat2);
    static LogicMatrix or(LogicMatrix mat1, LogicMatrix mat2, int num);
    static LogicMatrix not(LogicMatrix mat1);
    static LogicMatrix not(LogicMatrix mat1, int num);

    static Matrix condition(LogicMatrix mat1, Matrix mat2, Matrix mat3);
    static Matrix condition(LogicMatrix mat1, Matrix mat2, Matrix mat3. int num);

    static double rankFirst(const double* vec, int highIndex, int lowIndex);
    static double* rank(const double* vec);
    static Matrix rank(Matrix mat);
    static Matrix rank(Matrix mat, int num);
    
    static Matrix round(Matrix mat);
    static Matrix round(Matrix mat, int num);
    static Matrix floor(Matrix mat);
    static Matrix floor(Matrix mat, int num);
    static Matrix abs(Matrix mat);
    static Matrix abs(Matrix mat, int num);
    static Matrix minus(Matrix mat);
    static Matrix minus(Matrix mat, int num);
    static double* minus(double* vec);
    static Matrix sqrt(Matrix mat);
    static Matrix sqrt(Matrix mat, int num);
    static Matrix log(Matrix mat);
    static Matrix log(Matrix mat, int num);
    static Matrix exp(Matrix mat);
    static Matrix exp(Matrix mat, int num);
    static Matrix sign(Matrix mat);
    static Matrix sign(Matrix mat, int num);
    static Matrix inverse(Matrix mat);
    static Matrix inverse(Matrix mat, int num);
    static Matric signedpow(Matrix mat, double index);
    static Matrix signedpow(Matrix mat, double index, int num);

    static Matrix shift(Matrix mat, int n);
    static Matrix delay(Matrix mat, int n);
    static Matrix delay(Matrix mat, int n, int num);
    static Matrix delta(Matrix mat, int n);
    static Matrix delta(Matrix mat, int n, int num);
    static Matrix ratio(Matrix mat, int n);
    static Matrix ratio(Matrix mat, int n, int num);
    static Matrix sum(Matrix mat, int n);
    static Matrix sum(Matrix mat, int n, int num);
    static Matrix product(Matrix mat, int n);
    static Matrix product(Matrix mat, int n, int num);
    


    

    //time series related
    //tsmax, tsmin, tsargmax, tsargmin, tsrank, tsmean, tsstd, tsskewness, tsskurtosis, tscov, tscorr, tscountnan, tscounttrue, tscountconsecutivetrue

    //delaylinear, decayexponential, smoothByDecayLinear, imputeNaN, 

    //activate, normalize, normalizeBySpec, neutralize, mean, unify, unifyByL2, evalValidPct, evalAbsSum, evalMean, evalVariance, evalInnerProduct, evalCov, evalCorr

    //Det, Inverse, inv, treat, diag, inversediag

    //evalInnerProductionByLongVector, evalCorrByLongVector, evalBeta

    //cumSum, summaryMean, summaryVar, summarySkewness, summaryKurtosis, summaryCov, summaryCorr,summarySum, summaryGini



    //product
}