/*
 * @Author: lianke.qin@gmail.com 
 * @Date: 2019-07-08 10:53:17 
 * @Last Modified by: lianke.qin@gmail.com
 * @Last Modified time: 2019-07-09 11:24:15
 */

#ifndef _MATRIXCALCULATOR_H
#define _MATRIXCALCULATOR_H

#include<iostream>
#include <cblas.h>
#include<vector>
#include<cstdlib>
#include"LogicMatrix.h"
#include "MatrixFactory.h"
#include "Matrix.h"

class MatrixCalculator{
private:

    // static double* doubleIntArr;
    // static double* doubleIntDivideArr;

    static double int2Double(int n);
    static double intDoubleDivide(int a, int  b);
    

public:
    MatrixCalculator(){
        // doubleIntArr = new double[1000];
        // for(int i = 0; i < 1000; i++){
        //     doubleIntArr[i] = (double)i;
        // }
        // doubleIntDivideArr = new double[10000];
        // for(int i = 0; i < 100; i++){
        //     for(int j = 0; j < 100; j++){
        //         if(i == 0){
        //             doubleIntDivideArr[i * 100 + j] = 0.0;
        //         }else{
        //             doubleIntDivideArr[i * 100 + j] = (double)j / (double)i;
        //         }
        //     }
        // }
    }
    const static double VALIDITY_PERCENTAGE_REQUIREMENT;
    const static int MAX_HISTORY_LENGTH;
    static Matrix* add(Matrix* mat1, Matrix* mat2);
    static Matrix* add(Matrix* mat1, Matrix* mat2, int num);
    static Matrix* sub(Matrix* mat1, Matrix* mat2);
    static Matrix* sub(Matrix* mat1, Matrix* mat2, int num);
    static Matrix* div(Matrix* mat1, Matrix* mat2);    
    static Matrix* div(Matrix* mat1, Matrix* mat2, int num);   
    static Matrix* mul(Matrix* mat1, Matrix* mat2);    
    static Matrix* mul(Matrix* mat1, Matrix* mat2, int num); 
    static Matrix* mul(double val1, Matrix* mat2);
    static Matrix* mul(double val1, Matrix* mat2, int num);
    static double* mul(double val1, double* vec2);
    static Matrix* matrixMul(double* vec1, Matrix* mat2);
    static Matrix* matrixMul(Matrix* mat1, double* vec2);
    static Matrix* matrixMul(Matrix* mat1, Matrix* mat2);
//TODO : unify all parameter types as Matrix*

    static Matrix* max(Matrix* mat1, Matrix* mat2);
    static Matrix* max(Matrix* mat1, Matrix* mat2, int num);
    static Matrix* min(Matrix* mat1, Matrix* mat2);
    static Matrix* min(Matrix* mat1, Matrix* mat2, int num);

    static LogicMatrix* bigger(Matrix* mat1, Matrix* mat2);
    static LogicMatrix* bigger(Matrix* mat1, Matrix* mat2, int num);
    static LogicMatrix* bigger(Matrix* mat1, double val);
    static LogicMatrix* bigger(Matrix* mat1, double val, int num);

    static LogicMatrix* smaller(Matrix* mat1, Matrix* mat2);
    static LogicMatrix* smaller(Matrix* mat1, Matrix* mat2, int num);
    static LogicMatrix* smaller(Matrix* mat1, double val);
    static LogicMatrix* smaller(Matrix* mat1, double val, int num);

    static LogicMatrix* equal(Matrix* mat1, Matrix* mat2);
    static LogicMatrix* equal(Matrix* mat1, Matrix* mat2, int num);
    static LogicMatrix* equal(Matrix* mat1, double val);
    static LogicMatrix* equal(Matrix* mat1, double val, int num);

    static LogicMatrix* between(Matrix* mat1, double lowerbound, double upperbound);
    static Matrix* betweenValue(Matrix* mat1, double lowerbound, double upperbound);
    static LogicMatrix* matAnd(LogicMatrix* mat1, LogicMatrix* mat2);
    static LogicMatrix* matAnd(LogicMatrix* mat1, LogicMatrix* mat2, int num);
    static LogicMatrix* matOr(LogicMatrix* mat1, LogicMatrix* mat2);
    static LogicMatrix* matOr(LogicMatrix* mat1, LogicMatrix* mat2, int num);
    static LogicMatrix* matNot(LogicMatrix* mat1);
    static LogicMatrix* matNot(LogicMatrix* mat1, int num);

    static Matrix* condition(LogicMatrix* mat1, Matrix* mat2, Matrix* mat3);
    static Matrix* condition(LogicMatrix* mat1, Matrix* mat2, Matrix* mat3, int num);

    static double rankFirst(Matrix* vec, int highIndex, int lowIndex);
    static double* rank(double* vec, int n);
    static Matrix* rank(Matrix* mat);
    static Matrix* rank(Matrix* mat, int num);
    
    static Matrix* round(Matrix* mat);
    static Matrix* round(Matrix* mat, int num);
    static Matrix* floor(Matrix* mat);
    static Matrix* floor(Matrix* mat, int num);
    static Matrix* abs(Matrix* mat);
    static Matrix* abs(Matrix* mat, int num);
    static Matrix* minus(Matrix* mat);
    static Matrix* minus(Matrix* mat, int num);
    static double* minus(double* vec);//DELETE no need for this function
    static Matrix* sqrt(Matrix* mat);
    static Matrix* sqrt(Matrix* mat, int num);
    static Matrix* log(Matrix* mat);
    static Matrix* log(Matrix* mat, int num);
    static Matrix* exp(Matrix* mat);
    static Matrix* exp(Matrix* mat, int num);
    static Matrix* sign(Matrix* mat);
    static Matrix* sign(Matrix* mat, int num);
    static Matrix* inverse(Matrix* mat);
    static Matrix* inverse(Matrix* mat, int num);
    static Matrix* signedpow(Matrix* mat, double index);
    static Matrix* signedpow(Matrix* mat, double index, int num);

    static Matrix* shift(Matrix* mat, int n);
    static Matrix* delay(Matrix* mat, int n);
    static Matrix* delay(Matrix* mat, int n, int num);
    static Matrix* delta(Matrix* mat, int n);
    static Matrix* delta(Matrix* mat, int n, int num);
    static Matrix* ratio(Matrix* mat, int n);
    static Matrix* ratio(Matrix* mat, int n, int num);
    static Matrix* sum(Matrix* mat, int n);
    static Matrix* sum(Matrix* mat, int n, int num);
    static Matrix* product(Matrix* mat, int n);
    static Matrix* product(Matrix* mat, int n, int num);
    
    //NOTE here are some optimized functions with "_op" suffix
    static Matrix* sum_op(Matrix* mat, int n);
    static Matrix* product_op(Matrix* mat, int n);
    static Matrix* tsMax_op(Matrix* mat, int n);
    static Matrix* tsMin_op(Matrix* mat, int n);
    static Matrix* tsArgmax_op(Matrix* mat, int n);
    static Matrix* tsArgmin_op(Matrix* mat, int n);
    static Matrix* tsCountNaN_op(Matrix* mat, int n);
    static Matrix* tsCountTrue_op(LogicMatrix* mat, int n);
    static Matrix* tsCountConsecutiveTrue_op(LogicMatrix* mat, int n);
/*
    TODO optimize these later
    static Matrix* tsMax(Matrix* mat, int n, int num);
    static Matrix* tsMin(Matrix* mat, int n, int num);
    static Matrix* tsArgmax(Matrix* mat, int n, int num);
    static Matrix* tsArgmin(Matrix* mat, int n, int num);
    static Matrix* tsCountNaN(Matrix* mat, int n, int num);
    static Matrix* tsCountTrue(LogicMatrix* mat, int n, int num);
    static Matrix* tsCountConsecutiveTrue(LogicMatrix* mat, int n, int num);
*/


    //DELETE some of these timeseries functions are optimized above
    static Matrix* tsMax(Matrix* mat, int n);
    static Matrix* tsMax(Matrix* mat, int n, int num);
    static Matrix* tsMin(Matrix* mat, int n);
    static Matrix* tsMin(Matrix* mat, int n, int num);
    static Matrix* tsArgmax(Matrix* mat, int n);
    static Matrix* tsArgmax(Matrix* mat, int n, int num);
    static Matrix* tsArgmin(Matrix* mat, int n);
    static Matrix* tsArgmin(Matrix* mat, int n, int num);
    static Matrix* tsRank(Matrix* mat, int n);
    static Matrix* tsRank(Matrix* mat, int n, int num);
    static Matrix* tsMean(Matrix* mat, int n);
    static Matrix* tsMean(Matrix* mat, int n, int num);
    static Matrix* tsStd(Matrix* mat, int n);
    static Matrix* tsStd(Matrix* mat, int n, int num);
    static Matrix* tsSkewness(Matrix* mat, int n);
    static Matrix* tsSkewness(Matrix* mat, int n, int num);
    static Matrix* tsKurtosis(Matrix* mat, int n);
    static Matrix* tsKurtosis(Matrix* mat, int n, int num);
    static Matrix* tsCov(Matrix* mat1, Matrix* mat2, int n);
    static Matrix* tsCov(Matrix* mat1, Matrix* mat2, int n, int num);
    static Matrix* tsCorr(Matrix* mat1, Matrix* mat2, int n);
    static Matrix* tsCorr(Matrix* mat1, Matrix* mat2, int n, int num);
    static Matrix* tsCountNaN(Matrix* mat, int n);
    static Matrix* tsCountNaN(Matrix* mat, int n, int num);
    static Matrix* tsCountTrue(LogicMatrix* mat, int n);
    static Matrix* tsCountTrue(LogicMatrix* mat, int n, int num);
    static Matrix* tsCountConsecutiveTrue(LogicMatrix* mat, int n);
    static Matrix* tsCountConsecutiveTrue(LogicMatrix* mat, int n, int num);
    static Matrix* decayLinear(Matrix* mat, int n);
    static Matrix* decayLinear(Matrix* mat, int n, int num);
    static Matrix* decayExponential(Matrix* mat, int n);
    static Matrix* decayExponential(Matrix* mat, int n, int num);
    static void smoothByDecayLinear(Matrix* mat, int n);
    static void inputNaN(Matrix* mat, double val);

//TODO these are not tested yet
    static void activate(Matrix* mat, double threshold);
    static Matrix* normalize(Matrix* mat, double scale, double mean, double bound);
    static Matrix* normalize(Matrix* mat, double scale, double mean, double bound, int num);
    static void normalizeBySpec(Matrix* mat, double scale, double mean, double bound);





    static Matrix* neutralize(Matrix* mat);
    static Matrix* neutralize(Matrix* mat, int num);
    static Matrix* mean(Matrix* mat);
    static Matrix* unify(Matrix* mat);
    static Matrix* unify(Matrix* mat, int num);
    static Matrix* unifyByL2(Matrix* mat);
    static Matrix* evalValidPct(Matrix* alpha);
    static Matrix* evalAbsSum(Matrix* alpha);
    static Matrix* evalMean(Matrix* alpha);
    static Matrix* evalVariance(Matrix* alpha);
    static Matrix* evalInnerProduction(Matrix* alpha, Matrix* target);
    static Matrix* evalCovariance(Matrix* alpha, Matrix* target);
    static Matrix* evalCorrelation(Matrix* alpha, Matrix* target);


    static double Det(Matrix* mat, int N);
    static double Inverse(Matrix* mat1, int N, Matrix* mat3);
    static Matrix* inv(Matrix* mat);
    static double treat(Matrix* mat);
    static Matrix* diag(Matrix* mat);
    static Matrix* inverseDiag(Matrix* mat);
    static double evalInnerProductionByLongVector(Matrix* alpha, Matrix* target);
    static double evalCorrelationByLongVector(Matrix* alpha, Matrix* target);
    static Matrix* evalBeta(Matrix* alpha, Matrix* target);
    static double evalBetaByLongVector(Matrix* alpha, Matrix* target);
    static Matrix* cumSum(Matrix* ts);
    //static Matrix* cumProd(Matrix* ts);
    static double summaryMean(Matrix* ts, int highIndex, int lowIndex);
    static double summaryVariance(Matrix* ts, int highIndex, int lowIndex);
    static double summarySkewness(Matrix* ts, int highIndex, int lowIndex);
    static double summaryKurtosis(Matrix* ts, int highIndex, int lowIndex);
    static double summaryCovariance(Matrix* ts1, Matrix* ts2, int highIndex, int lowIndex);
    static double summaryCorrelation(Matrix* ts1, Matrix* ts2, int highIndex, int lowIndex);
    static double summarySum(Matrix* ts);
    static double summaryMean(Matrix* ts);
    static double summaryVariance(Matrix* ts);
    static double summarySkewness(Matrix* ts);
    static double summaryKurtosis(Matrix* ts);
    static double summaryCorrelation(Matrix* ts1, Matrix* ts2);
    static double summaryMaxDrawDown(Matrix* ts);
    static double summaryGini(Matrix* ts, int ngroup);


    //TODO:reshape, SVD, QR factorization, transpose, egienvalue, trace, 
};
#endif