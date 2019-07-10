/*
 * @Author: lianke.qin@gmail.com 
 * @Date: 2019-07-08 10:53:17 
 * @Last Modified by: lianke.qin@gmail.com
 * @Last Modified time: 2019-07-09 11:24:15
 */


#include<iostream>
#include"cblas.h"
#include<vector>
//#include<unistd.h>
#include<cstdlib>
#include"LogicMatrix.h"



class MatrixCalculator{
private:

    static double doubleIntArr[];
    static double doubleIntDivideArr[][];

    static double int2Double(int n);
    static double intDoubleDivide(int a, int  b);
    

public:
    static{
        doubleIntArr = new double[1000];
        for(int i = 0; i < 1000; i++){
            doubleIntArr[i] = (double)i;
        }
        doubleIntDivideArr = new double[100][100];
        for(int i = 0; i < 100; i++){
            for(int j = 0; j < 100; j++){
                if(i == 0){
                    doubleIntDivideArr[i][j] = 0.0;
                }else{
                    doubleIntDivideArr[i][j] = (double)j / (double)i;
                }
            }
        }
    }
    const static double VALIDITY_PERCENTAGE_REQUIREMENT = 0.3;
    const static int MAX_HISTORY_LENGTH = 90;
    MatrixCalculator(){}
    static double int2Double(int n);
    static double intDoubleDivide(int a, int b);
    static Matrix* add(Matrix mat1, Matrix mat2);
    static Matrix* add(Matrix mat1, Matrix mat2, int num);
    static Matrix* sub(Matrix mat1, Matrix mat2);
    static Matrix* sub(Matrix mat1, Matrix mat2, int num);
    static Matrix* div(Matrix mat1, Matrix mat2);    
    static Matrix* div(Matrix mat1, Matrix mat2, int num);   
    static Matrix* mul(Matrix mat1, Matrix mat2);    
    static Matrix* mul(Matrix mat1, Matrix mat2, int num); 
    static Matrix* mul(double val1, Matrix mat2);
    static Matrix* mul(double val1, Matrix mat2, int num);
    static double* mul(double val1, double* vec2);
    static Matrix* matrixMul(double* vec1, Matrix mat2);
    static Matrix* matrixMul(Matrix mat1, double* vec2);
    static Matrix* matrixMul(Matrix mat1, Matrix mat2);


    static Matrix* max(Matrix mat1, Matrix mat2);
    static Matrix* max(Matrix mat1, Matrix mat2, int num);
    static Matrix* min(Matrix mat1, Matrix mat2);
    static Matrix* min(Matrix mat1, Matrix mat2, int num);

    static LogicMatrix* bigger(Matrix mat1, Matrix mat2);
    static LogicMatrix* bigger(Matrix mat1, Matrix mat2, int num);
    static LogicMatrix* bigger(Matrix mat1, double val);
    static LogicMatrix* bigger(Matrix mat1, double val, int num);

    static LogicMatrix* smaller(Matrix mat1, Matrix mat2);
    static LogicMatrix* smaller(Matrix mat1, Matrix mat2, int num);
    static LogicMatrix* smaller(Matrix mat1, double val);
    static LogicMatrix* smaller(Matrix mat1, double val, int num);

    static LogicMatrix* equal(Matrix mat1, Matrix mat2);
    static LogicMatrix* equal(Matrix mat1, Matrix mat2, int num);
    static LogicMatrix* equal(Matrix mat1, double val);
    static LogicMatrix* equal(Matrix mat1, double val, int num);

    static LogicMatrix* between(Matrix mat1, double lowerbound, double upperbound);
    static Matrix* betweenValue(Matrix mat1, double lowerbound, double upperbound);
    static LogicMatrix* and(LogicMatrix mat1, LogicMatrix mat2);
    static LogicMatrix* and(LogicMatrix mat1, LogicMatrix mat2, int num);
    static LogicMatrix* or(LogicMatrix mat1, LogicMatrix mat2);
    static LogicMatrix* or(LogicMatrix mat1, LogicMatrix mat2, int num);
    static LogicMatrix* not(LogicMatrix mat1);
    static LogicMatrix* not(LogicMatrix mat1, int num);

    static Matrix* condition(LogicMatrix mat1, Matrix mat2, Matrix mat3);
    static Matrix* condition(LogicMatrix mat1, Matrix mat2, Matrix mat3. int num);

    static double rankFirst(const double* vec, int highIndex, int lowIndex);
    static double* rank(const double* vec);
    static Matrix* rank(Matrix mat);
    static Matrix* rank(Matrix mat, int num);
    
    static Matrix* round(Matrix mat);
    static Matrix* round(Matrix mat, int num);
    static Matrix* floor(Matrix mat);
    static Matrix* floor(Matrix mat, int num);
    static Matrix* abs(Matrix mat);
    static Matrix* abs(Matrix mat, int num);
    static Matrix* minus(Matrix mat);
    static Matrix* minus(Matrix mat, int num);
    static double* minus(double* vec);
    static Matrix* sqrt(Matrix mat);
    static Matrix* sqrt(Matrix mat, int num);
    static Matrix* log(Matrix mat);
    static Matrix* log(Matrix mat, int num);
    static Matrix* exp(Matrix mat);
    static Matrix* exp(Matrix mat, int num);
    static Matrix* sign(Matrix mat);
    static Matrix* sign(Matrix mat, int num);
    static Matrix* inverse(Matrix mat);
    static Matrix* inverse(Matrix mat, int num);
    static Matrix* signedpow(Matrix mat, double index);
    static Matrix* signedpow(Matrix mat, double index, int num);

    static Matrix* shift(Matrix mat, int n);
    static Matrix* delay(Matrix mat, int n);
    static Matrix* delay(Matrix mat, int n, int num);
    static Matrix* delta(Matrix mat, int n);
    static Matrix* delta(Matrix mat, int n, int num);
    static Matrix* ratio(Matrix mat, int n);
    static Matrix* ratio(Matrix mat, int n, int num);
    static Matrix* sum(Matrix mat, int n);
    static Matrix* sum(Matrix mat, int n, int num);
    static Matrix* product(Matrix mat, int n);
    static Matrix* product(Matrix mat, int n, int num);
    

    static Matrix* tsMax(Matrix mat, int n);
    static Matrix* tsMax(Matrix mat, int n, int num);
    static Matrix* tsMin(Matrix mat, int n);
    static Matrix* tsMin(Matrix mat, int n, int num);
    static Matrix* tsArgmax(Matrix mat, int n);
    static Matrix* tsArgmax(Matrix mat, int n, int num);
    static Matrix* tsArgmin(Matrix mat, int n);
    static Matrix* tsArgmin(Matrix mat, int n, int num);
    static Matrix* tsRank(Matrix mat, int n);
    static Matrix* tsRank(Matrix mat, int n, int num);
    static Matrix* tsMean(Matrix mat, int n);
    static Matrix* tsMean(Matrix mat, int n, int num);
    static Matrix* tsStd(Matrix mat, int n);
    static Matrix* tsStd(Matrix mat, int n, int num);
    static Matrix* tsSkewness(Matrix mat, int n);
    static Matrix* tsSkewness(Matrix mat, int n, int num);
    static Matrix* tsKurtosis(Matrix mat, int n);
    static Matrix* tsKurtosis(Matrix mat, int n, int num);
    static Matrix* tsCov(Matrix mat1, Matrix mat2, int n);
    static Matrix* tsCov(Matrix mat1, Matrix mat2, int n, int num);
    static Matrix* tsCorr(Matrix mat1, Matrix mat2, int n);
    static Matrix* tsCorr(Matrix mat1, Matrix mat2, int n, int num);
    static Matrix* tsCountNaN(Matrix mat, int n);
    static Matrix* tsCountNaN(Matrix mat, int n, int num);
    static Matrix* tsCountTrue(LogicMatrix mat, int n);
    static Matrix* tsCountTrue(LogicMatrix mat, int n, int num);
    static Matrix* tsCountConsecutiveTrue(LogicMatrix mat, int n);
    static Matrix* tsCountConsecutiveTrue(LogicMatrix mat, int n, int num);

    static Matrix* decayLinear(Matrix mat, int n);
    static Matrix* decayLinear(Matrix mat, int n, int num);
    static Matrix* decayExponential(Matrix mat, int n);
    static Matrix* decayExponential(Matrix mat, int n, int num);
    static void smoothByDecayLinear(Matrix mat, int n);
    static void inputNaN(Matrix mat, double val);


    static void activate(Matrix mat, double threshold);
    static Matrix* normalize(Matrix mat, double scale, double mean, double bound);
    static Matrix* normalize(Matrix mat, double scale, double mean, double bound, int num);
    static void normalizeBySpec(Matrix mat, double scale, double mean, double bound);
    static Matrix* neutralize(Matrix mat);
    static Matrix* neutralize(Matrix mat, int num);
    static Matrix* mean(Matrix mat);
    static Matrix* unify(Matrix mat);
    static Matrix* unify(Matrix mat, int num);
    static Matrix* unifyByL2(Matrix mat);
    static double* evalValidPct(Matrix alpha);
    static double* evalAbsSum(Matrix alpha);
    static double* evalMean(Matrix alpha);
    static double* evalVariance(Matrix alpha);
    static double* evalInnerProduction(Matrix alpha, Matrix target);
    static double* evalCovariance(Matrix alpha, Matrix target);
    static double* evalCorrelation(Matrix alpha, Matrix target);


    static double Det(double* mat, int N);
    static double Inverse(double* mat1, int N, double* mat3);
    static Matrix* inv(Matrix mat);
    static double treat(Matrix mat);
    static double* diag(Matrix mat);
    static double* inverseDiag(Matrix mat);
    static double evalInnerProductionByLongVector(Matrix alpha, Matrix target);
    static double evalCorrelationByLongVector(Matrix alpha, Matrix target);
    static double* evalBeta(Matrix alpha, Matrix target);
    static double evalBetaByLongVector(Matrix alpha, Matrix target);
    static double* cumSum(double* ts);
    static double summaryMean(double* ts, int highIndex, int lowIndex);
    static double summaryVariance(double* ts, int highIndex, int lowIndex);
    static double summarySkewness(double* ts, int highIndex, int lowIndex);
    static double summaryKurtosis(double* ts, int highIndex, int lowIndex);
    static double summaryCovariance(double* ts, int highIndex, int lowIndex);
    static double summaryCorrelation(double* ts1, double* ts2, int highIndex, int lowIndex);
    static double summarySum(double* ts);
    static double summaryMean(double* ts);
    static double summaryVariance(double* ts);
    static double summarySkewness(double* ts);
    static double summaryKurtosis(double* ts);
    static double summaryCorrelation(double* ts1, double* ts2);
    static double summaryMaxDrawDown(doubel* ts);
    static double summaryGini(double* ts, int ngroup);

}