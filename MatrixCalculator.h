#include<iostream>
#include<sys/time.h>
#include<cblas.h>
#include<vector>
#include<unistd.h>
#include<cstdlib>




class MatrixCalculator{
private:

    static double doubleIntAr[];
    static double doubleIntDivideArr[][];

    static double int2Double(int n);
    static double intDoubleDivide(int a, int  b);

public:
    const static int MAXHISTORYLENGTH = 90;

    static Matrix


    //implement + - * /


    //implement bigger, smaller, equal, between(lower, upper)


    //implement and, or, not


    //implement round(round down), floor(round up), abs, sqrt, log, exp, signum, inverse, signedpow

    //iimplement delay, shift, delta, ratio, sum
    



    //time series related
    //tsmax, tsmin, tsargmax, tsargmin, tsrank, tsmean, tsstd, tsskewness, tsskurtosis, tscov, tscorr, tscountnan, tscounttrue, tscountconsecutivetrue

    //delaylinear, decayexponential, smoothByDecayLinear, imputeNaN, 

    //activate, normalize, normalizeBySpec, neutralize, mean, unify, unifyByL2, evalValidPct, evalAbsSum, evalMean, evalVariance, evalInnerProduct, evalCov, evalCorr

    //Det, Inverse, inv, treat, diag, inversediag

    //evalInnerProductionByLongVector, evalCorrByLongVector, evalBeta

    //cumSum, summaryMean, summaryVar, summarySkewness, summaryKurtosis, summaryCov, summaryCorr,summarySum, summaryGini



    //product
}