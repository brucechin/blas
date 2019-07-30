
//#include "gtest/gtest.h"
#include "cblas.h"
#include "../MatrixCalculator.h"
#include <fstream>
#include<iostream>
using namespace std;


//get different results from different Calculator.function and save them into disk for comparation

int rangeMin = -10000;
int rangeMax = 10000;
int matSize = 1000;
int lowerBound = -100;
int upperBound = 100;
int stepSize = 50; // used for shift/delta/delay/ratio/sum/product
MatrixFactory* factory = new MatrixFactory();
MatrixCalculator* calculator = new MatrixCalculator();
string matrix_a = "a.mat";
string matrix_b = "b.mat";
string logicMatrix_a = "la.mat";
string logicMatrix_b = "lb.mat";


int main()
{	
    Matrix* a = factory->getInstanceOfRandomMatrix(matSize, matSize, rangeMin, rangeMax);
	a->saveMatrix(matrix_a);
	Matrix* b = factory->getInstanceOfRandomMatrix(matSize, matSize, rangeMin, rangeMax);
	b->saveMatrix(matrix_b);
	LogicMatrix* la = factory->getInstanceOfRandomLogicMatrix(matSize, matSize);
	la->saveMatrix(logicMatrix_a);
	LogicMatrix* lb = factory->getInstanceOfRandomLogicMatrix(matSize, matSize);
	lb->saveMatrix(logicMatrix_b);
	
    Matrix* res;
    LogicMatrix* lres;

    res = calculator->add(a, b);
    res->saveMatrix("add.mat");
    delete[] res;
    
    res = calculator->sub(a, b);
    res->saveMatrix("sub.mat");
    delete[] res;

    res = calculator->div(a, b);
    res->saveMatrix("div.mat");
    delete[] res;

    res = calculator->mul(a, b);
    res->saveMatrix("mul.mat");
    delete[] res;

    res = calculator->matrixMul(a, b);
    res->saveMatrix("matrixMul.mat");
    delete[] res;

    res = calculator->max(a, b);
    res->saveMatrix("max.mat");
    delete[] res;

    res = calculator->min(a, b);
    res->saveMatrix("min.mat");
    delete[] res;

    lres = calculator->bigger(a, b);
    lres->saveMatrix("bigger.mat");
    delete[] lres;

    lres = calculator->smaller(a, b);
    lres->saveMatrix("smaller.mat");
    delete[] lres;


    lres = calculator->equal(a, b);
    lres->saveMatrix("equal.mat");
    delete[] lres;
	
    lres = calculator->between(a, lowerBound, upperBound);
    lres->saveMatrix("between.mat");
    delete[] lres;

    lres = calculator->matAnd(la, lb);
    lres->saveMatrix("matAnd.mat");
    delete[] lres;

    lres = calculator->matOr(la, lb);
    lres->saveMatrix("matOr.mat");
    delete[] lres;

    lres = calculator->matNot(la);
    lres->saveMatrix("matNot.mat");
    delete[] lres;

    res = calculator->condition(la, a, b);
    res->saveMatrix("condition.mat");
    delete[] res;

    res = calculator->rank(a);
    res->saveMatrix("rank.mat");
    delete[] res;

    res = calculator->round(a);
    res->saveMatrix("round.mat");
    delete[] res;

    res = calculator->floor(a);
    res->saveMatrix("floor.mat");
    delete[] res;

    res = calculator->abs(a);
    res->saveMatrix("abs.mat");
    delete[] res;

    res = calculator->minus(a);
    res->saveMatrix("minus.mat");
    delete[] res;

    res = calculator->sqrt(a);
    res->saveMatrix("sqrt.mat");
    delete[] res;

    res = calculator->log(a);
    res->saveMatrix("log.mat");
    delete[] res;

    res = calculator->exp(a);
    res->saveMatrix("exp.mat");
    delete[] res;

    res = calculator->sign(a);
    res->saveMatrix("sign.mat");
    delete[] res;

    res = calculator->inverse(a);
    res->saveMatrix("inverse.mat");
    delete[] res;

    res = calculator->signedpow(a, 2);
    res->saveMatrix("signedpow.mat");
    delete[] res;

    res = calculator->shift(a, stepSize);
    res->saveMatrix("shift.mat");
    delete[] res;

    res = calculator->delay(a, stepSize);
    res->saveMatrix("delay.mat");
    delete[] res;

    res = calculator->delta(a, stepSize);
    res->saveMatrix("delta.mat");
    delete[] res;

    res = calculator->ratio(a, stepSize);
    res->saveMatrix("ratio.mat");
    delete[] res;

    res = calculator->sum(a, stepSize);
    res->saveMatrix("sum.mat");
    delete[] res;

    res = calculator->product(a, stepSize);
    res->saveMatrix("product.mat");
    delete[] res;


    res = calculator->tsMax(a, stepSize);
    res->saveMatrix("tsmax.mat");
    delete[] res;

    res = calculator->tsMin(a, stepSize);
    res->saveMatrix("tsmin.mat");
    delete[] res;

    res = calculator->tsArgmax(a, stepSize);
    res->saveMatrix("tsargmax.mat");
    delete[] res;

    res = calculator->tsArgmin(a, stepSize);
    res->saveMatrix("tsargmin.mat");
    delete[] res;

    res = calculator->tsRank(a, stepSize);
    res->saveMatrix("tsrank.mat");
    delete[] res;

    res = calculator->tsMean(a, stepSize);
    res->saveMatrix("tsmean.mat");
    delete[] res;

    res = calculator->tsStd(a, stepSize);
    res->saveMatrix("tsstd.mat");
    delete[] res;

    res = calculator->tsSkewness(a, stepSize);
    res->saveMatrix("tsskewness.mat");
    delete[] res;

    res = calculator->tsKurtosis(a, stepSize);
    res->saveMatrix("tskurtosis.mat");
    delete[] res;

    res = calculator->tsCov(a, b, stepSize);
    res->saveMatrix("tscov.mat");
    delete[] res;

    res = calculator->tsCorr(a, b, stepSize);
    res->saveMatrix("tscorr.mat");
    delete[] res;

    res = calculator->tsCountTrue(la, stepSize);
    res->saveMatrix("tscounttrue.mat");
    delete[] res;

    res = calculator->tsCountConsecutiveTrue(la, stepSize);
    res->saveMatrix("tscountconsecutivetrue.mat");
    delete[] res;

    res = calculator->decayLinear(a, stepSize);
    res->saveMatrix("decaylinear.mat");
    delete[] res;

    res = calculator->decayExponential(a, stepSize);
    res->saveMatrix("decayexponential.mat");
    delete[] res;

    res = calculator->neutralize(a);
    res->saveMatrix("neutralize.mat");
    delete[] res;

    res = calculator->unify(a);
    res->saveMatrix("unify.mat");
    delete[] res;

    res = calculator->evalValidPct(a);
    res->saveMatrix("evalvalidpct.mat");
    delete[] res;

    res = calculator->evalAbsSum(a);
    res->saveMatrix("evalabssum.mat");
    delete[] res;

    res = calculator->evalMean(a);
    res->saveMatrix("evalmean.mat");
    delete[] res;

    res = calculator->evalVariance(a);
    res->saveMatrix("evalvariance.mat");
    delete[] res;

    res = calculator->evalInnerProduction(a, b);
    res->saveMatrix("evalinnerproduction.mat");
    delete[] res;

    res = calculator->evalCovariance(a, b);
    res->saveMatrix("evalcovariance.mat");
    delete[] res;

    res = calculator->evalCorrelation(a, b);
    res->saveMatrix("evalcorrelation.mat");
    delete[] res;

    res = calculator->evalBeta(a, b);
    res->saveMatrix("evalbeta.mat");
    delete[] res;

    res = calculator->diag(a);
    res->saveMatrix("diag.mat");
    delete[] res;

    res = calculator->cumSum(a);
    res->saveMatrix("cumsum.mat");
    delete[] res;


    return 0;
}	