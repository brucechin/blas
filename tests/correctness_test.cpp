
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
    res->clear();
    
    res = calculator->sub(a, b);
    res->saveMatrix("sub.mat");
    res->clear();

    res = calculator->div(a, b);
    res->saveMatrix("div.mat");
    res->clear();

    res = calculator->mul(a, b);
    res->saveMatrix("mul.mat");
    res->clear();

    res = calculator->matrixMul(a, b);
    res->saveMatrix("matrixMul.mat");
    res->clear();

    res = calculator->max(a, b);
    res->saveMatrix("max.mat");
    res->clear();

    res = calculator->min(a, b);
    res->saveMatrix("min.mat");
    res->clear();

    lres = calculator->bigger(a, b);
    lres->saveMatrix("bigger.mat");
    lres->clear();

    lres = calculator->smaller(a, b);
    lres->saveMatrix("smaller.mat");
    lres->clear();


    lres = calculator->equal(a, b);
    lres->saveMatrix("equal.mat");
    lres->clear();
	
    lres = calculator->between(a, lowerBound, upperBound);
    lres->saveMatrix("between.mat");
    lres->clear();

    lres = calculator->matAnd(la, lb);
    lres->saveMatrix("matAnd.mat");
    lres->clear();

    lres = calculator->matOr(la, lb);
    lres->saveMatrix("matOr.mat");
    lres->clear();

    lres = calculator->matNot(la);
    lres->saveMatrix("matNot.mat");
    lres->clear();

    res = calculator->condition(la, a, b);
    res->saveMatrix("condition.mat");
    res->clear();

    res = calculator->rank(a);
    res->saveMatrix("rank.mat");
    res->clear();

    res = calculator->round(a);
    res->saveMatrix("round.mat");
    res->clear();

    res = calculator->floor(a);
    res->saveMatrix("floor.mat");
    res->clear();

    res = calculator->abs(a);
    res->saveMatrix("abs.mat");
    res->clear();

    res = calculator->minus(a);
    res->saveMatrix("minus.mat");
    res->clear();

    res = calculator->sqrt(a);
    res->saveMatrix("sqrt.mat");
    res->clear();

    res = calculator->log(a);
    res->saveMatrix("log.mat");
    res->clear();

    res = calculator->exp(a);
    res->saveMatrix("exp.mat");
    res->clear();

    res = calculator->sign(a);
    res->saveMatrix("sign.mat");
    res->clear();

    res = calculator->inverse(a);
    res->saveMatrix("inverse.mat");
    res->clear();

    res = calculator->signedpow(a, 2);
    res->saveMatrix("signedpow.mat");
    res->clear();

    res = calculator->shift(a, stepSize);
    res->saveMatrix("shift.mat");
    res->clear();

    res = calculator->delay(a, stepSize);
    res->saveMatrix("delay.mat");
    res->clear();

    res = calculator->delta(a, stepSize);
    res->saveMatrix("delta.mat");
    res->clear();

    res = calculator->ratio(a, stepSize);
    res->saveMatrix("ratio.mat");
    res->clear();

    res = calculator->sum(a, stepSize);
    res->saveMatrix("sum.mat");
    res->clear();

    res = calculator->product(a, stepSize);
    res->saveMatrix("product.mat");
    res->clear();


    res = calculator->tsMax(a, stepSize);
    res->saveMatrix("tsmax.mat");
    res->clear();

    res = calculator->tsMin(a, stepSize);
    res->saveMatrix("tsmin.mat");
    res->clear();

    res = calculator->tsArgmax(a, stepSize);
    res->saveMatrix("tsargmax.mat");
    res->clear();

    res = calculator->tsArgmin(a, stepSize);
    res->saveMatrix("tsargmin.mat");
    res->clear();

    res = calculator->tsRank(a, stepSize);
    res->saveMatrix("tsrank.mat");
    res->clear();

    res = calculator->tsMean(a, stepSize);
    res->saveMatrix("tsmean.mat");
    res->clear();

    res = calculator->tsStd(a, stepSize);
    res->saveMatrix("tsstd.mat");
    res->clear();

    res = calculator->tsSkewness(a, stepSize);
    res->saveMatrix("tsskewness.mat");
    res->clear();

    res = calculator->tsKurtosis(a, stepSize);
    res->saveMatrix("tskurtosis.mat");
    res->clear();

    res = calculator->tsCov(a, b, stepSize);
    res->saveMatrix("tscov.mat");
    res->clear();

    res = calculator->tsCorr(a, b, stepSize);
    res->saveMatrix("tscorr.mat");
    res->clear();

    res = calculator->tsCountTrue(la, stepSize);
    res->saveMatrix("tscounttrue.mat");
    res->clear();

    res = calculator->tsCountConsecutiveTrue(la, stepSize);
    res->saveMatrix("tscountconsecutivetrue.mat");
    res->clear();

    res = calculator->decayLinear(a, stepSize);
    res->saveMatrix("decaylinear.mat");
    res->clear();

    res = calculator->decayExponential(a, stepSize);
    res->saveMatrix("decayexponential.mat");
    res->clear();

    res = calculator->neutralize(a);
    res->saveMatrix("neutralize.mat");
    res->clear();

    res = calculator->unify(a);
    res->saveMatrix("unify.mat");
    res->clear();

    res = calculator->evalValidPct(a);
    res->saveMatrix("evalvalidpct.mat");
    res->clear();

    res = calculator->evalAbsSum(a);
    res->saveMatrix("evalabssum.mat");
    res->clear();

    res = calculator->evalMean(a);
    res->saveMatrix("evalmean.mat");
    res->clear();

    res = calculator->evalVariance(a);
    res->saveMatrix("evalvariance.mat");
    res->clear();

    res = calculator->evalInnerProduction(a, b);
    res->saveMatrix("evalinnerproduction.mat");
    res->clear();

    res = calculator->evalCovariance(a, b);
    res->saveMatrix("evalcovariance.mat");
    res->clear();

    res = calculator->evalCorrelation(a, b);
    res->saveMatrix("evalcorrelation.mat");
    res->clear();

    res = calculator->evalBeta(a, b);
    res->saveMatrix("evalbeta.mat");
    res->clear();

    res = calculator->diag(a);
    res->saveMatrix("diag.mat");
    res->clear();

    res = calculator->cumSum(a);
    res->saveMatrix("cumsum.mat");
    res->clear();


    return 0;
}	