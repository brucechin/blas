#include "gtest/gtest.h"
#include "cblas.h"
#include "../MatrixCalculator.h"
#include <fstream>
#include<iostream>
using namespace std;

int rangeMin = -10000;
int rangeMax = 10000;
int matSize = 1000;
MatrixFactory* factory = new MatrixFactory();
MatrixCalculator* calculator = new MatrixCalculator();
string matrix_a = "a.mat";
string matrix_b = "b.mat";
string logicMatrix_a = "la.mat";
string logicMatrix_b = "lb.mat";

TEST(Matrix, MatrixIO) {
	Matrix* a = factory->getInstanceOfRandomMatrix(matSize, matSize, rangeMin, rangeMax);
	a->saveMatrix("MatrixIOTest.mat");
	Matrix* b = new Matrix();
	b->readMatrix("MatrixIOTest.mat");
	EXPECT_TRUE(a->compareMatrix(b));
	delete a;
	delete b;
}

//TODO save the res matrix as binary file to compare with python and java

//TODO all functions with "num" parameter are not tested here

TEST(MatrixCalculator, Add) {
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res = calculator->add(a, b);
	EXPECT_TRUE(res->compareMatrix(calculator->add(b, a)));
}

TEST(MatrixCalculator, Sub){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res = calculator->sub(a, b);
	EXPECT_TRUE(a->compareMatrix(calculator->add(b, res)));
	a->clear();
	b->clear();
	res->clear();
}


TEST(MatrixCalculator, Mul){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res = calculator->mul(a, b);
	EXPECT_TRUE(a->compareMatrix(calculator->div(res, b)));//FIX this is not exactly the same because of precision issue.
	a->clear();
	b->clear();
	res->clear();
}


TEST(MatrixCalculator, Div){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res = calculator->div(a, b);
	EXPECT_TRUE(a->compareMatrix(calculator->mul(res, b)));
	a->clear();
	b->clear();
	res->clear();
}

TEST(MatrixCalculator, MatrixMul){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res = calculator->matrixMul(a, b);
	a->clear();
	b->clear();
	res->clear();
}


TEST(MatrixCalculator, Max){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res = calculator->max(a, b);
	EXPECT_TRUE(res->compareMatrix(calculator->max(b, a)));
	a->clear();
	b->clear();
	res->clear();
}

TEST(MatrixCalculator, Min){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res = calculator->min(a, b);
	EXPECT_TRUE(res->compareMatrix(calculator->min(b, a)));
	a->clear();
	b->clear();
	res->clear();
}

TEST(MatrixCalculator, BiggerAndSmaller){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	LogicMatrix* res1 = calculator->bigger(a, b);
	LogicMatrix* res2 = calculator->smaller(b, a);
	EXPECT_TRUE(res1->compareMatrix(res2));
	a->clear();
	b->clear();
	res->clear();
}

TEST(MatrixCalculator, Equal){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	LogicMatrix* res1 = calculator->equal(a, b);
	LogicMatrix* res2 = calculator->equal(a, a);
	a->clear();
	b->clear();
	res->clear();
}

TEST(MatrixCalculator, Condition){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	LogicMatrix* la = new LogicMatrix();
	la->readMatrix(logicMatrix_a);
	LogicMatrix* lb = new LogicMatrix();
	lb->readMatrix(logicMatrix_b);

	Matrix* res1 = calculator->condition(la, a, b);
	Matrix* res2 = calculator->condition(lb, a, b);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

// TEST(MatrixCalculator, Rank){
	
// 	Matrix* a = new Matrix();
// 	a->readMatrix(matrix_a);
// 	Matrix* b = new Matrix();
// 	b->readMatrix(matrix_b);
// 	Matrix* res1 = calculator->rank(a);
// 	Matrix* res2 = calculator->rank(b);
	
// }
//FIX core dump detected

TEST(MatrixCalculator, Round){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->round(a);
	Matrix* res2 = calculator->round(b);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, Floor){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->floor(a);
	Matrix* res2 = calculator->floor(b);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, Abs){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->abs(a);
	Matrix* res2 = calculator->abs(b);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, Minus){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->minus(a);
	Matrix* res2 = calculator->minus(b);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, Sqrt){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->sqrt(a);
	Matrix* res2 = calculator->sqrt(b);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, Log){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->log(a);
	Matrix* res2 = calculator->log(b);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, Exp){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->exp(a);
	Matrix* res2 = calculator->exp(b);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, Sign){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->sign(a);
	Matrix* res2 = calculator->sign(b);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, Inverse){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->inverse(a);
	Matrix* res2 = calculator->inverse(b);
	EXPECT_TRUE(a->compareMatrix(calculator->inverse(res1))); //FIX after two inverse we can not get the original one
	EXPECT_TRUE(b->compareMatrix(calculator->inverse(res2)));
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, SignedPow){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->signedpow(a, 3);
	Matrix* res2 = calculator->signedpow(b, 3);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, Shift){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->shift(a, 100);
	Matrix* res2 = calculator->shift(b ,100);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, Delay){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->delay(a, 100);
	Matrix* res2 = calculator->delay(b ,100);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, Delta){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->delta(a, 100);
	Matrix* res2 = calculator->delta(b ,100);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, Ratio){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->ratio(a, 100);
	Matrix* res2 = calculator->ratio(b ,100);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, Sum){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->sum(a, 100);
	Matrix* res2 = calculator->sum(b ,100);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, Product){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->product(a, 100);
	Matrix* res2 = calculator->product(b ,100);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, TSMax){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->tsMax(a, 100);
	Matrix* res2 = calculator->tsMax(b ,100);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, TSMin){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->tsMin(a, 100);
	Matrix* res2 = calculator->tsMin(b ,100);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, TSArgMax){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->tsArgmax(a, 100);
	Matrix* res2 = calculator->tsArgmax(b ,100);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, TSArgMin){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->tsArgmin(a, 100);
	Matrix* res2 = calculator->tsArgmin(b ,100);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, TSRank){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->tsRank(a, 100);
	Matrix* res2 = calculator->tsRank(b ,100);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, TSMean){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->tsMean(a, 100);
	Matrix* res2 = calculator->tsMean(b ,100);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, TSStd){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->tsStd(a, 100);
	Matrix* res2 = calculator->tsStd(b ,100);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, TSSkewness){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->tsSkewness(a, 100);
	Matrix* res2 = calculator->tsSkewness(b ,100);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, TSKurtosis){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->tsKurtosis(a, 100);
	Matrix* res2 = calculator->tsKurtosis(b ,100);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}	

TEST(MatrixCalculator, TSCov){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->tsCov(a, b, 100);
	Matrix* res2 = calculator->tsCov(b, a, 100);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, TSCorr){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->tsCorr(a, b, 100);
	Matrix* res2 = calculator->tsCorr(b, a, 100);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, TSCountTrue){
	
	LogicMatrix* a = new LogicMatrix();
	a->readMatrix(matrix_a);
	LogicMatrix* b = new LogicMatrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->tsCountTrue(a, 100);
	Matrix* res2 = calculator->tsCountTrue(b, 100);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, TSCountConsecutiveTrue){
	
	LogicMatrix* a = new LogicMatrix();
	a->readMatrix(matrix_a);
	LogicMatrix* b = new LogicMatrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->tsCountConsecutiveTrue(a, 100);
	Matrix* res2 = calculator->tsCountConsecutiveTrue(b, 100);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, DecayLinear){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->decayLinear(a, 100);
	Matrix* res2 = calculator->decayLinear(b, 100);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, DecayExp){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->decayExponential(a, 100);
	Matrix* res2 = calculator->decayExponential(b, 100);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, Neutralize){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->neutralize(a);
	Matrix* res2 = calculator->neutralize(b);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, Mean){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->mean(a);
	Matrix* res2 = calculator->mean(b);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, Unify){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->unify(a);
	Matrix* res2 = calculator->unify(b);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, EvalValidPct){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->evalValidPct(a);
	Matrix* res2 = calculator->evalValidPct(b);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, EvalAbsSum){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->evalAbsSum(a);
	Matrix* res2 = calculator->evalAbsSum(b);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, EvalMean){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->evalMean(a);
	Matrix* res2 = calculator->evalMean(b);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, EvalVariance){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->evalVariance(a);
	Matrix* res2 = calculator->evalVariance(b);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, EvalCov){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->evalCovariance(a, b);
	Matrix* res2 = calculator->evalCovariance(b, a);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, EvalCorr){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->evalCorrelation(a, b);
	Matrix* res2 = calculator->evalCorrelation(b, a);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, EvalInnerProduction){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->evalInnerProduction(a, b);
	Matrix* res2 = calculator->evalInnerProduction(b, a);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, Det){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	double res1 = calculator->Det(a, matSize);
	double res2 = calculator->Det(b, matSize);
	a->clear();
	b->clear();

}


TEST(MatrixCalculator, Treat){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	double res1 = calculator->treat(a);
	double res2 = calculator->treat(b);
	a->clear();
	b->clear();

}

TEST(MatrixCalculator, Inv){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->inv(a);
	Matrix* res2 = calculator->inv(b);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, Diag){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->diag(a);
	Matrix* res2 = calculator->diag(b);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, InverseDiag){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->inverseDiag(a);
	Matrix* res2 = calculator->inverseDiag(b);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, EvalBeta){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->evalBeta(a, b);
	Matrix* res2 = calculator->evalBeta(b, a);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

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
	testing::InitGoogleTest();
	RUN_ALL_TESTS();
    return 0;
}		    
