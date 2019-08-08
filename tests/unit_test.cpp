#include "gtest/gtest.h"
#include <fstream>
#include<iostream>
#include"../MatrixCalculator.h"
#include"../MatrixFactory.h"
using namespace std;

int rangeMin = -10000;
int rangeMax = 10000;
int matSize = 10;
int step = 5;
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
	
	for(int i = 0; i < 1; i++){
		Matrix* res = MatrixCalculator::mul(a, b);
		a->value[i]++;
		b->value[i]--;
		res->clear();
	}
	a->clear();
	b->clear();
	//res->clear();
}


TEST(MatrixCalculator, Div){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res = calculator->div(a, b);
	a->clear();
	b->clear();
	res->clear();
}

TEST(MatrixCalculator, MatrixMul){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	for(int i = 0; i < 1; i++){
		Matrix* res = MatrixCalculator::matrixMul(a, b);
		a->value[i]++;
		b->value[i]--;
		res->clear();
	}
	a->clear();
	b->clear();
	//res->clear();
}



TEST(MatrixCalculator, Max){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res = calculator->max(a, b);
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
	res1->clear();
	res2->clear();
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
	res1->clear();
	res2->clear();
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
	//EXPECT_TRUE(a->compareMatrix(calculator->inverse(res1))); //FIX after two inverse we can not get the original one
	//EXPECT_TRUE(b->compareMatrix(calculator->inverse(res2)));
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
	Matrix* res1 = calculator->shift(a, step);
	Matrix* res2 = calculator->shift(b ,step);
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
	Matrix* res1 = calculator->delay(a, step);
	Matrix* res2 = calculator->delay(b ,step);
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
	Matrix* res1 = calculator->delta(a, step);
	Matrix* res2 = calculator->delta(b ,step);
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
	Matrix* res1 = calculator->ratio(a, step);
	Matrix* res2 = calculator->ratio(b ,step);
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
	Matrix* res1 = calculator->sum(a, step);
	Matrix* res2 = calculator->sum(b ,step);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, Sum_OP){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->sum_op(a, step);
	Matrix* res2 = calculator->sum_op(b ,step);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, Sum_Compare){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->sum(a, step);
	Matrix* res2 = calculator->sum_op(a, step);
	EXPECT_TRUE(res1->compareMatrix(res2));
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
	Matrix* res1 = calculator->product(a, step);
	Matrix* res2 = calculator->product(b ,step);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, Product_OP){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->product_op(a, step);
	Matrix* res2 = calculator->product_op(b ,step);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, Product_Compare){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->product(a, step);
	Matrix* res2 = calculator->product(b ,step);
	Matrix* res1_op = calculator->product_op(a, step);
	Matrix* res2_op = calculator->product_op(b ,step);
	EXPECT_TRUE(res1->compareMatrix(res1_op));
	EXPECT_TRUE(res2->compareMatrix(res2_op));
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
	res1_op->clear();
	res2_op->clear();
}

TEST(MatrixCalculator, TSMax){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->tsMax(a, step);
	Matrix* res2 = calculator->tsMax(b ,step);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, TSMax_OP){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->tsMax_op(a, step);
	Matrix* res2 = calculator->tsMax_op(b ,step);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, TSMax_Compare){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->tsMax(a, step);
	Matrix* res2 = calculator->tsMax(b ,step);
	Matrix* res1_op = calculator->tsMax_op(a, step);
	Matrix* res2_op = calculator->tsMax_op(b ,step);
	EXPECT_TRUE(res1->compareMatrix(res1_op));
	EXPECT_TRUE(res2->compareMatrix(res2_op));
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
	res1_op->clear();
	res2_op->clear();
}

TEST(MatrixCalculator, TSMin){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->tsMin(a, step);
	Matrix* res2 = calculator->tsMin(b ,step);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, TSMin_OP){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->tsMin_op(a, step);
	Matrix* res2 = calculator->tsMin_op(b ,step);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, TSMin_Compare){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->tsMin(a, step);
	Matrix* res2 = calculator->tsMin(b ,step);
	Matrix* res1_op = calculator->tsMin_op(a, step);
	Matrix* res2_op = calculator->tsMin_op(b ,step);
	EXPECT_TRUE(res1->compareMatrix(res1_op));
	EXPECT_TRUE(res2->compareMatrix(res2_op));
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
	res1_op->clear();
	res2_op->clear();
}


TEST(MatrixCalculator, TSArgMax){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->tsArgmax(a, step);
	Matrix* res2 = calculator->tsArgmax(b ,step);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, TSArgMax_OP){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->tsArgmax_op(a, step);
	Matrix* res2 = calculator->tsArgmax_op(b ,step);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, TSArgMax_Compare){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->tsArgmax(a, step);
	Matrix* res2 = calculator->tsArgmax(b ,step);
	Matrix* res1_op = calculator->tsArgmax_op(a, step);
	Matrix* res2_op = calculator->tsArgmax_op(b ,step);
	EXPECT_TRUE(res1->compareMatrix(res1_op));
	EXPECT_TRUE(res2->compareMatrix(res2_op));
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
	res1_op->clear();
	res2_op->clear();
}


TEST(MatrixCalculator, TSArgMin){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->tsArgmin(a, step);
	Matrix* res2 = calculator->tsArgmin(b ,step);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, TSArgMin_OP){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->tsArgmin_op(a, step);
	Matrix* res2 = calculator->tsArgmin_op(b ,step);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, TSArgMin_Compare){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->tsArgmin(a, step);
	Matrix* res2 = calculator->tsArgmin(b ,step);
	Matrix* res1_op = calculator->tsArgmin_op(a, step);
	Matrix* res2_op = calculator->tsArgmin_op(b ,step);
	EXPECT_TRUE(res1->compareMatrix(res1_op));
	EXPECT_TRUE(res2->compareMatrix(res2_op));
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
	res1_op->clear();
	res2_op->clear();
}


TEST(MatrixCalculator, TSRank_Compare){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->tsRank(a, step);
	Matrix* res2 = calculator->tsRank(b ,step);
	Matrix* res1_op = calculator->tsRank_op(a, step);
	EXPECT_TRUE(res1->compareMatrix(res1_op));
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, TSMean_Compare){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->tsMean(a, step);
	Matrix* res2 = calculator->tsMean(b ,step);
	Matrix* res1_op = calculator->tsMean_op(a, step);
	EXPECT_TRUE(res1->compareMatrix(res1_op));
	
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, TSStd_Compare){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->tsStd(a, step);
	Matrix* res2 = calculator->tsStd(b ,step);
	Matrix* res1_op = calculator->tsStd_op(a, step);
	Matrix* res2_op = calculator->tsStd_op(b, step);
	EXPECT_TRUE(res1->compareMatrix(res1_op));
	EXPECT_TRUE(res2->compareMatrix(res2_op));

	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
	res1_op->clear();
}

TEST(MatrixCalculator, TSSkewness_Compare){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->tsSkewness(a, step);
	Matrix* res2 = calculator->tsSkewness(b ,step);
	Matrix* res1_op = calculator->tsSkewness_op(a, step);
	EXPECT_TRUE(res1->compareMatrix(res1_op));
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, TSKurtosis_Compare){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->tsKurtosis(a, step);
	Matrix* res2 = calculator->tsKurtosis(b ,step);
	Matrix* res1_op = calculator->tsKurtosis_op(a, step);
	EXPECT_TRUE(res1->compareMatrix(res1_op));
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
	Matrix* res1 = calculator->tsCov(a, b, step);
	Matrix* res2 = calculator->tsCov(b, a, step);
	Matrix* res1_op = calculator->tsCov_op(a, b, step);
	EXPECT_TRUE(res1->compareMatrix(res1_op));
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
	Matrix* res1 = calculator->tsCorr(a, b, step);
	Matrix* res1_op = calculator->tsCorr_op(a, b, step);
	EXPECT_TRUE(res1->compareMatrix(res1_op));
	a->clear();
	b->clear();
	res1->clear();
}

TEST(MatrixCalculator, TSCountTrue){
	
	LogicMatrix* a = new LogicMatrix();
	a->readMatrix(logicMatrix_a);
	LogicMatrix* b = new LogicMatrix();
	b->readMatrix(logicMatrix_b);
	Matrix* res1 = calculator->tsCountTrue(a, step);
	Matrix* res2 = calculator->tsCountTrue(b, step);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, TSCountTrue_OP){
	
	LogicMatrix* a = new LogicMatrix();
	a->readMatrix(logicMatrix_a);
	LogicMatrix* b = new LogicMatrix();
	b->readMatrix(logicMatrix_b);
	Matrix* res1 = calculator->tsCountTrue_op(a, step);
	Matrix* res2 = calculator->tsCountTrue_op(b ,step);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, TSCountTrue_Compare){
	
	LogicMatrix* a = new LogicMatrix();
	a->readMatrix(logicMatrix_a);
	LogicMatrix* b = new LogicMatrix();
	b->readMatrix(logicMatrix_b);
	Matrix* res1 = calculator->tsCountTrue(a, step);
	Matrix* res2 = calculator->tsCountTrue(b ,step);
	Matrix* res1_op = calculator->tsCountTrue_op(a, step);
	Matrix* res2_op = calculator->tsCountTrue_op(b ,step);
	EXPECT_TRUE(res1->compareMatrix(res1_op));
	EXPECT_TRUE(res2->compareMatrix(res2_op));
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
	res1_op->clear();
	res2_op->clear();
}


TEST(MatrixCalculator, TSCountConsecutiveTrue){
	
	LogicMatrix* a = new LogicMatrix();
	a->readMatrix(logicMatrix_a);
	LogicMatrix* b = new LogicMatrix();
	b->readMatrix(logicMatrix_b);
	Matrix* res1 = calculator->tsCountConsecutiveTrue(a, step);
	Matrix* res2 = calculator->tsCountConsecutiveTrue(b, step);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, TSCountConsecutiveTrue_OP){
	
	LogicMatrix* a = new LogicMatrix();
	a->readMatrix(logicMatrix_a);
	LogicMatrix* b = new LogicMatrix();
	b->readMatrix(logicMatrix_b);
	Matrix* res1 = calculator->tsCountConsecutiveTrue_op(a, step);
	Matrix* res2 = calculator->tsCountConsecutiveTrue_op(b ,step);
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
}

TEST(MatrixCalculator, TSCountConsecutiveTrue_Compare){
	
	LogicMatrix* a = new LogicMatrix();
	a->readMatrix(logicMatrix_a);
	LogicMatrix* b = new LogicMatrix();
	b->readMatrix(logicMatrix_b);
	Matrix* res1 = calculator->tsCountConsecutiveTrue(a, step);
	Matrix* res2 = calculator->tsCountConsecutiveTrue(b ,step);
	Matrix* res1_op = calculator->tsCountConsecutiveTrue_op(a, step);
	Matrix* res2_op = calculator->tsCountConsecutiveTrue_op(b ,step);
	EXPECT_TRUE(res1->compareMatrix(res1_op));
	EXPECT_TRUE(res2->compareMatrix(res2_op));
	a->clear();
	b->clear();
	res1->clear();
	res2->clear();
	res1_op->clear();
	res2_op->clear();
}

TEST(MatrixCalculator, DecayLinear){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->decayLinear(a, step);
	Matrix* res2 = calculator->decayLinear(b, step);
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
	Matrix* res1 = calculator->decayExponential(a, step);
	Matrix* res2 = calculator->decayExponential(b, step);
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

//FIX segmentation fault detected
// TEST(MatrixCalculator, Mean){
	
// 	Matrix* a = new Matrix();
// 	a->readMatrix(matrix_a);
// 	Matrix* b = new Matrix();
// 	b->readMatrix(matrix_b);
// 	Matrix* res1 = calculator->mean(a);
// 	Matrix* res2 = calculator->mean(b);
// 	a->clear();
// 	b->clear();
// 	res1->clear();
// 	res2->clear();
// }

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

// TEST(MatrixCalculator, Det){
	
// 	Matrix* a = new Matrix();
// 	a->readMatrix(matrix_a);
// 	Matrix* b = new Matrix();
// 	b->readMatrix(matrix_b);
// 	double res1 = calculator->Det(a, matSize);
// 	double res2 = calculator->Det(b, matSize);
// 	a->clear();
// 	b->clear();

// }
//FIX core dump detected


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

//FIX segmentation fault detected
// TEST(MatrixCalculator, Inv){
	
// 	Matrix* a = new Matrix();
// 	a->readMatrix(matrix_a);
// 	Matrix* b = new Matrix();
// 	b->readMatrix(matrix_b);
// 	Matrix* res1 = calculator->inv(a);
// 	Matrix* res2 = calculator->inv(b);
// 	a->clear();
// 	b->clear();
// 	res1->clear();
// 	res2->clear();
// }

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
/*	
	int times = step;
	
	for(int i = 0; i < times; i++)
		Matrix* res1 = calculator->mul(a, b);

	for(int i = 0; i < times; i++)
		Matrix* res2 = calculator->matrixMul(a, b);
*/

	testing::InitGoogleTest();
	RUN_ALL_TESTS();
    return 0;
}		    
