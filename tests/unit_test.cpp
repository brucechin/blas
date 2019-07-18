#include "gtest/gtest.h"
#include "cblas.h"
#include "../MatrixCalculator.h"
#include <fstream>
#include<iostream>
using namespace std;

int rangeMin = -10000;
int rangeMax = 10000;
int matSize = 2000;
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

}

TEST(MatrixCalculator, MulAndDiv){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res = calculator->mul(a, b);
	EXPECT_TRUE(a->compareMatrix(calculator->div(res, b)));

}

TEST(MatrixCalculator, Mul){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res = calculator->mul(a, b);
	EXPECT_TRUE(a->compareMatrix(calculator->div(res, b)));
	
}


TEST(MatrixCalculator, Div){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res = calculator->div(a, b);
	EXPECT_TRUE(a->compareMatrix(calculator->mul(res, b)));
	
}

TEST(MatrixCalculator, MatrixMul){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res = calculator->matrixMul(a, b);
	
}


TEST(MatrixCalculator, Max){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res = calculator->max(a, b);
	EXPECT_TRUE(res->compareMatrix(calculator->max(b, a)));
	
}

TEST(MatrixCalculator, Min){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res = calculator->min(a, b);
	EXPECT_TRUE(res->compareMatrix(calculator->min(b, a)));
	
}

TEST(MatrixCalculator, BiggerAndSmaller){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	LogicMatrix* res1 = calculator->bigger(a, b);
	LogicMatrix* res2 = calculator->smaller(b, a);
	EXPECT_TRUE(res1->compareMatrix(res2));
	
}

TEST(MatrixCalculator, Equal){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	LogicMatrix* res1 = calculator->equal(a, b);
	LogicMatrix* res2 = calculator->equal(a, a);

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

}

TEST(MatrixCalculator, Rank){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->rank(a);
	Matrix* res2 = calculator->rank(b);

}

TEST(MatrixCalculator, Round){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->round(a);
	Matrix* res2 = calculator->round(b);

}

TEST(MatrixCalculator, Floor){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->floor(a);
	Matrix* res2 = calculator->floor(b);

}

TEST(MatrixCalculator, Abs){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->abs(a);
	Matrix* res2 = calculator->abs(b);

}

TEST(MatrixCalculator, Minus){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->minus(a);
	Matrix* res2 = calculator->minus(b);

}

TEST(MatrixCalculator, Sqrt){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->sqrt(a);
	Matrix* res2 = calculator->sqrt(b);

}

TEST(MatrixCalculator, Log){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->log(a);
	Matrix* res2 = calculator->log(b);

}

TEST(MatrixCalculator, Exp){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->exp(a);
	Matrix* res2 = calculator->exp(b);

}

TEST(MatrixCalculator, Sign){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->sign(a);
	Matrix* res2 = calculator->sign(b);

}

TEST(MatrixCalculator, Inverse){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->inverse(a);
	Matrix* res2 = calculator->inverse(b);
	EXPECT_TRUE(a->compareMatrix(calculator->inverse(res1)));
	EXPECT_TRUE(b->compareMatrix(calculator->inverse(res2)));
}

TEST(MatrixCalculator, SignedPow){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->signedpow(a, 3);
	Matrix* res2 = calculator->signedpow(b, 3);

}

TEST(MatrixCalculator, Shift){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->shift(a, 100);
	Matrix* res2 = calculator->shift(b ,100);

}

TEST(MatrixCalculator, Delay){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->delay(a, 100);
	Matrix* res2 = calculator->delay(b ,100);

}

TEST(MatrixCalculator, Delta){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->delta(a, 100);
	Matrix* res2 = calculator->delta(b ,100);

}

TEST(MatrixCalculator, Ratio){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->ratio(a, 100);
	Matrix* res2 = calculator->ratio(b ,100);

}

TEST(MatrixCalculator, Sum){
	
	Matrix* a = new Matrix();
	a->readMatrix(matrix_a);
	Matrix* b = new Matrix();
	b->readMatrix(matrix_b);
	Matrix* res1 = calculator->sum(a, 100);
	Matrix* res2 = calculator->sum(b ,100);

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
