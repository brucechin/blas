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
string file_a = "a.mat";
string file_b = "b.mat";

TEST(Matrix, MatrixIO) {
	Matrix* a = factory->getInstanceOfRandomMatrix(matSize, matSize, rangeMin, rangeMax);
	a->saveMatrix("MatrixIOTest.mat");
	Matrix* b = new Matrix();
	b->readMatrix("MatrixIOTest.mat");
	EXPECT_TRUE(a->compareMatrix(b));
	delete a;
	delete b;
}


TEST(MatrixCalculator, Add) {
	Matrix* a = new Matrix();
	a->readMatrix(file_a);
	Matrix* b = new Matrix();
	b->readMatrix(file_b);
	Matrix* res = calculator->add(a, b);
	EXPECT_TRUE(res->compareMatrix(calculator->add(b, a)));
}
int main()
{	
        Matrix* a = factory->getInstanceOfRandomMatrix(matSize, matSize, rangeMin, rangeMax);
	a->saveMatrix(file_a);
	Matrix* b = factory->getInstanceOfRandomMatrix(matSize, matSize, rangeMin, rangeMax);
	b->saveMatrix(file_b);
	testing::InitGoogleTest();
	RUN_ALL_TESTS();
        return 0;
}
						    
