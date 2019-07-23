#include "cblas.h"
#include "MatrixCalculator.h"
#include <fstream>
#include<iostream>
using namespace std;

int rangeMin = -100;
int rangeMax = 100;
int matSize = 10;
MatrixFactory* factory = new MatrixFactory();
MatrixCalculator* calculator = new MatrixCalculator();
string file_a = "a.mat";
string file_b = "b.mat";


int main()
{	
	Matrix* a = factory->getInstanceOfRandomMatrix(matSize, matSize, rangeMin, rangeMax);
	a->saveMatrix(file_a);
	Matrix* b = factory->getInstanceOfRandomMatrix(matSize, matSize, rangeMin, rangeMax);
	b->saveMatrix(file_b);

	double res = calculator->Det(a, matSize);

	system("pause");
    return 0;
}
