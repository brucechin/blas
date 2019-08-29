#include <fstream>
#include<iostream>
#include"../src/MatrixCalculator.h"
#include"../src/MatrixFactory.h"
#include <chrono> 
using namespace std;
using namespace chrono;
int rangeMin = -10000;
int rangeMax = 10000;
int matSize = 3000;
int step = 50;
int times = 10;

MatrixFactory* factory = new MatrixFactory();
MatrixCalculator* calculator = new MatrixCalculator();
string matrix_a = "a.mat";
string matrix_b = "b.mat";
string logicMatrix_a = "la.mat";
string logicMatrix_b = "lb.mat";

int main(){
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

    auto start = system_clock::now();
    for(int i = 0; i < times; i++){
        res = MatrixCalculator::add(a, b);
        delete res;
    }
    auto end   = system_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    cout <<  "add : " << double(duration.count()) * 1000 *microseconds::period::num / microseconds::period::den << " ms"<< endl;

    start = system_clock::now();
    for(int i = 0; i < times; i++){
        res = MatrixCalculator::sub(a, b);
        delete res;
    }
    end   = system_clock::now();
    duration = duration_cast<microseconds>(end - start);
    cout <<  "sub : " << double(duration.count()) * 1000 *microseconds::period::num / microseconds::period::den << " ms"<< endl;

    start = system_clock::now();
    for(int i = 0; i < times; i++){
        res = MatrixCalculator::mul(a, b);
        delete res;
    }
    end   = system_clock::now();
    duration = duration_cast<microseconds>(end - start);
    cout <<  "mul : " << double(duration.count()) * 1000 *microseconds::period::num / microseconds::period::den << " ms"<< endl;

    start = system_clock::now();
    for(int i = 0; i < times; i++){
        res = MatrixCalculator::div(a, b);
        delete res;
    }
    end   = system_clock::now();
    duration = duration_cast<microseconds>(end - start);
    cout <<  "div : " << double(duration.count()) * 1000 *microseconds::period::num / microseconds::period::den << " ms"<< endl;

    start = system_clock::now();
    for(int i = 0; i < times; i++){
        res = MatrixCalculator::max(a, b);
        delete res;
    }
    end   = system_clock::now();
    duration = duration_cast<microseconds>(end - start);
    cout <<  "max : " << double(duration.count()) * 1000 *microseconds::period::num / microseconds::period::den << " ms"<< endl;


    start = system_clock::now();
    for(int i = 0; i < times; i++){
        lres = MatrixCalculator::bigger(a, b);
        delete lres;
    }
    end   = system_clock::now();
    duration = duration_cast<microseconds>(end - start);
    cout <<  "bigger : " << double(duration.count()) * 1000 *microseconds::period::num / microseconds::period::den << " ms"<< endl;


    start = system_clock::now();
    for(int i = 0; i < times; i++){
        lres = MatrixCalculator::matAnd(la, lb);
        delete lres;
    }
    end   = system_clock::now();
    duration = duration_cast<microseconds>(end - start);
    cout <<  "and : " << double(duration.count()) * 1000 *microseconds::period::num / microseconds::period::den << " ms"<< endl;


    start = system_clock::now();
    for(int i = 0; i < times; i++){
        res = MatrixCalculator::condition(la, a, b);
        delete res;
    }
    end   = system_clock::now();
    duration = duration_cast<microseconds>(end - start);
    cout <<  "condition : " << double(duration.count()) * 1000 *microseconds::period::num / microseconds::period::den << " ms"<< endl;


    start = system_clock::now();
    for(int i = 0; i < times; i++){
        res = MatrixCalculator::abs(a);
        delete res;
    }
    end   = system_clock::now();
    duration = duration_cast<microseconds>(end - start);
    cout <<  "abs : " << double(duration.count()) * 1000 *microseconds::period::num / microseconds::period::den << " ms"<< endl;

    start = system_clock::now();
    for(int i = 0; i < times; i++){
        res = MatrixCalculator::sqrt(a);
        delete res;
    }
    end   = system_clock::now();
    duration = duration_cast<microseconds>(end - start);
    cout <<  "sqrt : " << double(duration.count()) * 1000 *microseconds::period::num / microseconds::period::den << " ms"<< endl;


    start = system_clock::now();
    for(int i = 0; i < times; i++){
        res = MatrixCalculator::shift(a, step);
        delete res;
    }
    end   = system_clock::now();
    duration = duration_cast<microseconds>(end - start);
    cout <<  "shift : " << double(duration.count()) * 1000 *microseconds::period::num / microseconds::period::den << " ms"<< endl;

    start = system_clock::now();
    for(int i = 0; i < times; i++){
        res = MatrixCalculator::ratio(a, step);
        delete res;
    }
    end   = system_clock::now();
    duration = duration_cast<microseconds>(end - start);
    cout <<  "ratio : " << double(duration.count()) * 1000 *microseconds::period::num / microseconds::period::den << " ms"<< endl;

        start = system_clock::now();
    for(int i = 0; i < times; i++){
        res = MatrixCalculator::sum(a, step);
        delete res;
    }
    end   = system_clock::now();
    duration = duration_cast<microseconds>(end - start);
    cout <<  "sum : " << double(duration.count()) * 1000 *microseconds::period::num / microseconds::period::den << " ms"<< endl;

    start = system_clock::now();
    for(int i = 0; i < times; i++){
        res = MatrixCalculator::tsMax(a, step);
        delete res;
    }
    end   = system_clock::now();
    duration = duration_cast<microseconds>(end - start);
    cout <<  "tsMax : " << double(duration.count()) * 1000 *microseconds::period::num / microseconds::period::den << " ms"<< endl;

        start = system_clock::now();
    for(int i = 0; i < times; i++){
        res = MatrixCalculator::tsRank(a, step);
        delete res;
    }
    end   = system_clock::now();
    duration = duration_cast<microseconds>(end - start);
    cout <<  "tsRank : " << double(duration.count()) * 1000 *microseconds::period::num / microseconds::period::den << " ms"<< endl;

    start = system_clock::now();
    for(int i = 0; i < times; i++){
        res = MatrixCalculator::tsStd(a, step);
        delete res;
    }
    end   = system_clock::now();
    duration = duration_cast<microseconds>(end - start);
    cout <<  "tsStd : " << double(duration.count()) * 1000 *microseconds::period::num / microseconds::period::den << " ms"<< endl;


    start = system_clock::now();
    for(int i = 0; i < times; i++){
        res = MatrixCalculator::tsMean(a, step);
        delete res;
    }
    end   = system_clock::now();
    duration = duration_cast<microseconds>(end - start);
    cout <<  "tsMean : " << double(duration.count()) * 1000 *microseconds::period::num / microseconds::period::den << " ms"<< endl;

    start = system_clock::now();
    for(int i = 0; i < times; i++){
        res = MatrixCalculator::tsCov(a, b, step);
        delete res;
    }
    end   = system_clock::now();
    duration = duration_cast<microseconds>(end - start);
    cout <<  "tscov : " << double(duration.count()) * 1000 *microseconds::period::num / microseconds::period::den << " ms"<< endl;

    start = system_clock::now();
    for(int i = 0; i < times; i++){
        res = MatrixCalculator::tsSkewness(a, step);
        delete res;
    }
    end   = system_clock::now();
    duration = duration_cast<microseconds>(end - start);
    cout <<  "tsSkewness : " << double(duration.count()) * 1000 *microseconds::period::num / microseconds::period::den << " ms"<< endl;

    start = system_clock::now();
    for(int i = 0; i < times; i++){
        res = MatrixCalculator::tsKurtosis(a, step);
        delete res;
    }
    end   = system_clock::now();
    duration = duration_cast<microseconds>(end - start);
    cout <<  "tsKurtosis : " << double(duration.count()) * 1000 *microseconds::period::num / microseconds::period::den << " ms"<< endl;




}