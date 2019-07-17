

//#include "../LogicMatrix.h"
//#include "../Matrix.h"
#include "../MatrixCalculator.h"
//#include "../MatrixFactory.h"

#include<vector>
#include<iostream>
#include<fstream>
#include<string>

using namespace std;

Matrix* getInstanceOfRandomMatrix(int n, int m, int min, int max){
        srand((unsigned)time(NULL));
        Matrix* res = new Matrix(n, m);
        double* p = &res->value[0];
        for(int i = 0; i < n; i++){
            for(int j = 0; j < m; j++){
                double m1 = (double)(rand()%101)/101.0;
                double m2 = (double)(rand()%(max - min + 1) + min);
                m2 -= 1;
                (*p++) = m1 + m2;
            }
        }
        return res;
    }

int main(){

    Matrix* a = getInstanceOfRandomMatrix(10, 10, -1000, 1000);
    Matrix* b = getInstanceOfRandomMatrix(10, 10, -1000, 1000);
    a->print();
    // ofstream f1, f2;
    // f1.open("1k-a.txt",ios::binary);
    // f1.write((char *)&a, sizeof(a));
    // f2.open("1k-b.txt",ios::binary);
    // f2.write((char *)&b, sizeof(b));
    // f1.close();
    // f2.close();
    system("pause");
}