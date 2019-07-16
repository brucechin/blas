/*
 * @Author: lianke.qin@gmail.com 
 * @Date: 2019-07-16 11:03:57 
 * @Goal: matrix computation speed test
 * @Last Modified by: lianke.qin@gmail.com
 * @Last Modified time: 2019-07-16 11:31:15
 */


#include "../LogicMatrix.h"
#include "../Matrix.h"
#include "../MatrixCalculator.h"
#include "../MatrixFactory.h"
#include<vector>
#include<iostream>
#include<fstream>
#include<string>

using namespace std;

int main(){
    Matrix a = new Matrix(1000, 1000, -1000, 1000);
    Matrix b = new Matrix(1000, 1000, -1000, 1000);
    ofstream f1, f2;
    f1.open("1k-a.txt",ios::binary);
    f1.write((char *)&a, sizeof(a));
    f2.open("1k-b.txt",ios::binary);
    f2.write((char *)&b, sizeof(b));
    f1.close();
    f2.close();
    
}