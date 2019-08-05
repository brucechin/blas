#########################################################################
# File Name: correctness_test.sh
# Author: ma6174
# mail: ma6174@163.com
# Created Time: Tue Jul 30 18:35:02 2019
#########################################################################
#!/bin/bash
g++  correctness_test.cpp ../MatrixCalculator.cc -O3 -std=c++11 -o test -lpthread -lopenblas
./test
cp a.mat b.mat la.mat lb.mat ../java
