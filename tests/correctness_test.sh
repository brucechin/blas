#########################################################################
# File Name: correctness_test.sh
# Author: lianke qin
# mail: lianke.qin@gmail.com
# Created Time: Tue Jul 30 18:35:02 2019
#########################################################################
#!/bin/bash
g++  correctness_test.cpp ../src/MatrixCalculator.cc -O3 -std=c++11 -o test -lpthread -lopenblas
./test
