#########################################################################
# File Name: correctness_test.sh
# Author: lianke qin
# mail: lianke.qin@gmail.com
# Created Time: Tue Jul 30 18:35:02 2019
#########################################################################
#!/bin/bash
INCLUDE_PATH='/usr/local/include/'
LIBRARY_PATH='/usr/local/lib/'
g++ -g -std=c++11 -O3 unit_test.cpp ../src/MatrixCalculator.cc -o unit_test -I $INCLUDE_PATH -L $LIBRARY_PATH -lgtest -lpthread -lopenblas -lgfortran
./unit_test
