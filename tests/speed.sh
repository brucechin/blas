INCLUDE_PATH='/usr/local/include/'
LIBRARY_PATH='/usr/local/lib/'
g++ -g -std=c++11 -O3 speedtest.cpp ../src/MatrixCalculator.cc -o speedtest -I $INCLUDE_PATH -L $LIBRARY_PATH -lgtest -lpthread -lopenblas -lgfortran
./speedtest
