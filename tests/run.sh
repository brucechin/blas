g++ -std=c++11 -O3 unit_test.cpp ../MatrixCalculator.cc -o unit_test -I /usr/local/include -L /usr/local/lib -lgtest -lpthread -lopenblas -lgfortran
./unit_test
rm unit_test
