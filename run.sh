g++ matrix.cc -O3 -o matrix -I/opt/OpenBLAS/include/ -L/opt/OpenBLAS/lib -lopenblas -lpthread -lgfortran
OMP_NUM_THREAD=1 ./matrix