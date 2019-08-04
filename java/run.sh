#########################################################################
# File Name: run.sh
# Author: ma6174
# mail: ma6174@163.com
# Created Time: Mon Jul 22 11:44:26 2019
#########################################################################
#!/bin/bash

javac MatrixCalculator.java
javah -cp /home/liankeqin/ blas.java.MatrixCalculator
g++ -I /opt/Java/jdk/include/ -I /opt/Java/jdk/include/linux/ -shared -fPIC blas_java_MatrixCalculator.cc -o libMatrixCalculator.so
sudo cp libMatrixCalculator.so /usr/lib/
javac MatrixCalculator.java
java -Xmx128M blas.java.MatrixCalculator
