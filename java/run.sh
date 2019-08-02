#########################################################################
# File Name: run.sh
# Author: ma6174
# mail: ma6174@163.com
# Created Time: Mon Jul 22 11:44:26 2019
#########################################################################
#!/bin/bash

javac Matrix.java
javah -cp /home/liankeqin/ blas.java.Matrix
g++ -I /opt/Java/jdk/include/ -I /opt/Java/jdk/include/linux/ -shared -fPIC blas_java_Matrix.cc -o libMatrix.so
sudo cp libMatrix.so /usr/lib/
javac Matrix.java
java -Xmx128M blas.java.Matrix
