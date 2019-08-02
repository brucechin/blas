#########################################################################
# File Name: run.sh
# Author: ma6174
# mail: ma6174@163.com
# Created Time: Mon Jul 22 11:44:26 2019
#########################################################################
#!/bin/bash

javac LogicMatrix.java
javah -cp /home/liankeqin/ blas.java.LogicMatrix
g++ -I /opt/Java/jdk/include/ -I /opt/Java/jdk/include/linux/ -shared -fPIC blas_java_LogicMatrix.cc -o libLogicMatrix.so
sudo cp libLogicMatrix.so /usr/lib/
javac LogicMatrix.java
java -Xmx128M blas.java.LogicMatrix
