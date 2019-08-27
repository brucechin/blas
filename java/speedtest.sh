#cd ../tests
#sh correctness_test.sh
#cp *.mat ../java
#rm *.mat
#cd ../java

#run each API in MatrixCalculator and print their time usage
JAVA_HOME='/opt/Java/jdk'
OPENBLAS_HOME=''


g++ -std=c++11 -I $JAVA_HOME/include -I $JAVA_HOME/include/linux -shared -fPIC blas_java_Matrix.cc -o libMatrix.so 
g++ -std=c++11 -I $JAVA_HOME/include -I $JAVA_HOME/include/linux -shared -fPIC blas_java_LogicMatrix.cc -o libLogicMatrix.so 
g++ -std=c++11 -I $JAVA_HOME/include -I $JAVA_HOME/include/linux -shared -fPIC blas_java_MatrixFactory.cc -o libMatrixFactory.so 
g++ -std=c++11 -I $JAVA_HOME/include -I $JAVA_HOME/include/linux -shared -fPIC ../src/MatrixCalculator.cc blas_java_MatrixCalculator.cc -o libMatrixCalculator.so -lpthread -lopenblas

sudo cp lib*.so /usr/lib/

javac Matrix.java
javac LogicMatrix.java
javac MatrixFactory.java
javac MatrixCalculator.java
javac MatrixSpeedTest.java

java blas.java.MatrixSpeedTest
rm *.mat
