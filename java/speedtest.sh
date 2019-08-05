cd ../tests
sh correctness_test.sh
cp *.mat ../java
rm *.mat
cd ../java

g++ -I /opt/Java/jdk/include -I /opt/Java/jdk/include/linux -shared -fPIC blas_java_Matrix.cc -o libMatrix.so 
g++ -I /opt/Java/jdk/include -I /opt/Java/jdk/include/linux -shared -fPIC blas_java_LogicMatrix.cc -o libLogicMatrix.so 
g++ -I /opt/Java/jdk/include -I /opt/Java/jdk/include/linux -shared -fPIC ../MatrixCalculator.cc blas_java_MatrixCalculator.cc -o libMatrixCalculator.so -lpthread -lopenblas

sudo cp lib*.so /usr/lib/

javac Matrix.java
javac LogicMatrix.java
javac MatrixCalculator.java
javac MatrixSpeedTest.java

java blas.java.MatrixSpeedTest
rm *.mat