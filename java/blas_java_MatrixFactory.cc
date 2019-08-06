#include "../Matrix.h"
#include "blas_java_MatrixFactory.h"
#include "../LogicMatrix.h"
#include "../MatrixFactory.h"
#include<iostream>
#include<stdexcept>
#include<unistd.h>
#include<cstring>

JNIEXPORT jlong JNICALL Java_blas_java_MatrixFactory_getInstanceOfEmptyMatrixNative(JNIEnv * env, jclass obj){
    Matrix* res = MatrixFactory::getInstanceOfEmptyMatrix();
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixFactory_getInstanceOfNaNMatrixNative(JNIEnv * env, jclass obj, jint size){
    Matrix* res = MatrixFactory::getInstanceOfNaNMatrix(size);
    return (jlong)res;
}



JNIEXPORT jlong JNICALL Java_blas_java_MatrixFactory_getInstanceOfRandomMatrixNative(JNIEnv * env, jclass obj, jint n, jint m, jint lower, jint upper){
    Matrix* res = MatrixFactory::getInstanceOfRandomMatrix(n, m, lower, upper);
    return (jlong)res;
}

JNIEXPORT jlong JNICALL Java_blas_java_MatrixFactory_getInstanceOfRandomLogicMatrixNative(JNIEnv * env, jclass obj, jint n, jint m){
    LogicMatrix* res = MatrixFactory::getInstanceOfRandomLogicMatrix(n, m);
    return (jlong)res;
}

JNIEXPORT jlong JNICALL Java_blas_java_MatrixFactory_getInstanceOfRandomMatrixWithAbnormalValueNative(JNIEnv * env, jclass obj, jint n, jint m, jint lower, jint upper, jint frequency){
    Matrix* res = MatrixFactory::getInstanceOfRandomMatrixWithAbnormalValue(n, m, lower, upper, frequency);
    return (jlong)res;
}