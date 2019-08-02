#include "../LogicMatrix.h"
#include "blas_java_LogicMatrix.h"
#include<iostream>
#include<stdexcept>
#include<unistd.h>
#include<cstring>

JNIEXPORT void JNICALL Java_blas_java_LogicMatrix_setElementNative(JNIEnv *env, jobject jb, jint i, jint j, jboolean val, jlong ptr){
	LogicMatrix* mat = (LogicMatrix*)ptr;
	mat->setElement(i, j, val);
}


JNIEXPORT jboolean JNICALL Java_blas_java_LogicMatrix_getElementNative(JNIEnv *env, jobject jb, jint i, jint j, jlong ptr){
	LogicMatrix* mat = (LogicMatrix*)ptr;
	jboolean res = mat->getElement(i, j);
	return res;
}


JNIEXPORT jlong JNICALL Java_blas_java_LogicMatrix_getRowVectorNative(JNIEnv *env, jobject jb, jint i, jlong ptr){
	LogicMatrix* mat = (LogicMatrix*)ptr;
	LogicMatrix* res = mat->getRowVector(i);
	return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_LogicMatrix_getColVectorNative(JNIEnv *env, jobject jb, jint j, jlong ptr){
	LogicMatrix* mat = (LogicMatrix*)ptr;
	LogicMatrix* res = mat->getColVector(j);
	return (jlong)res;
	
}


JNIEXPORT jint JNICALL Java_blas_java_LogicMatrix_getNRowNative(JNIEnv *env, jobject jb, jlong ptr){
	LogicMatrix* mat = (LogicMatrix*)ptr;
	return (jint)mat->getNRow();
}


JNIEXPORT jint JNICALL Java_blas_java_LogicMatrix_getNColNative(JNIEnv *env, jobject jb, jlong ptr){
	LogicMatrix* mat = (LogicMatrix*)ptr;
	return (jint)mat->getNCol();
}


JNIEXPORT void JNICALL Java_blas_java_LogicMatrix_clearNative(JNIEnv *env, jobject jb, jlong ptr){
	delete((LogicMatrix*)ptr);
	return ;
}


JNIEXPORT void JNICALL Java_blas_java_LogicMatrix_printNative(JNIEnv *env, jobject jb, jlong ptr){
	LogicMatrix* mat = (LogicMatrix*)ptr;
	mat->print();
}


JNIEXPORT void JNICALL Java_blas_java_LogicMatrix_saveMatrixNative(JNIEnv *env, jobject jb, jstring file, jlong ptr){
	LogicMatrix* mat = (LogicMatrix*)ptr;
	std::string filename = env->GetStringUTFChars(file, 0);
	mat->saveMatrix(filename);
}


JNIEXPORT void JNICALL Java_blas_java_LogicMatrix_readMatrixNative(JNIEnv *env, jobject jb, jstring file, jlong ptr){
	LogicMatrix* mat = (LogicMatrix*)ptr;
	std::string filename = env->GetStringUTFChars(file, 0);
	mat->readMatrix(filename);
}



JNIEXPORT jlong JNICALL Java_blas_java_LogicMatrix_ccLogicMatrixNative(JNIEnv *env, jobject jb, jint n, jint m){
	jint size = n * m;
	LogicMatrix* mat = new LogicMatrix(n, m);
	mat->setValue(true);
	std::cout << "C++ constructor done"<<std::endl;
	return (jlong) mat;

}
