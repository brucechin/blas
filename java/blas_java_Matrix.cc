#include "../Matrix.h"
#include "blas_java_Matrix.h"
#include<iostream>
#include<stdexcept>
#include<unistd.h>
#include<cstring>
static JavaVM *cached_jvm = 0;

/*
JNIEXPORT jint JNICALL JNI_OnLoad(JavaVM* jvm, void* reserved){
	cached_jvm = jvm;
	return JNI_VERSION_1_2;
}

static JNIEnv* JNU_GetEnv(){
	JNIEnv* env;
	jint rc = cached_jvm->GetEnv((void **)&env, JNI_VERSION_1_2);
	return env;
}

*/

JNIEXPORT void JNICALL Java_blas_java_Matrix_setElementNative(JNIEnv *env, jobject jb, jint i, jint j, jdouble val, jlong ptr){
	Matrix* mat = (Matrix*)ptr;
	mat->setElement(i, j, val);
}


JNIEXPORT jdouble JNICALL Java_blas_java_Matrix_getElementNative(JNIEnv *env, jobject jb, jint i, jint j, jlong ptr){
	Matrix* mat = (Matrix*)ptr;
	jdouble res = mat->getElement(i, j);
	return res;
}


JNIEXPORT jlong JNICALL Java_blas_java_Matrix_getRowVectorNative(JNIEnv *env, jobject jb, jint i, jlong ptr){
	Matrix* mat = (Matrix*)ptr;
	Matrix* res = mat->getRowVector(i);
	return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_Matrix_getColVectorNative(JNIEnv *env, jobject jb, jint j, jlong ptr){
	Matrix* mat = (Matrix*)ptr;
	Matrix* res = mat->getColVector(j);
	return (jlong)res;
	
}


JNIEXPORT jint JNICALL Java_blas_java_Matrix_getNRowNative(JNIEnv *env, jobject jb, jlong ptr){
	Matrix* mat = (Matrix*)ptr;
	return (jint)mat->getNRow();
}


JNIEXPORT jint JNICALL Java_blas_java_Matrix_getNColNative(JNIEnv *env, jobject jb, jlong ptr){
	Matrix* mat = (Matrix*)ptr;
	return (jint)mat->getNCol();
}


JNIEXPORT void JNICALL Java_blas_java_Matrix_clearNative(JNIEnv *env, jobject jb, jlong ptr){
	delete((Matrix*)ptr);
	return ;
}


JNIEXPORT void JNICALL Java_blas_java_Matrix_printNative(JNIEnv *env, jobject jb, jlong ptr){
	Matrix* mat = (Matrix*)ptr;
	mat->print();
}


JNIEXPORT void JNICALL Java_blas_java_Matrix_saveMatrixNative(JNIEnv *env, jobject jb, jstring file, jlong ptr){
	Matrix* mat = (Matrix*)ptr;
	std::string filename = env->GetStringUTFChars(file, 0);
	mat->saveMatrix(filename);
}


JNIEXPORT void JNICALL Java_blas_java_Matrix_readMatrixNative(JNIEnv *env, jobject jb, jstring file, jlong ptr){
	Matrix* mat = (Matrix*)ptr;
	std::string filename = env->GetStringUTFChars(file, 0);
	mat->readMatrix(filename);
}



JNIEXPORT jlong JNICALL Java_blas_java_Matrix_ccMatrixNative(JNIEnv *env, jobject jb, jint n, jint m){
	jint size = n * m;
	Matrix* mat = new Matrix(n, m);
	mat->setValue(10);
	std::cout << "C++ constructor done"<<std::endl;
	return (jlong) mat;


	/*
	std::cout <<"start allocate"<<std::endl;
	usleep(10000000);
	jdoubleArray buffer = env->NewDoubleArray(size);
	std::cout <<"end allocate";
	//memset(buffer,0.0, size*sizeof(double));
	//double* res = new double[size];
	Matrix* mat = new Matrix(0, 0);
	//for(int i = 0; i < size; i++) res[i] = 99.99;
	mat->setNRow(n);
	mat->setNCol(m);
	jboolean isCopy = JNI_FALSE;
	std::cout << "start copying"<<std::endl;
	usleep(10000000);
	void* jbuffer = static_cast<void *>(env->GetDoubleArrayElements(buffer, &isCopy));
	std::cout <<"end copying"<<std::endl;
	if(isCopy == JNI_TRUE){
		std::cout<<"copy!!!"<<std::endl;
	}else{
		std::cout<<"not a copy"<<std::endl;
	}
	//mat->value = static_cast<double *>(jbuffer);
	mat->value = (double*)buffer;
	std::cout << "jdoubleArray allocated"<<std::endl;
	//mat->print();
	//mat->setValue(10);
	//env->SetDoubleArrayRegion(buffer, 0, size, res);
	std::cout << "jdoubleArray written"<<std::endl;
	//mat->print();
	//mat->setValue(10);
	//mat->print();
	return (jlong) mat;
*/
}
