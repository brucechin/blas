#include "../Matrix.h"
#include "blas_java_Matrix.h"
#include<iostream>
#include<stdexcept>
#include<unistd.h>
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



JNIEXPORT void JNICALL Java_blas_java_Matrix_clear(JNIEnv *env, jobject jb){
	

}

JNIEXPORT void JNICALL Java_blas_java_Matrix_getNCol(JNIEnv *env, jobject jb){


}

JNIEXPORT void JNICALL Java_blas_java_Matrix_getNRow(JNIEnv *env, jobject jb){


}
JNIEXPORT void JNICALL Java_blas_java_Matrix_print(JNIEnv *env, jobject jb){


}
/*

JNIEXPORT void JNICALL Java_blas_java_Matrix_saveMatrix(JNIEnv *env, jobject jb, jstring filename){


}


JNIEXPORT void JNICALL Java_blas_java_Matrix_readMatrix(JNIEnv *env, jobject jb, jstring filename){

}
*/

JNIEXPORT jlong JNICALL Java_blas_java_Matrix_ccMatrix(JNIEnv *env, jobject jb, jint n, jint m){
	//long ptr;
	jint size = n * m;
	
	Matrix* mat = new Matrix(n, m);
	mat->setValue(10);
	std::cout << "constructor done"<<std::endl;
	//usleep(3000000);
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