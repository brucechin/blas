#include<iostream>
#include"blas_java_HelloNative.h"
#include"Test.h"

JNIEXPORT jlong JNICALL Java_blas_java_HelloNative_createCNative(JNIEnv *env, jclass jc){

	Test* res = new Test();
	return (jlong) res;

}

JNIEXPORT jlong JNICALL Java_blas_java_HelloNative_freeCNative(JNIEnv *env, jclass jc, jlong ptr){
	delete((Test*) ptr);
	return 0;
}

JNIEXPORT jlong JNICALL Java_blas_java_HelloNative_sayHello(JNIEnv *env, jclass jc) //难点二

{
	    printf("Hello Native\n");
		return 0;
}
