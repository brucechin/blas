/* DO NOT EDIT THIS FILE - it is machine generated */
:q
#include <jni.h>
/* Header for class blas_java_HelloNative */

#ifndef _Included_blas_java_HelloNative
#define _Included_blas_java_HelloNative
#ifdef __cplusplus
extern "C" {
#endif
/*
 * Class:     blas_java_HelloNative
 * Method:    createCNative
 * Signature: ()J
 */
JNIEXPORT jlong JNICALL Java_blas_java_HelloNative_createCNative
  (JNIEnv *, jclass);

/*
 * Class:     blas_java_HelloNative
 * Method:    freeCNative
 * Signature: (J)J
 */
JNIEXPORT jlong JNICALL Java_blas_java_HelloNative_freeCNative
  (JNIEnv *, jclass, jlong);

/*
 * Class:     blas_java_HelloNative
 * Method:    sayHello
 * Signature: ()J
 */
JNIEXPORT jlong JNICALL Java_blas_java_HelloNative_sayHello
  (JNIEnv *, jclass);

#ifdef __cplusplus
}
#endif
#endif
