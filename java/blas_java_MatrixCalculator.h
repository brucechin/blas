
#include <jni.h>


#ifndef _Included_blas_java_MatrixCalculator
#define _Included_blas_java_MatrixCalculator
#ifdef __cplusplus
extern "C" {
#endif




JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_addNative__JJ
  (JNIEnv *, jclass, jlong, jlong);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_subNative__JJ
  (JNIEnv *, jclass, jlong, jlong);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_divNative__JJ
  (JNIEnv *, jclass, jlong, jlong);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_mulNative__JJ
  (JNIEnv *, jclass, jlong, jlong);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_mulNative__DJ
  (JNIEnv *, jclass, jdouble, jlong);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_matrixMulNative__JJ
  (JNIEnv *, jclass, jlong, jlong);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_maxNative__JJ
  (JNIEnv *, jclass, jlong, jlong);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_minNative__JJ
  (JNIEnv *, jclass, jlong, jlong);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_biggerNative__JJ
  (JNIEnv *, jclass, jlong, jlong);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_smallerNative__JJ
  (JNIEnv *, jclass, jlong, jlong);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_biggerNative__JD
  (JNIEnv *, jclass, jlong, jdouble);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_smallerNative__JD
  (JNIEnv *, jclass, jlong, jdouble);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_equalNative__JJ
  (JNIEnv *, jclass, jlong, jlong);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_equalNative__JD
  (JNIEnv *, jclass, jlong, jdouble);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_orNative__JJ
  (JNIEnv *, jclass, jlong, jlong);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_andNative__JJ
  (JNIEnv *, jclass, jlong, jlong);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_notNative__J
  (JNIEnv *, jclass, jlong);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_conditionNative__JJJ
  (JNIEnv *, jclass, jlong, jlong, jlong);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_betweenNative
  (JNIEnv *, jclass, jlong, jdouble, jdouble);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_betweenValueNative
  (JNIEnv *, jclass, jlong, jdouble, jdouble);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_rankNative__J
  (JNIEnv *, jclass, jlong);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_roundNative__J
  (JNIEnv *, jclass, jlong);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_floorNative__J
  (JNIEnv *, jclass, jlong);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_absNative__J
  (JNIEnv *, jclass, jlong);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_minusNative__J
  (JNIEnv *, jclass, jlong);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_sqrtNative__J
  (JNIEnv *, jclass, jlong);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_logNative__J
  (JNIEnv *, jclass, jlong);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_expNative__J
  (JNIEnv *, jclass, jlong);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_signNative__J
  (JNIEnv *, jclass, jlong);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_inverseNative__J
  (JNIEnv *, jclass, jlong);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_signedpowNative__JD
  (JNIEnv *, jclass, jlong, jdouble);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_shiftNative__JI
  (JNIEnv *, jclass, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_delayNative__JI
  (JNIEnv *, jclass, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_deltaNative__JI
  (JNIEnv *, jclass, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_ratioNative__JI
  (JNIEnv *, jclass, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_sumNative__JI
  (JNIEnv *, jclass, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_productNative__JI
  (JNIEnv *, jclass, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsMaxNative__JI
  (JNIEnv *, jclass, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsMinNative__JI
  (JNIEnv *, jclass, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsArgmaxNative__JI
  (JNIEnv *, jclass, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsArgminNative__JI
  (JNIEnv *, jclass, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsRankNative__JI
  (JNIEnv *, jclass, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsMeanNative__JI
  (JNIEnv *, jclass, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsStdNative__JI
  (JNIEnv *, jclass, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsSkewnessNative__JI
  (JNIEnv *, jclass, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsKurtosisNative__JI
  (JNIEnv *, jclass, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsCovNative__JJI
  (JNIEnv *, jclass, jlong, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsCorrNative__JJI
  (JNIEnv *, jclass, jlong, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsCountNaNNative__JI
  (JNIEnv *, jclass, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsCountTrueNative__JI
  (JNIEnv *, jclass, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsCountConsecutiveTrueNative__JI
  (JNIEnv *, jclass, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_decayLinearNative__JI
  (JNIEnv *, jclass, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_decayExponentialNative__JI
  (JNIEnv *, jclass, jlong, jint);


JNIEXPORT void JNICALL Java_blas_java_MatrixCalculator_smoothByDecayLinearNative
  (JNIEnv *, jclass, jlong, jint);


JNIEXPORT void JNICALL Java_blas_java_MatrixCalculator_inputNaNNative
  (JNIEnv *, jclass, jlong, jdouble);


JNIEXPORT void JNICALL Java_blas_java_MatrixCalculator_activateNative
  (JNIEnv *, jclass, jlong, jdouble);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_normalizeNative__JDDD
  (JNIEnv *, jclass, jlong, jdouble, jdouble, jdouble);


JNIEXPORT void JNICALL Java_blas_java_MatrixCalculator_normalizeBySpecNative
  (JNIEnv *, jclass, jlong, jdouble, jdouble, jdouble);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_neutralizeNative__J
  (JNIEnv *, jclass, jlong);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_meanNative
  (JNIEnv *, jclass, jlong);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_unifyNative__J
  (JNIEnv *, jclass, jlong);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_unifyByL2Native
  (JNIEnv *, jclass, jlong);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_evalValidPctNative
  (JNIEnv *, jclass, jlong);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_evalAbsSumNative
  (JNIEnv *, jclass, jlong);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_evalMeanNative
  (JNIEnv *, jclass, jlong);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_evalVarianceNative
  (JNIEnv *, jclass, jlong);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_evalInnerProductionNative
  (JNIEnv *, jclass, jlong, jlong);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_evalCovarianceNative
  (JNIEnv *, jclass, jlong, jlong);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_evalCorrelationNative
  (JNIEnv *, jclass, jlong, jlong);


JNIEXPORT jdouble JNICALL Java_blas_java_MatrixCalculator_DetNative
  (JNIEnv *, jclass, jlong, jint);


JNIEXPORT jdouble JNICALL Java_blas_java_MatrixCalculator_InverseNative
  (JNIEnv *, jclass, jlong, jint, jlong);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_invNative
  (JNIEnv *, jclass, jlong);


JNIEXPORT jdouble JNICALL Java_blas_java_MatrixCalculator_treatNative
  (JNIEnv *, jclass, jlong);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_diagNative
  (JNIEnv *, jclass, jlong);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_inverseDiagNative
  (JNIEnv *, jclass, jlong);


JNIEXPORT jdouble JNICALL Java_blas_java_MatrixCalculator_evalInnerProductionByLongVectorNative
  (JNIEnv *, jclass, jlong, jlong);


JNIEXPORT jdouble JNICALL Java_blas_java_MatrixCalculator_evalCorrelationByLongVectorNative
  (JNIEnv *, jclass, jlong, jlong);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_evalBetaNative
  (JNIEnv *, jclass, jlong, jlong);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_cumSumNative
  (JNIEnv *, jclass, jlong);


JNIEXPORT jdouble JNICALL Java_blas_java_MatrixCalculator_evalBetaByLongVectorNative
  (JNIEnv *, jclass, jlong, jlong);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_addNative__JJI
  (JNIEnv *, jclass, jlong, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_subNative__JJI
  (JNIEnv *, jclass, jlong, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_divNative__JJI
  (JNIEnv *, jclass, jlong, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_mulNative__JJI
  (JNIEnv *, jclass, jlong, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_mulNative__DJI
  (JNIEnv *, jclass, jdouble, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_matrixMulNative__JJI
  (JNIEnv *, jclass, jlong, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_maxNative__JJI
  (JNIEnv *, jclass, jlong, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_minNative__JJI
  (JNIEnv *, jclass, jlong, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_biggerNative__JJI
  (JNIEnv *, jclass, jlong, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_smallerNative__JJI
  (JNIEnv *, jclass, jlong, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_biggerNative__JDI
  (JNIEnv *, jclass, jlong, jdouble, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_smallerNative__JDI
  (JNIEnv *, jclass, jlong, jdouble, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_equalNative__JJI
  (JNIEnv *, jclass, jlong, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_equalNative__JDI
  (JNIEnv *, jclass, jlong, jdouble, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_orNative__JJI
  (JNIEnv *, jclass, jlong, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_andNative__JJI
  (JNIEnv *, jclass, jlong, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_notNative__JI
  (JNIEnv *, jclass, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_conditionNative__JJJI
  (JNIEnv *, jclass, jlong, jlong, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_rankNative__JI
  (JNIEnv *, jclass, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_roundNative__JI
  (JNIEnv *, jclass, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_floorNative__JI
  (JNIEnv *, jclass, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_absNative__JI
  (JNIEnv *, jclass, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_minusNative__JI
  (JNIEnv *, jclass, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_sqrtNative__JI
  (JNIEnv *, jclass, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_logNative__JI
  (JNIEnv *, jclass, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_expNative__JI
  (JNIEnv *, jclass, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_signNative__JI
  (JNIEnv *, jclass, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_inverseNative__JI
  (JNIEnv *, jclass, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_signedpowNative__JDI
  (JNIEnv *, jclass, jlong, jdouble, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_shiftNative__JII
  (JNIEnv *, jclass, jlong, jint, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_delayNative__JII
  (JNIEnv *, jclass, jlong, jint, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_deltaNative__JII
  (JNIEnv *, jclass, jlong, jint, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_ratioNative__JII
  (JNIEnv *, jclass, jlong, jint, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_sumNative__JII
  (JNIEnv *, jclass, jlong, jint, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_productNative__JII
  (JNIEnv *, jclass, jlong, jint, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsMaxNative__JII
  (JNIEnv *, jclass, jlong, jint, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsMinNative__JII
  (JNIEnv *, jclass, jlong, jint, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsArgmaxNative__JII
  (JNIEnv *, jclass, jlong, jint, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsArgminNative__JII
  (JNIEnv *, jclass, jlong, jint, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsRankNative__JII
  (JNIEnv *, jclass, jlong, jint, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsMeanNative__JII
  (JNIEnv *, jclass, jlong, jint, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsStdNative__JII
  (JNIEnv *, jclass, jlong, jint, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsSkewnessNative__JII
  (JNIEnv *, jclass, jlong, jint, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsKurtosisNative__JII
  (JNIEnv *, jclass, jlong, jint, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsCovNative__JJII
  (JNIEnv *, jclass, jlong, jlong, jint, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsCorrNative__JJII
  (JNIEnv *, jclass, jlong, jlong, jint, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsCountNaNNative__JII
  (JNIEnv *, jclass, jlong, jint, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsCountTrueNative__JII
  (JNIEnv *, jclass, jlong, jint, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsCountConsecutiveTrueNative__JII
  (JNIEnv *, jclass, jlong, jint, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_decayLinearNative__JII
  (JNIEnv *, jclass, jlong, jint, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_decayExponentialNative__JII
  (JNIEnv *, jclass, jlong, jint, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_normalizeNative__JDDDI
  (JNIEnv *, jclass, jlong, jdouble, jdouble, jdouble, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_neutralizeNative__JI
  (JNIEnv *, jclass, jlong, jint);


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_unifyNative__JI
  (JNIEnv *, jclass, jlong, jint);

#ifdef __cplusplus
}
#endif
#endif
