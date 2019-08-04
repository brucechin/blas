#include<iostream>
#include<stdexcept>
#include<unistd.h>
#include<cstring>
#include"../Matrix.h"
#include"../LogicMatrix.h"
#include"../MatrixCalculator.h"
#include"blas_java_MatrixCalculator.h"


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_addNative__JJ(JNIEnv * env, jclass obj, jlong p1, jlong p2){
    Matrix* res = MatrixCalculator::add((Matrix*) p1, (Matrix*) p2);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_subNative__JJ(JNIEnv * env, jclass obj, jlong p1, jlong p2){
    Matrix* res = MatrixCalculator::sub((Matrix*) p1, (Matrix*) p2);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_divNative__JJ(JNIEnv * env, jclass obj, jlong p1, jlong p2){
    Matrix* res = MatrixCalculator::div((Matrix*) p1, (Matrix*) p2);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_mulNative__JJ(JNIEnv * env, jclass obj, jlong p1, jlong p2){
    Matrix* res = MatrixCalculator::mul((Matrix*) p1, (Matrix*) p2);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_mulNative__DJ(JNIEnv * env, jclass obj, jdouble val, jlong p2){
    Matrix* res = MatrixCalculator::mul(val, (Matrix*) p2);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_matrixMulNative__JJ(JNIEnv * env, jclass obj, jlong p1, jlong p2){
    Matrix* res = MatrixCalculator::matrixMul((Matrix*) p1, (Matrix*) p2);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_maxNative__JJ(JNIEnv * env, jclass obj, jlong p1, jlong p2){
    Matrix* res = MatrixCalculator::max((Matrix*) p1, (Matrix*) p2);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_minNative__JJ(JNIEnv * env, jclass obj, jlong p1, jlong p2){
    Matrix* res = MatrixCalculator::min((Matrix*) p1, (Matrix*) p2);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_biggerNative__JJ(JNIEnv * env, jclass obj, jlong p1, jlong p2){
    LogicMatrix* res = MatrixCalculator::bigger((Matrix*) p1, (Matrix*) p2);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_smallerNative__JJ(JNIEnv * env, jclass obj, jlong p1, jlong p2){
    LogicMatrix* res = MatrixCalculator::smaller((Matrix*) p1, (Matrix*) p2);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_biggerNative__JD(JNIEnv * env, jclass obj, jlong p1, jdouble val){
    LogicMatrix* res = MatrixCalculator::bigger((Matrix*) p1, val);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_smallerNative__JD(JNIEnv * env, jclass obj, jlong p1, jdouble val){
    LogicMatrix* res = MatrixCalculator::smaller((Matrix*) p1, val);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_equalNative__JJ(JNIEnv * env, jclass obj, jlong p1, jlong p2){
    LogicMatrix* res = MatrixCalculator::equal((Matrix*) p1, (Matrix*) p2);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_equalNative__JD(JNIEnv * env, jclass obj, jlong p1, jdouble val){
    LogicMatrix* res = MatrixCalculator::equal((Matrix*) p1, val);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_orNative__JJ(JNIEnv * env, jclass obj, jlong p1, jlong p2){
    LogicMatrix* res = MatrixCalculator::matOr((LogicMatrix*) p1, (LogicMatrix*) p2);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_andNative__JJ(JNIEnv * env, jclass obj, jlong p1, jlong p2){
    LogicMatrix* res = MatrixCalculator::matAnd((LogicMatrix*) p1, (LogicMatrix*) p2);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_notNative__J(JNIEnv * env, jclass obj, jlong p1){
    LogicMatrix* res = MatrixCalculator::matNot((LogicMatrix*) p1);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_conditionNative__JJJ(JNIEnv * env, jclass obj, jlong p1, jlong p2, jlong p3){
	Matrix* res = MatrixCalculator::condition((LogicMatrix*) p1, (Matrix*) p2, (Matrix*) p3);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_betweenNative(JNIEnv * env, jclass obj, jlong p1, jdouble lowerbound, jdouble upperbound){
    LogicMatrix* res = MatrixCalculator::between((Matrix*) p1, lowerbound, upperbound);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_betweenValueNative(JNIEnv * env, jclass obj, jlong p1, jdouble lowerbound, jdouble upperbound){
    Matrix* res = MatrixCalculator::betweenValue((Matrix*) p1, lowerbound, upperbound);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_rankNative__J(JNIEnv * env, jclass obj, jlong p1){
    Matrix* res = MatrixCalculator::rank((Matrix*) p1);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_roundNative__J(JNIEnv * env, jclass obj, jlong p1){
    Matrix* res = MatrixCalculator::round((Matrix*) p1);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_floorNative__J(JNIEnv * env, jclass obj, jlong p1){
    Matrix* res = MatrixCalculator::floor((Matrix*) p1);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_absNative__J(JNIEnv * env, jclass obj, jlong p1){
    Matrix* res = MatrixCalculator::abs((Matrix*) p1);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_minusNative__J(JNIEnv * env, jclass obj, jlong p1){
    Matrix* res = MatrixCalculator::minus((Matrix*) p1);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_sqrtNative__J(JNIEnv * env, jclass obj, jlong p1){
    Matrix* res = MatrixCalculator::sqrt((Matrix*) p1);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_logNative__J(JNIEnv * env, jclass obj, jlong p1){
    Matrix* res = MatrixCalculator::log((Matrix*) p1);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_expNative__J(JNIEnv * env, jclass obj, jlong p1){
    Matrix* res = MatrixCalculator::exp((Matrix*) p1);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_signNative__J(JNIEnv * env, jclass obj, jlong p1){
    Matrix* res = MatrixCalculator::sign((Matrix*) p1);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_inverseNative__J(JNIEnv * env, jclass obj, jlong p1){
    Matrix* res = MatrixCalculator::inverse((Matrix*) p1);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_signedpowNative__JD(JNIEnv * env, jclass obj, jlong p1, jdouble index){
    Matrix* res = MatrixCalculator::signedpow((Matrix*) p1, index);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_shiftNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n){
    Matrix* res = MatrixCalculator::shift((Matrix*) p1, n);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_delayNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n){
    Matrix* res = MatrixCalculator::delay((Matrix*) p1, n);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_deltaNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n){
    Matrix* res = MatrixCalculator::delta((Matrix*) p1, n);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_ratioNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n){
    Matrix* res = MatrixCalculator::ratio((Matrix*) p1, n);
    return (jlong)res;
}   


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_sumNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n){
    Matrix* res = MatrixCalculator::sum((Matrix*) p1, n);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_productNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n){
    Matrix* res = MatrixCalculator::product((Matrix*) p1, n);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsMaxNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n){
    Matrix* res = MatrixCalculator::tsMax((Matrix*) p1, n);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsMinNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n){
    Matrix* res = MatrixCalculator::tsMin((Matrix*) p1, n);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsArgmaxNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n){
    Matrix* res = MatrixCalculator::tsArgmax((Matrix*) p1, n);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsArgminNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n){
    Matrix* res = MatrixCalculator::tsArgmin((Matrix*) p1, n);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsRankNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n){
    Matrix* res = MatrixCalculator::tsRank((Matrix*) p1, n);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsMeanNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n){
    Matrix* res = MatrixCalculator::tsMean((Matrix*) p1, n);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsStdNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n){
    Matrix* res = MatrixCalculator::tsStd((Matrix*) p1, n);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsSkewnessNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n){
    Matrix* res = MatrixCalculator::tsSkewness((Matrix*) p1, n);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsKurtosisNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n){
    Matrix* res = MatrixCalculator::tsKurtosis((Matrix*) p1, n);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsCovNative__JJI(JNIEnv * env, jclass obj, jlong p1, jlong p2, jint n){
    Matrix* res = MatrixCalculator::tsCov((Matrix*) p1, (Matrix*) p2, n);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsCorrNative__JJI(JNIEnv * env, jclass obj, jlong p1, jlong p2, jint n){
    Matrix* res = MatrixCalculator::tsCorr((Matrix*) p1, (Matrix*) p2, n);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsCountNaNNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n){
    Matrix* res = MatrixCalculator::tsCountNaN((Matrix*) p1, n);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsCountTrueNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n){
    Matrix* res = MatrixCalculator::tsCountTrue((LogicMatrix*) p1, n);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsCountConsecutiveTrueNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n){
    Matrix* res = MatrixCalculator::tsCountConsecutiveTrue((LogicMatrix*) p1, n);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_decayLinearNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n){
    Matrix* res = MatrixCalculator::decayLinear((Matrix*) p1, n);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_decayExponentialNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n){
    Matrix* res = MatrixCalculator::decayExponential((Matrix*) p1, n);
    return (jlong)res;
}


JNIEXPORT void JNICALL Java_blas_java_MatrixCalculator_smoothByDecayLinearNative(JNIEnv * env, jclass obj, jlong p1, jint n){
    MatrixCalculator::smoothByDecayLinear((Matrix*) p1, n);
}


JNIEXPORT void JNICALL Java_blas_java_MatrixCalculator_inputNaNNative(JNIEnv * env, jclass obj, jlong p1, jdouble val){
    MatrixCalculator::inputNaN((Matrix*) p1, val);
}


JNIEXPORT void JNICALL Java_blas_java_MatrixCalculator_activateNative(JNIEnv * env, jclass obj, jlong p1, jdouble val){
    MatrixCalculator::activate((Matrix*) p1, val);
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_normalizeNative__JDDD(JNIEnv * env, jclass obj, jlong p1, jdouble scale, jdouble mean, jdouble bound){
    Matrix* res = MatrixCalculator::normalize((Matrix*) p1, scale, mean, bound);
    return (jlong)res;
}


JNIEXPORT void JNICALL Java_blas_java_MatrixCalculator_normalizeBySpecNative(JNIEnv * env, jclass obj, jlong p1, jdouble scale, jdouble mean, jdouble bound){
    MatrixCalculator::normalizeBySpec((Matrix*) p1, scale, mean, bound);
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_neutralizeNative__J(JNIEnv * env, jclass obj, jlong p1){
    Matrix* res = MatrixCalculator::neutralize((Matrix*) p1);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_meanNative(JNIEnv * env, jclass obj, jlong p1){
    Matrix* res = MatrixCalculator::mean((Matrix*) p1);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_unifyNative__J(JNIEnv * env, jclass obj, jlong p1){
    Matrix* res = MatrixCalculator::unify((Matrix*) p1);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_unifyByL2Native(JNIEnv * env, jclass obj, jlong p1){
    Matrix* res = MatrixCalculator::unifyByL2((Matrix*) p1);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_evalValidPctNative(JNIEnv * env, jclass obj, jlong p1){
    Matrix* res = MatrixCalculator::evalValidPct((Matrix*) p1);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_evalAbsSumNative(JNIEnv * env, jclass obj, jlong p1){
    Matrix* res = MatrixCalculator::evalAbsSum((Matrix*) p1);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_evalMeanNative(JNIEnv * env, jclass obj, jlong p1){
    Matrix* res = MatrixCalculator::evalMean((Matrix*) p1);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_evalVarianceNative(JNIEnv * env, jclass obj, jlong p1){
    Matrix* res = MatrixCalculator::evalVariance((Matrix*) p1);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_evalInnerProductionNative(JNIEnv * env, jclass obj, jlong p1, jlong p2){
    Matrix* res = MatrixCalculator::evalInnerProduction((Matrix*) p1, (Matrix*) p2);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_evalCovarianceNative(JNIEnv * env, jclass obj, jlong p1, jlong p2){
    Matrix* res = MatrixCalculator::evalCovariance((Matrix*) p1, (Matrix*) p2);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_evalCorrelationNative(JNIEnv * env, jclass obj, jlong p1, jlong p2){
    Matrix* res = MatrixCalculator::evalCorrelation((Matrix*) p1, (Matrix*) p2);
    return (jlong)res;
}


JNIEXPORT jdouble JNICALL Java_blas_java_MatrixCalculator_DetNative(JNIEnv * env, jclass obj, jlong p1, jint N){
    jdouble res = MatrixCalculator::Det((Matrix*) p1, N);
    return res;
}


JNIEXPORT jdouble JNICALL Java_blas_java_MatrixCalculator_InverseNative(JNIEnv * env, jclass obj, jlong p1, jint N, jlong p3){
    jdouble res = MatrixCalculator::Inverse((Matrix*) p1, N, (Matrix*) p3);
    return res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_invNative(JNIEnv * env, jclass obj, jlong p1){
    Matrix* res = MatrixCalculator::inv((Matrix*) p1);
    return (jlong)res;
}


JNIEXPORT jdouble JNICALL Java_blas_java_MatrixCalculator_treatNative(JNIEnv * env, jclass obj, jlong p1){
    jdouble res = MatrixCalculator::treat((Matrix*) p1);
    return res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_diagNative(JNIEnv * env, jclass obj, jlong p1){
    Matrix* res = MatrixCalculator::diag((Matrix*) p1);
    return (jlong)res;
}   


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_inverseDiagNative(JNIEnv * env, jclass obj, jlong p1){
    Matrix* res = MatrixCalculator::inverseDiag((Matrix*) p1);
    return (jlong)res;
}


JNIEXPORT jdouble JNICALL Java_blas_java_MatrixCalculator_evalInnerProductionByLongVectorNative(JNIEnv * env, jclass obj, jlong p1, jlong p2){
    double res = MatrixCalculator::evalInnerProductionByLongVector((Matrix*) p1, (Matrix*) p2);
    return res;
}


JNIEXPORT jdouble JNICALL Java_blas_java_MatrixCalculator_evalCorrelationByLongVectorNative(JNIEnv * env, jclass obj, jlong p1, jlong p2){
    double res = MatrixCalculator::evalCorrelationByLongVector((Matrix*) p1, (Matrix*) p2);
    return res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_evalBetaNative(JNIEnv * env, jclass obj, jlong p1, jlong p2){
    Matrix* res = MatrixCalculator::evalBeta((Matrix*) p1, (Matrix*) p2);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_cumSumNative(JNIEnv * env, jclass obj, jlong p1){
    Matrix* res = MatrixCalculator::cumSum((Matrix*) p1);
    return (jlong)res;
}


JNIEXPORT jdouble JNICALL Java_blas_java_MatrixCalculator_evalBetaByLongVectorNative(JNIEnv * env, jclass obj, jlong p1, jlong p2){
    jdouble res = MatrixCalculator::evalBetaByLongVector((Matrix*) p1, (Matrix*) p2);
    return res;
}



JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_addNative__JJ(JNIEnv * env, jclass obj, jlong p1, jlong p2, jint num){
    Matrix* res = MatrixCalculator::add((Matrix*) p1, (Matrix*) p2, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_subNative__JJ(JNIEnv * env, jclass obj, jlong p1, jlong p2, jint num){
    Matrix* res = MatrixCalculator::sub((Matrix*) p1, (Matrix*) p2, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_divNative__JJ(JNIEnv * env, jclass obj, jlong p1, jlong p2, jint num){
    Matrix* res = MatrixCalculator::div((Matrix*) p1, (Matrix*) p2, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_mulNative__JJ(JNIEnv * env, jclass obj, jlong p1, jlong p2, jint num){
    Matrix* res = MatrixCalculator::mul((Matrix*) p1, (Matrix*) p2, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_mulNative__DJ(JNIEnv * env, jclass obj, jdouble val, jlong p2, jint num){
    Matrix* res = MatrixCalculator::mul(val, (Matrix*) p2, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_maxNative__JJ(JNIEnv * env, jclass obj, jlong p1, jlong p2, jint num){
    Matrix* res = MatrixCalculator::max((Matrix*) p1, (Matrix*) p2, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_minNative__JJ(JNIEnv * env, jclass obj, jlong p1, jlong p2, jint num){
    Matrix* res = MatrixCalculator::min((Matrix*) p1, (Matrix*) p2, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_biggerNative__JJ(JNIEnv * env, jclass obj, jlong p1, jlong p2, jint num){
    LogicMatrix* res = MatrixCalculator::bigger((Matrix*) p1, (Matrix*) p2, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_smallerNative__JJ(JNIEnv * env, jclass obj, jlong p1, jlong p2, jint num){
    LogicMatrix* res = MatrixCalculator::smaller((Matrix*) p1, (Matrix*) p2, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_biggerNative__JD(JNIEnv * env, jclass obj, jlong p1, jdouble val, jint num){
    LogicMatrix* res = MatrixCalculator::bigger((Matrix*) p1, val, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_smallerNative__JD(JNIEnv * env, jclass obj, jlong p1, jdouble val, jint num){
    LogicMatrix* res = MatrixCalculator::smaller((Matrix*) p1, val, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_equalNative__JJ(JNIEnv * env, jclass obj, jlong p1, jlong p2, jint num){
    LogicMatrix* res = MatrixCalculator::equal((Matrix*) p1, (Matrix*) p2, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_equalNative__JD(JNIEnv * env, jclass obj, jlong p1, jdouble val, jint num){
    LogicMatrix* res = MatrixCalculator::equal((Matrix*) p1, val, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_orNative__JJ(JNIEnv * env, jclass obj, jlong p1, jlong p2, jint num){
    LogicMatrix* res = MatrixCalculator::matOr((LogicMatrix*) p1, (LogicMatrix*) p2, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_andNative__JJ(JNIEnv * env, jclass obj, jlong p1, jlong p2, jint num){
    LogicMatrix* res = MatrixCalculator::matAnd((LogicMatrix*) p1, (LogicMatrix*) p2, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_notNative__J(JNIEnv * env, jclass obj, jlong p1, jint num){
    LogicMatrix* res = MatrixCalculator::matNot((LogicMatrix*) p1, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_conditionNative__JJJ(JNIEnv * env, jclass obj, jlong p1, jlong p2, jlong p3, jint num){
    Matrix* res = MatrixCalculator::condition((LogicMatrix*) p1, (Matrix*) p2, (Matrix*) p3, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_rankNative__J(JNIEnv * env, jclass obj, jlong p1, jint num){
    Matrix* res = MatrixCalculator::rank((Matrix*) p1, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_roundNative__J(JNIEnv * env, jclass obj, jlong p1, jint num){
    Matrix* res = MatrixCalculator::round((Matrix*) p1, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_floorNative__J(JNIEnv * env, jclass obj, jlong p1, jint num){
    Matrix* res = MatrixCalculator::floor((Matrix*) p1, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_absNative__J(JNIEnv * env, jclass obj, jlong p1, jint num){
    Matrix* res = MatrixCalculator::abs((Matrix*) p1, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_minusNative__J(JNIEnv * env, jclass obj, jlong p1, jint num){
    Matrix* res = MatrixCalculator::minus((Matrix*) p1, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_sqrtNative__J(JNIEnv * env, jclass obj, jlong p1, jint num){
    Matrix* res = MatrixCalculator::sqrt((Matrix*) p1, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_logNative__J(JNIEnv * env, jclass obj, jlong p1, jint num){
    Matrix* res = MatrixCalculator::log((Matrix*) p1, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_expNative__J(JNIEnv * env, jclass obj, jlong p1, jint num){
    Matrix* res = MatrixCalculator::exp((Matrix*) p1, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_signNative__J(JNIEnv * env, jclass obj, jlong p1, jint num){
    Matrix* res = MatrixCalculator::sign((Matrix*) p1, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_inverseNative__J(JNIEnv * env, jclass obj, jlong p1, jint num){
    Matrix* res = MatrixCalculator::inverse((Matrix*) p1, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_signedpowNative__JD(JNIEnv * env, jclass obj, jlong p1, jdouble index, jint num){
    Matrix* res = MatrixCalculator::signedpow((Matrix*) p1, index, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_shiftNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n, jint num){
    Matrix* res = MatrixCalculator::shift((Matrix*) p1, n, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_delayNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n, jint num){
    Matrix* res = MatrixCalculator::delay((Matrix*) p1, n, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_deltaNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n, jint num){
    Matrix* res = MatrixCalculator::delta((Matrix*) p1, n, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_ratioNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n, jint num){
    Matrix* res = MatrixCalculator::ratio((Matrix*) p1, n, num);
    return (jlong)res;
}   


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_sumNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n, jint num){
    Matrix* res = MatrixCalculator::sum((Matrix*) p1, n, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_productNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n, jint num){
    Matrix* res = MatrixCalculator::product((Matrix*) p1, n, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsMaxNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n, jint num){
    Matrix* res = MatrixCalculator::tsMax((Matrix*) p1, n, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsMinNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n, jint num){
    Matrix* res = MatrixCalculator::tsMin((Matrix*) p1, n, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsArgmaxNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n, jint num){
    Matrix* res = MatrixCalculator::tsArgmax((Matrix*) p1, n, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsArgminNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n, jint num){
    Matrix* res = MatrixCalculator::tsArgmin((Matrix*) p1, n, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsRankNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n, jint num){
    Matrix* res = MatrixCalculator::tsRank((Matrix*) p1, n, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsMeanNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n, jint num){
    Matrix* res = MatrixCalculator::tsMean((Matrix*) p1, n, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsStdNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n, jint num){
    Matrix* res = MatrixCalculator::tsStd((Matrix*) p1, n, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsSkewnessNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n, jint num){
    Matrix* res = MatrixCalculator::tsSkewness((Matrix*) p1, n, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsKurtosisNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n, jint num){
    Matrix* res = MatrixCalculator::tsKurtosis((Matrix*) p1, n, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsCovNative__JJI(JNIEnv * env, jclass obj, jlong p1, jlong p2, jint n, jint num){
    Matrix* res = MatrixCalculator::tsCov((Matrix*) p1, (Matrix*) p2, n, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsCorrNative__JJI(JNIEnv * env, jclass obj, jlong p1, jlong p2, jint n, jint num){
    Matrix* res = MatrixCalculator::tsCorr((Matrix*) p1, (Matrix*) p2, n, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsCountNaNNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n, jint num){
    Matrix* res = MatrixCalculator::tsCountNaN((Matrix*) p1, n, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsCountTrueNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n, jint num){
    Matrix* res = MatrixCalculator::tsCountTrue((LogicMatrix*) p1, n, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsCountConsecutiveTrueNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n, jint num){
    Matrix* res = MatrixCalculator::tsCountConsecutiveTrue((LogicMatrix*) p1, n, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_decayLinearNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n, jint num){
    Matrix* res = MatrixCalculator::decayLinear((Matrix*) p1, n, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_decayExponentialNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n, jint num){
    Matrix* res = MatrixCalculator::decayExponential((Matrix*) p1, n, num);
    return (jlong)res;
}



JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_normalizeNative__JDDD(JNIEnv * env, jclass obj, jlong p1, jdouble scale, jdouble mean, jdouble bound, jint num){
    Matrix* res = MatrixCalculator::normalize((Matrix*) p1, scale, mean, bound, num);
    return (jlong)res;
}



JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_neutralizeNative__J(JNIEnv * env, jclass obj, jlong p1, jint num){
    Matrix* res = MatrixCalculator::neutralize((Matrix*) p1, num);
    return (jlong)res;
}


JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_unifyNative__J(JNIEnv * env, jclass obj, jlong p1, jint num){
    Matrix* res = MatrixCalculator::unify((Matrix*) p1, num);
    return (jlong)res;
}


































// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_addNative__JJI(JNIEnv * env, jclass obj, jlong, jlong p1, jint n){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_subNative__JJI(JNIEnv * env, jclass obj, jlong, jlong p1, jint n){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_divNative__JJI(JNIEnv * env, jclass obj, jlong, jlong p1, jint n){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_mulNative__JJI(JNIEnv * env, jclass obj, jlong, jlong p1, jint n){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_mulNative__DJI(JNIEnv * env, jclass obj, jdouble, jlong p1, jint n){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_matrixMulNative__JJI(JNIEnv * env, jclass obj, jlong, jlong p1, jint n){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_maxNative__JJI(JNIEnv * env, jclass obj, jlong, jlong p1, jint n){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_minNative__JJI(JNIEnv * env, jclass obj, jlong, jlong p1, jint n){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_biggerNative__JJI(JNIEnv * env, jclass obj, jlong, jlong p1, jint n){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_smallerNative__JJI(JNIEnv * env, jclass obj, jlong, jlong p1, jint n){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_biggerNative__JDI(JNIEnv * env, jclass obj, jlong, jdouble, jint){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_smallerNative__JDI(JNIEnv * env, jclass obj, jlong, jdouble, jint){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_equalNative__JJI(JNIEnv * env, jclass obj, jlong, jlong p1, jint n){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_equalNative__JDI(JNIEnv * env, jclass obj, jlong, jdouble, jint){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_orNative__JJI(JNIEnv * env, jclass obj, jlong, jlong p1, jint n){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_andNative__JJI(JNIEnv * env, jclass obj, jlong, jlong p1, jint n){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_notNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_conditionNative__JJJI(JNIEnv * env, jclass obj, jlong, jlong, jlong p1, jint n){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_rankNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_roundNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_floorNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_absNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_minusNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_sqrtNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_logNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_expNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_signNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_inverseNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_signedpowNative__JDI(JNIEnv * env, jclass obj, jlong, jdouble, jint){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_shiftNative__JII(JNIEnv * env, jclass obj, jlong, jint, jint){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_delayNative__JII(JNIEnv * env, jclass obj, jlong, jint, jint){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_deltaNative__JII(JNIEnv * env, jclass obj, jlong, jint, jint){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_ratioNative__JII(JNIEnv * env, jclass obj, jlong, jint, jint){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_sumNative__JII(JNIEnv * env, jclass obj, jlong, jint, jint){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_productNative__JII(JNIEnv * env, jclass obj, jlong, jint, jint){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsMaxNative__JII(JNIEnv * env, jclass obj, jlong, jint, jint){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsMinNative__JII(JNIEnv * env, jclass obj, jlong, jint, jint){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsArgmaxNative__JII(JNIEnv * env, jclass obj, jlong, jint, jint){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsArgminNative__JII(JNIEnv * env, jclass obj, jlong, jint, jint){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsRankNative__JII(JNIEnv * env, jclass obj, jlong, jint, jint){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsMeanNative__JII(JNIEnv * env, jclass obj, jlong, jint, jint){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsStdNative__JII(JNIEnv * env, jclass obj, jlong, jint, jint){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsSkewnessNative__JII(JNIEnv * env, jclass obj, jlong, jint, jint){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsKurtosisNative__JII(JNIEnv * env, jclass obj, jlong, jint, jint){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsCovNative__JJII(JNIEnv * env, jclass obj, jlong, jlong, jint, jint){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsCorrNative__JJII(JNIEnv * env, jclass obj, jlong, jlong, jint, jint){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsCountNaNNative__JII(JNIEnv * env, jclass obj, jlong, jint, jint){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsCountTrueNative__JII(JNIEnv * env, jclass obj, jlong, jint, jint){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_tsCountConsecutiveTrueNative__JII(JNIEnv * env, jclass obj, jlong, jint, jint){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_decayLinearNative__JII(JNIEnv * env, jclass obj, jlong, jint, jint){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_decayExponentialNative__JII(JNIEnv * env, jclass obj, jlong, jint, jint){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_normalizeNative__JDDDI(JNIEnv * env, jclass obj, jlong, jdouble, jdouble, jdouble, jint){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_neutralizeNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n){

// }


// JNIEXPORT jlong JNICALL Java_blas_java_MatrixCalculator_unifyNative__JI(JNIEnv * env, jclass obj, jlong p1, jint n){

// }
