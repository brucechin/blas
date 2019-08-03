package blas.java;
public class MatrixCalculator{
	
	static{
		System.loadLibrary("MatrixCalculator");
	}

	long nativeMatrixCalculator;
	
	public MatrixCalculator(){
		nativeMatrixCalculator = createMatrixCalculator();
	}



	public static Matrix add(Matrix mat1, Matrix mat2){
		Matrix res = new Matrix();
		res.setPtr(addNative(mat1.getPtr(), mat2.getPtr()));
		return res;
	}
	public static Matrix sub(Matrix mat1, Matrix mat2){
		Matrix res = new Matrix();
		res.setPtr(subNative(mat1.getPtr(), mat2.getPtr()));
		return res;
	}
	public static Matrix div(Matrix mat1, Matrix mat2){
		Matrix res = new Matrix();
		res.setPtr(divNative(mat1.getPtr(), mat2.getPtr()));
		return res;
	}
	public static Matrix mul(Matrix mat1, Matrix mat2){
		Matrix res = new Matrix();
		res.setPtr(mulNative(mat1.getPtr(), mat2.getPtr()));
		return res;
	}
	public static Matrix mul(double val, Matrix mat){
		Matrix res = new Matrix();
		res.setPtr(addNative(val, mat.getPtr()));
		return res;
	}
	public static Matrix matrixMul(Matrix mat1, Matrix mat2){
		Matrix res = new Matrix();
		res.setPtr(matrixMulNative(mat1.getPtr(), mat2.getPtr()));
		return res;
	}

	public static Matrix max(Matrix mat1, Matrix mat2){
		Matrix res = new Matrix();
		res.setPtr(maxNative(mat1.getPtr(), mat2.getPtr()));
		return res;
	}
	public static Matrix min(Matrix mat1, Matrix mat2){
		Matrix res = new Matrix();
		res.setPtr(minNative(mat1.getPtr(), mat2.getPtr()));
		return res;
	}
	public static LogicMatrix bigger(Matrix mat1, Matrix mat2){
		LogicMatrix res = new LogicMatrix();
		res.setPtr(biggerNative(mat1.getPtr(), mat2.getPtr()));
		return res;
	}
	public static LogicMatrix smaller(Matrix mat1, Matrix mat2){
		LogicMatrix res = new LogicMatrix();
		res.setPtr(smallerNative(mat1.getPtr(), mat2.getPtr()));
		return res;
	}
	public static LogicMatrix bigger(Matrix mat1, double val){
		LogicMatrix res = new LogicMatrix();
		res.setPtr(biggerrNative(mat1.getPtr(), val));
		return res;
	}
	public static LogicMatrix smaller(Matrix mat1, double val){
		LogicMatrix res = new LogicMatrix();
		res.setPtr(smallerNative(mat1.getPtr(), val));
		return res;
	}
	public static LogicMatrix equal(Matrix mat1, Matrix mat2){
		LogicMatrix res = new LogicMatrix();
		res.setPtr(equalNative(mat1.getPtr(), mat2.getPtr()));
		return res;
	}
	public static LogicMatrix equal(Matrix mat1, double val){
		LogicMatrix res = new LogicMatrix();
		res.setPtr(equalrNative(mat1.getPtr(), val));
		return res;
	}
	public static LogicMatrix or(LogicMatrix mat1, LogicMatrix mat2){
		LogicMatrix res = new LogicMatrix();
		res.setPtr(orNative(mat1.getPtr(), mat2.getPtr()));
		return res;
	}
	public static LogicMatrix and(LogicMatrix mat1, LogicMatrix mat2){
		LogicMatrix res = new LogicMatrix();
		res.setPtr(andNative(mat1.getPtr(), mat2.getPtr()));
		return res;
	}
	public static LogicMatrix not(LogicMatrix mat){
		LogicMatrix res = new LogicMatrix();
		res.setPtr(notNative(mat.getPtr()));
		return res;
	}
	public static Matrix condition(LogicMatrix mat1, Matrix mat2, Matrix mat3){
		Matrix res = new Matrix();
		res.setPtr(conditionNative(mat1.getPtr(), mat2.getPtr(), mat3.getPtr()));
		return res;
	}
	public static LogicMatrix between(Matrix mat, double lowerbound, double upperbound){
		LogicMatrix res = new LogicMatrix();
		res.setPtr(betweenNative(mat.getPtr(), lowerbound, upperbound);
		return res;
	}
	public static Matrix betweenValue(Matrix mat, double lowerbound, double upperbound){
		Matrix res = new Matrix();
		res.setPtr(betweenValueNative(mat.getPtr(), lowerbound, upperbound));
		return res;
	}
	public static Matrix rank(Matrix mat){
		Matrix res = new Matrix();
		res.setPtr(rankNative(mat.getPtr()));
		return res;
	}
	public static Matrix round(Matrix mat){
		Matrix res = new Matrix();
		res.setPtr(roundNative(mat.getPtr()));
		return res;
	}
	public static Matrix floor(Matrix mat){
		Matrix res = new Matrix();
		res.setPtr(floorNative(mat.getPtr()));
		return res;
	}
	public static Matrix abs(Matrix mat){
		Matrix res = new Matrix();
		res.setPtr(absNative(mat.getPtr()));
		return res;
	}
	public static Matrix minus(Matrix mat){
		Matrix res = new Matrix();
		res.setPtr(minusNative(mat.getPtr()));
		return res;
	}
	public static Matrix sqrt(Matrix mat){
		Matrix res = new Matrix();
		res.setPtr(sqrtNative(mat.getPtr()));
		return res;
	}
	public static Matrix log(Matrix mat){
		Matrix res = new Matrix();
		res.setPtr(logNative(mat.getPtr()));
		return res;
	}
	public static Matrix exp(Matrix mat){
		Matrix res = new Matrix();
		res.setPtr(expNative(mat.getPtr()));
		return res;
	}
	public static Matrix sign(Matrix mat){
		Matrix res = new Matrix();
		res.setPtr(signNative(mat.getPtr()));
		return res;
	}
	public static Matrix inverse(Matrix mat){
		Matrix res = new Matrix();
		res.setPtr(inverseNative(mat.getPtr()));
		return res;
	}
	public static Matrix signedpow(Matrix mat, double index){
		Matrix res = new Matrix();
		res.setPtr(signedpowNative(mat.getPtr(), index));
		return res;
	}
	public static Matrix shift(Matrix mat, int n){
		Matrix res = new Matrix();
		res.setPtr(shiftNative(mat.getPtr(), n));
		return res;
	}
	public static Matrix delay(Matrix mat, int n){
		Matrix res = new Matrix();
		res.setPtr(delayNative(mat.getPtr(), n));
		return res;
	}
	public static Matrix delta(Matrix mat, int n){
		Matrix res = new Matrix();
		res.setPtr(deltaNative(mat.getPtr(), n));
		return res;
	}
	public static Matrix ratio(Matrix mat, int n){
		Matrix res = new Matrix();
		res.setPtr(ratioNative(mat.getPtr(), n));
		return res;
	}
	public static Matrix sum(Matrix mat, int n){
		Matrix res = new Matrix();
		res.setPtr(sumNative(mat.getPtr(), n));
		return res;
	}
	public static Matrix product(Matrix mat, int n){
		Matrix res = new Matrix();
		res.setPtr(productNative(mat.getPtr(), n));
		return res;
	}

	public static Matrix tsMax(Matrix mat, int n){
		Matrix res = new Matrix();
		res.setPtr(tsMaxNative(mat.getPtr(), n));
		return res;
	}
	public static Matrix tsMin(Matrix mat, int n){
		Matrix res = new Matrix();
		res.setPtr(tsMinNative(mat.getPtr(), n));
		return res;
	}
	public static Matrix tsArgmax(Matrix mat, int n){
		Matrix res = new Matrix();
		res.setPtr(tsArgmaxNative(mat.getPtr(), n));
		return res;
	}
	public static Matrix tsArgmin(Matrix mat, int n){
		Matrix res = new Matrix();
		res.setPtr(tsArgminNative(mat.getPtr(), n));
		return res;
	}
	public static Matrix tsRank(Matrix mat, int n){
		Matrix res = new Matrix();
		res.setPtr(tsRankNative(mat.getPtr(), n));
		return res;
	}
	public static Matrix tsMean(Matrix mat, int n){
		Matrix res = new Matrix();
		res.setPtr(tsMeanNative(mat.getPtr(), n));
		return res;
	}
	public static Matrix tsStd(Matrix mat, int n){
		Matrix res = new Matrix();
		res.setPtr(tsStdNative(mat.getPtr(), n));
		return res;
	}
	public static Matrix tsSkewness(Matrix mat, int n){
		Matrix res = new Matrix();
		res.setPtr(tsSkewnessNative(mat.getPtr(), n));
		return res;
	}
	public static Matrix tsKurtosis(Matrix mat, int n){
		Matrix res = new Matrix();
		res.setPtr(tsKurtosisNative(mat.getPtr(), n));
		return res;
	}
	public static Matrix tsCov(Matrix mat1, Matrix mat2, int n){
		Matrix res = new Matrix();
		res.setPtr(tsCovNative(mat1.getPtr(), mat2.getPtr(), n));
		return res;
	}
	public static Matrix tsCorr(Matrix mat1, Matrix mat2, int n);{
		Matrix res = new Matrix();
		res.setPtr(tsCorrNative(mat1.getPtr(), mat2.getPtr(), n));
		return res;
	}
	public static Matrix tsCountTrue(LogicMatrix mat, int n){
		Matrix res = new Matrix();
		res.setPtr(tsCountTrueNative(mat.getPtr(), n));
		return res;
	}
	public static Matrix tsCountConsecutiveTrue(LogicMatrix mat, int n){
		Matrix res = new Matrix();
		res.setPtr(tsCountConsecutiveTrueNative(mat.getPtr(), n));
		return res;
	}

	public static Matrix tsCountNaN(Matrix mat, int n){
		Matrix res = new Matrix();
		res.setPtr(tsCountNaNNative(mat.getPtr(), n));
		return res;
	}
	public static Matrix decayLinear(Matrix mat, int n){
		Matrix res = new Matrix();
		res.setPtr(decayLinearNative(mat.getPtr(), n));
		return res;
	}
	public static Matrix decayExponential(Matrix mat, int n){
		Matrix res = new Matrix();
		res.setPtr(decayExponentialNative(mat.getPtr(), n));
		return res;
	}


	public static void smoothByDecayLinear(Matrix mat, int n){
		smoothByDecayLinearNative(mat.getPtr(), n);
	}
	public static void inputNaN(Matrix mat, double val){
		inputNaNNative(mat.getPtr(), val);
	}
	public static void activate(Matrix mat, double threshold){
		activateNative(mat.getPtr(), threshold)
	}

	public static Matrix normalize(Matrix mat, double scale, double mean, double bound){
		Matrix res = new Matrix();
		res.setPtr(normalizeNative(mat.getPtr(), scale, mean, bound));
		return res;
	}

	public static void normalizeBySpec(Matrix mat, double scale, double mean, double bound){
		normalizeBySpecNative(mat.getPtr(), scale, mean, bound);
	}

	public static Matrix neutralize(Matrix mat){
		Matrix res = new Matrix();
		res.setPtr(neutralizeNative(mat.getPtr()));
		return res;
	}
	public static Matrix mean(Matrix mat){
		Matrix res = new Matrix();
		res.setPtr(meanNative(mat.getPtr()));
		return res;
	}
	public static Matrix unify(Matrix mat){
		Matrix res = new Matrix();
		res.setPtr(unifyNative(mat.getPtr()));
		return res;
	}
	public static Matrix unifyByL2(Matrix mat){
		Matrix res = new Matrix();
		res.setPtr(unifyByL2Native(mat.getPtr()));
		return res;
	}
	public static Matrix evalValidPct(Matrix mat){
		Matrix res = new Matrix();
		res.setPtr(evalValidPctNative(mat.getPtr()));
		return res;
	}
	public static Matrix evalAbsSum(Matrix mat){
		Matrix res = new Matrix();
		res.setPtr(evalAbsSumNative(mat.getPtr()));
		return res;
	}
	public static Matrix evalMean(Matrix mat){
		Matrix res = new Matrix();
		res.setPtr(evalMeanNative(mat.getPtr()));
		return res;
	}
	public static Matrix evalVariance(Matrix mat){
		Matrix res = new Matrix();
		res.setPtr(evalVarianceNative(mat.getPtr()));
		return res;
	}
	public static Matrix evalInnerProduction(Matrix alpha, Matrix target){
		Matrix res = new Matrix();
		res.setPtr(evalInnerProductionNative(alpha.getPtr()), target.getPtr());
		return res;
	}
	public static Matrix evalCovariance(Matrix alpha, Matrix target){
		Matrix res = new Matrix();
		res.setPtr(evalCovarianceNative(alpha.getPtr()), target.getPtr());
		return res;
	}
	public static Matrix evalCorrelation(Matrix alpha, Matrix target){
		Matrix res = new Matrix();
		res.setPtr(evalCorrelationNative(alpha.getPtr()), target.getPtr());
		return res;
	}

	public static double Det(Matrix mat, int N){////REMINDER : N = mat.getNRow() - 1
		double res = Det(mat.getPtr(), N));
		return res;
	}
	public static double Inverse(Matrix mat1, int N, Matrix mat3){
		Matrix res = new Matrix();
		res.setPtr(InverseNative(mat1.getPtr()), N, mat3.getPtr());
		return res;
	}
	public static Matrix inv(Matrix mat){
		Matrix res = new Matrix();
		res.setPtr(invNative(mat.getPtr()));
		return res;
	}
	public static double treat(Matrix mat){
		double res = treatNative(mat.getPtr()));
		return res;
	}
	public static Matrix diag(Matrix mat){
		Matrix res = new Matrix();
		res = diagNative(mat.getPtr());
		return res;
	}
	public static Matrix inverseDiag(Matrix mat){
		Matrix res = new Matrix();
		res = inverseDiagNative(mat.getPtr());
		return res;
	}
	public static double evalInnerProductionByLongVector(Matrix alpha, Matrix target){
		double res = evalInnerProductionByLongVectorNative(alpha.getPtr()), target.getPtr());
		return res;
	}
	public static double evalCorrelationByLongVector(Matrix alpha, Matrix target){
		double res = evalCorrelationByLongVectorNative(alpha.getPtr()), target.getPtr());
		return res;
	}
	public static Matrix evalBeta(Matrix alpha, Matrix target){
		Matrix res = new Matrix();
		res = evalBetaNative(alpha.getPtr(), target.getPtr());
		return res;
	}
	public static Matrix cumSum(Matrix ts){
		Matrix res = new Matrix();
		res = cumSumNative(ts.getPt());
		return res;
	}
	public static double evalBetaByLongVector(Matrix alpha, Matrix target){
		double res = evalBetaByLongVectorNative(alpha.getPtr()), target.getPtr());
		return res;
	}

	public static Matrix add(Matrix mat1, Matrix mat2, int num){
		Matrix res = new Matrix();
		res.setPtr(addNative(mat1.getPtr(), mat2.getPtr(), num));
		return res;
	}
	public static Matrix sub(Matrix mat1, Matrix mat2, int num){
		Matrix res = new Matrix();
		res.setPtr(subNative(mat1.getPtr(), mat2.getPtr(), num));
		return res;
	}
	public static Matrix div(Matrix mat1, Matrix mat2, int num){
		Matrix res = new Matrix();
		res.setPtr(divNative(mat1.getPtr(), mat2.getPtr(), num));
		return res;
	}
	public static Matrix mul(Matrix mat1, Matrix mat2, int num){
		Matrix res = new Matrix();
		res.setPtr(mulNative(mat1.getPtr(), mat2.getPtr(), num));
		return res;
	}
	public static Matrix mul(double val, Matrix mat, int num){
		Matrix res = new Matrix();
		res.setPtr(addNative(val, mat.getPtr(), num));
		return res;
	}
	public static Matrix matrixMul(Matrix mat1, Matrix mat2, int num){
		Matrix res = new Matrix();
		res.setPtr(matrixMulNative(mat1.getPtr(), mat2.getPtr(), num));
		return res;
	}

	public static Matrix max(Matrix mat1, Matrix mat2, int num){
		Matrix res = new Matrix();
		res.setPtr(maxNative(mat1.getPtr(), mat2.getPtr(), num));
		return res;
	}
	public static Matrix min(Matrix mat1, Matrix mat2, int num){
		Matrix res = new Matrix();
		res.setPtr(minNative(mat1.getPtr(), mat2.getPtr(), num));
		return res;
	}
	public static LogicMatrix bigger(Matrix mat1, Matrix mat2, int num){
		LogicMatrix res = new LogicMatrix();
		res.setPtr(biggerNative(mat1.getPtr(), mat2.getPtr(), num));
		return res;
	}
	public static LogicMatrix smaller(Matrix mat1, Matrix mat2, int num){
		LogicMatrix res = new LogicMatrix();
		res.setPtr(smallerNative(mat1.getPtr(), mat2.getPtr(), num));
		return res;
	}
	public static LogicMatrix bigger(Matrix mat1, double val, int num){
		LogicMatrix res = new LogicMatrix();
		res.setPtr(biggerrNative(mat1.getPtr(), val, num));
		return res;
	}
	public static LogicMatrix smaller(Matrix mat1, double val, int num){
		LogicMatrix res = new LogicMatrix();
		res.setPtr(smallerNative(mat1.getPtr(), val, num));
		return res;
	}
	public static LogicMatrix equal(Matrix mat1, Matrix mat2, int num){
		LogicMatrix res = new LogicMatrix();
		res.setPtr(equalNative(mat1.getPtr(), mat2.getPtr(), num));
		return res;
	}
	public static LogicMatrix equal(Matrix mat1, double val, int num){
		LogicMatrix res = new LogicMatrix();
		res.setPtr(equalrNative(mat1.getPtr(), val, num));
		return res;
	}
	public static LogicMatrix or(LogicMatrix mat1, LogicMatrix mat2, int num){
		LogicMatrix res = new LogicMatrix();
		res.setPtr(orNative(mat1.getPtr(), mat2.getPtr(), num));
		return res;
	}
	public static LogicMatrix and(LogicMatrix mat1, LogicMatrix mat2, int num){
		LogicMatrix res = new LogicMatrix();
		res.setPtr(andNative(mat1.getPtr(), mat2.getPtr(), num));
		return res;
	}
	public static LogicMatrix not(LogicMatrix mat, int num){
		LogicMatrix res = new LogicMatrix();
		res.setPtr(notNative(mat.getPtr(), num));
		return res;
	}
	public static Matrix condition(LogicMatrix mat1, Matrix mat2, Matrix mat3, int num){
		Matrix res = new Matrix();
		res.setPtr(conditionNative(mat1.getPtr(), mat2.getPtr(), mat3.getPtr(), num));
		return res;
	}
	public static Matrix rank(Matrix mat, int num){
		Matrix res = new Matrix();
		res.setPtr(rankNative(mat.getPtr(), num));
		return res;
	}
	public static Matrix round(Matrix mat, int num){
		Matrix res = new Matrix();
		res.setPtr(roundNative(mat.getPtr(), num));
		return res;
	}
	public static Matrix floor(Matrix mat, int num){
		Matrix res = new Matrix();
		res.setPtr(floorNative(mat.getPtr(), num));
		return res;
	}
	public static Matrix abs(Matrix mat, int num){
		Matrix res = new Matrix();
		res.setPtr(absNative(mat.getPtr(), num));
		return res;
	}
	public static Matrix minus(Matrix mat, int num){
		Matrix res = new Matrix();
		res.setPtr(minusNative(mat.getPtr(), num));
		return res;
	}
	public static Matrix sqrt(Matrix mat, int num){
		Matrix res = new Matrix();
		res.setPtr(sqrtNative(mat.getPtr(), num));
		return res;
	}
	public static Matrix log(Matrix mat, int num){
		Matrix res = new Matrix();
		res.setPtr(logNative(mat.getPtr(), num));
		return res;
	}
	public static Matrix exp(Matrix mat, int num){
		Matrix res = new Matrix();
		res.setPtr(expNative(mat.getPtr(), num));
		return res;
	}
	public static Matrix sign(Matrix mat, int num){
		Matrix res = new Matrix();
		res.setPtr(signNative(mat.getPtr(), num));
		return res;
	}
	public static Matrix inverse(Matrix mat, int num){
		Matrix res = new Matrix();
		res.setPtr(inverseNative(mat.getPtr(), num));
		return res;
	}
	public static Matrix signedpow(Matrix mat, double index, int num){
		Matrix res = new Matrix();
		res.setPtr(signedpowNative(mat.getPtr(), index, num));
		return res;
	}
	public static Matrix shift(Matrix mat, int n, int num){
		Matrix res = new Matrix();
		res.setPtr(shiftNative(mat.getPtr(), n, num));
		return res;
	}
	public static Matrix delay(Matrix mat, int n, int num){
		Matrix res = new Matrix();
		res.setPtr(delayNative(mat.getPtr(), n, num));
		return res;
	}
	public static Matrix delta(Matrix mat, int n, int num){
		Matrix res = new Matrix();
		res.setPtr(deltaNative(mat.getPtr(), n, num));
		return res;
	}
	public static Matrix ratio(Matrix mat, int n, int num){
		Matrix res = new Matrix();
		res.setPtr(ratioNative(mat.getPtr(), n, num));
		return res;
	}
	public static Matrix sum(Matrix mat, int n, int num){
		Matrix res = new Matrix();
		res.setPtr(sumNative(mat.getPtr(), n, num));
		return res;
	}
	public static Matrix product(Matrix mat, int n, int num){
		Matrix res = new Matrix();
		res.setPtr(productNative(mat.getPtr(), n, num));
		return res;
	}

	public static Matrix tsMax(Matrix mat, int n, int num){
		Matrix res = new Matrix();
		res.setPtr(tsMaxNative(mat.getPtr(), n, num));
		return res;
	}
	public static Matrix tsMin(Matrix mat, int n, int num){
		Matrix res = new Matrix();
		res.setPtr(tsMinNative(mat.getPtr(), n, num));
		return res;
	}
	public static Matrix tsArgmax(Matrix mat, int n, int num){
		Matrix res = new Matrix();
		res.setPtr(tsArgmaxNative(mat.getPtr(), n, num));
		return res;
	}
	public static Matrix tsArgmin(Matrix mat, int n, int num){
		Matrix res = new Matrix();
		res.setPtr(tsArgminNative(mat.getPtr(), n, num));
		return res;
	}
	public static Matrix tsRank(Matrix mat, int n, int num){
		Matrix res = new Matrix();
		res.setPtr(tsRankNative(mat.getPtr(), n, num));
		return res;
	}
	public static Matrix tsMean(Matrix mat, int n, int num){
		Matrix res = new Matrix();
		res.setPtr(tsMeanNative(mat.getPtr(), n, num));
		return res;
	}
	public static Matrix tsStd(Matrix mat, int n, int num){
		Matrix res = new Matrix();
		res.setPtr(tsStdNative(mat.getPtr(), n, num));
		return res;
	}
	public static Matrix tsSkewness(Matrix mat, int n, int num){
		Matrix res = new Matrix();
		res.setPtr(tsSkewnessNative(mat.getPtr(), n, num));
		return res;
	}
	public static Matrix tsKurtosis(Matrix mat, int n, int num){
		Matrix res = new Matrix();
		res.setPtr(tsKurtosisNative(mat.getPtr(), n, num));
		return res;
	}
	public static Matrix tsCov(Matrix mat1, Matrix mat2, int n, int num){
		Matrix res = new Matrix();
		res.setPtr(tsCovNative(mat1.getPtr(), mat2.getPtr(), n, num));
		return res;
	}
	public static Matrix tsCorr(Matrix mat1, Matrix mat2, int n){
		Matrix res = new Matrix();
		res.setPtr(tsCorrNative(mat1.getPtr(), mat2.getPtr(), n, num));
		return res;
	}
	public static Matrix tsCountTrue(LogicMatrix mat, int n, int num){
		Matrix res = new Matrix();
		res.setPtr(tsCountTrueNative(mat.getPtr(), n, num));
		return res;
	}
	public static Matrix tsCountConsecutiveTrue(LogicMatrix mat, int n, int num){
		Matrix res = new Matrix();
		res.setPtr(tsCountConsecutiveTrueNative(mat.getPtr(), n, num));
		return res;
	}

	public static Matrix tsCountNaN(Matrix mat, int n, int num){
		Matrix res = new Matrix();
		res.setPtr(tsCountNaNNative(mat.getPtr(), n, num));
		return res;
	}
	public static Matrix decayLinear(Matrix mat, int n, int num){
		Matrix res = new Matrix();
		res.setPtr(decayLinearNative(mat.getPtr(), n, num));
		return res;
	}
	public static Matrix decayExponential(Matrix mat, int n, int num){
		Matrix res = new Matrix();
		res.setPtr(decayExponentialNative(mat.getPtr(), n, num));
		return res;
	}


	public static Matrix normalize(Matrix mat, double scale, double mean, double bound, int num){
		Matrix res = new Matrix();
		res.setPtr(normalizeNative(mat.getPtr(), scale, mean, bound, num));
		return res;
	}


	public static Matrix neutralize(Matrix mat, int num){
		Matrix res = new Matrix();
		res.setPtr(neutralizeNative(mat.getPtr(), num));
		return res;
	}
	public static Matrix mean(Matrix mat, int num){
		Matrix res = new Matrix();
		res.setPtr(meanNative(mat.getPtr(), num));
		return res;
	}
	public static Matrix unify(Matrix mat, int num){
		Matrix res = new Matrix();
		res.setPtr(unifyNative(mat.getPtr(), num));
		return res;
	}
	public static Matrix unifyByL2(Matrix mat, int num){
		Matrix res = new Matrix();
		res.setPtr(unifyByL2Native(mat.getPtr(), num));
		return res;
	}


	private native long createMatrixCalculator();

	private static native long addNative(long mat1, long mat2);
	private static native long subNative(long mat1, long mat2);
	private static native long divNative(long mat1, long mat2);
	private static native long mulNative(long mat1, long mat2);	
	private static native long mulNative(double val, long mat);
	private static native long matrixMulNative(long mat1, long mat2);

	private static native long maxNative(long mat1, long mat2);
	private static native long minNative(long mat1, long mat2);
	private static native long biggerNative(long mat1, long mat2);
	private static native long smallerNative(long mat1, long mat2);
	private static native long biggerNative(long mat1, double val);
	private static native long smallerNative(long mat1, double val);
	private static native long equalNative(long mat1, long mat2);
	private static native long equalNative(long mat1, double val);
	private static native long orNative(long mat1, long mat2);
	private static native long andNative(long mat1, long mat2);
	private static native long notNative(long mat);
	private static native long conditionNative(long mat1, long mat2, long mat3);
	private static native long betweenNative(long mat, double lowerbound, double upperbound);
	private static native long betweenValueNative(long mat, double lowerbound, double upperbound);

	private static native long rankNative(long mat);
	private static native long roundNative(long mat);
	private static native long floorNative(long mat);
	private static native long absNative(long mat);
	private static native long minusNative(long mat);
	private static native long sqrtNative(long mat);
	private static native long logNative(long mat);
	private static native long expNative(long mat);
	private static native long signNative(long mat);
	private static native long inverseNative(long mat);
	private static native long signedpowNative(long mat, double index);
	private static native long shiftNative(long mat, int n);
	private static native long delayNative(long mat, int n);
	private static native long deltaNative(long mat, int n);
	private static native long ratioNative(long mat, int n);
	private static native long sumNative(long mat, int n);
	private static native long productNative(long mat, int n);

	private static native long tsMaxNative(long mat, int n);
	private static native long tsMinNative(long mat, int n);
	private static native long tsArgmaxNative(long mat, int n);
	private static native long tsArgminNative(long mat, int n);
	private static native long tsRankNative(long mat, int n);
	private static native long tsMeanNative(long mat, int n);
	private static native long tsStdNative(long mat, int n);
	private static native long tsSkewnessNative(long mat, int n);
	private static native long tsKurtosisNative(long mat, int n);
	private static native long tsCovNative(long mat, int n);
	private static native long tsCorrNative(long mat, int n);
	private static native long tsCountNaNNative(long mat, int n);
	private static native long tsCountTrueNative(long mat, int n);
	private static native long tsCountConsecutiveTrueNative(long mat, int n);
	private static native long decayLinearNative(long mat, int n);
	private static native long decayExponentialNative(long mat, int n);


	private static native void smoothByDecayLinearNative(long mat, int n);
	private static native void inputNaNNative(long mat, double val);
	private static native void activateNative(long mat, double threashold);

	private static native long normalizeNative(long mat, double scale, double mean, double bound);
	private static native long neutralizeNative(long mat);
	private static native long meanNative(long mat);
	private static native long unifyNative(long mat);
	private static native long unifyByL2Native(long mat);
	private static native long evalValidPctNative(long mat);
	private static native long evalAbsSumNative(long mat);
	private static native long evalMeanNative(long mat);
	private static native long evalVarianceNative(long mat);
	private static native long evalInnerProductionNative(long alpha, long target);
	private static native long evalCovarianceNative(long alpha, long target);
	private static native long evalCorrelationNatve(long alpha, long target);

	private static native double DetNative(long mat, int N);////REMINDER : N = mat.getNRow() - 1
	private static native double InverseNative(long mat1, int N, long mat3);
	private static native long invNative(long mat);
	private static native double treatNative(long mat);
	private static native long diagNative(long mat);
	private static native long inverseDiagNative(long mat);
	private static native double evalInnerProductionByLongVectorNative(long alpha, long target);
	private static native double evalCorrelationByLongVectorNative(long alpha, long target);
	private static native long evalBetaNative(long alpha, long target);
	private static native long cumSumNative(long ts);
	private static native double evalBetaByLongVectorNative(long alpha, long target);






	private static native long addNative(long mat1, long mat2, int num);
	private static native long subNative(long mat1, long mat2, int num);
	private static native long divNative(long mat1, long mat2, int num);
	private static native long mulNative(long mat1, long mat2, int num);	
	private static native long mulNative(double val, long mat, int num);
	private static native long matrixMulNative(long mat1, long mat2, int num);

	private static native long maxNative(long mat1, long mat2, int num);
	private static native long minNative(long mat1, long mat2, int num);
	private static native long biggerNative(long mat1, long mat2, int num);
	private static native long smallerNative(long mat1, long mat2, int num);
	private static native long biggerNative(long mat1, double val, int num);
	private static native long smallerNative(long mat1, double val, int num);
	private static native long equalNative(long mat1, long mat2, int num);
	private static native long equalNative(long mat1, double val, int num);
	private static native long orNative(long mat1, long mat2, int num);
	private static native long andNative(long mat1, long mat2, int num);
	private static native long notNative(long mat, int num);
	private static native long conditionNative(long mat1, long mat2, long mat3, int num);
	private static native long rankNative(long mat, int num);
	private static native long roundNative(long mat, int num);
	private static native long floorNative(long mat, int num);
	private static native long absNative(long mat, int num);
	private static native long minusNative(long mat, int num);
	private static native long sqrtNative(long mat, int num);
	private static native long logNative(long mat, int num);
	private static native long expNative(long mat, int num);
	private static native long signNative(long mat, int num);
	private static native long inverseNative(long mat, int num);
	private static native long signedpowNative(long mat, double index, int num);
	private static native long shiftNative(long mat, int n, int num);
	private static native long delayNative(long mat, int n, int num);
	private static native long deltaNative(long mat, int n, int num);
	private static native long ratioNative(long mat, int n, int num);
	private static native long sumNative(long mat, int n, int num);
	private static native long productNative(long mat, int n, int num);

	private static native long tsMaxNative(long mat, int n, int num);
	private static native long tsMinNative(long mat, int n, int num);
	private static native long tsArgmaxNative(long mat, int n, int num);
	private static native long tsArgminNative(long mat, int n, int num);
	private static native long tsRankNative(long mat, int n, int num);
	private static native long tsMeanNative(long mat, int n, int num);
	private static native long tsStdNative(long mat, int n, int num);
	private static native long tsSkewnessNative(long mat, int n, int num);
	private static native long tsKurtosisNative(long mat, int n, int num);
	private static native long tsCovNative(long mat, int n, int num);
	private static native long tsCorrNative(long mat, int n, int num);
	private static native long tsCountNaNNative(long mat, int n, int num);
	private static native long tsCountTrueNative(long mat, int n, int num);
	private static native long tsCountConsecutiveTrueNative(long mat, int n, int num);
	private static native long decayLinearNative(long mat, int n, int num);
	private static native long decayExponentialNative(long mat, int n, int num);
	private static native long normalizeNative(long mat, double scale, double mean, double bound, int num);
	private static native long neutralizeNative(long mat, int num);
	


}
