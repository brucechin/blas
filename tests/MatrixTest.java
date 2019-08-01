/**
 * @CLassName MatrixTest
 * @Description This class is used for testing if Java/C++ MatrixCalculator will output the same
 * @Author BruceChin
 * @Date 2019/7/30 5:55 PM
 * @Version 1.0
 **/
public class MatrixTest {
    public static void main(String[] args){

        int rangeMin = -10000;
        int rangeMax = 10000;
        int matSize = 1000;
        int lowerBound = -100;
        int upperBound = 100;
        int stepSize = 50; // used for shift/delta/delay/ratio/sum/product
        MatrixFactory factory = new MatrixFactory();
        MatrixCalculator calculator = new MatrixCalculator();
        String matrix_a = "a.mat";
        String matrix_b = "b.mat";
        String logicMatrix_a = "la.mat";
        String logicMatrix_b = "lb.mat";
        int counter = 0;

        Matrix a = new Matrix(matrix_a);
        Matrix b = new Matrix(matrix_b);
        LogicMatrix la = new LogicMatrix(logicMatrix_a);
        LogicMatrix lb = new LogicMatrix(logicMatrix_b);

        Matrix res;
        LogicMatrix lres;
        Matrix compareMatrix;
        LogicMatrix compareLogicMatrix;

        System.out.println("Det(a) : " + calculator.Det(a.value, a.getNRow()));

        System.out.println("treat(a) : " + calculator.treat(a));

        res = calculator.Normalize(a, 2.0, 10.0, 1000.0);
        compareMatrix.readCppMatrix("normalize.mat");
        res.compareMat(compareMatrix);
        res.clear();
        
        calculator->smoothByDecayLinear(a, 5);
        compareMatrix.readCppMatrix("smooth.mat");
        a.compareMat(compareMatrix);

        res = calculator.add(a, b);
        compareMatrix.readCppMatrix("add.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.sub(a, b);
        compareMatrix.readCppMatrix("sub.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.div(a, b);
        compareMatrix.readCppMatrix("div.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.mul(a, b);
        compareMatrix.readCppMatrix("mul.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.matrixMul(a, b);
        compareMatrix.readCppMatrix("matrixMul.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.max(a, b);
        compareMatrix.readCppMatrix("max.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.min(a, b);
        compareMatrix.readCppMatrix("min.mat");
        res.compareMat(compareMatrix);
        res.clear();

        lres = calculator.bigger(a, b);
        compareLogicMatrix.readCppMatrix("bigger.mat");
        lres.compareLogicMat(compareLogicMatrix);
        lres.clear();

        lres = calculator.smaller(a, b);
        compareLogicMatrix.readCppMatrix("smaller.mat");
        lres.compareLogicMat(compareLogicMatrix);
        lres.clear();


        lres = calculator.equal(a, b);
        compareLogicMatrix.readCppMatrix("equal.mat");
        lres.compareLogicMat(compareLogicMatrix);
        lres.clear();

        lres = calculator.between(a, lowerBound, upperBound);
        compareLogicMatrix.readCppMatrix("between.mat");
        lres.compareLogicMat(compareLogicMatrix);
        lres.clear();

        lres = calculator.matAnd(la, lb);
        compareLogicMatrix.readCppMatrix("matAnd.mat");
        lres.compareLogicMat(compareLogicMatrix);
        lres.clear();

        lres = calculator.matOr(la, lb);
        compareLogicMatrix.readCppMatrix("matOr.mat");
        lres.compareLogicMat(compareLogicMatrix);
        lres.clear();

        lres = calculator.matNot(la);
        compareLogicMatrix.readCppMatrix("matNot.mat");
        lres.compareLogicMat(compareLogicMatrix);
        lres.clear();

        res = calculator.condition(la, a, b);
        compareMatrix.readCppMatrix("condition.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.rank(a);
        compareMatrix.readCppMatrix("rank.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.round(a);
        compareMatrix.readCppMatrix("round.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.floor(a);
        compareMatrix.readCppMatrix("floor.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.abs(a);
        compareMatrix.readCppMatrix("abs.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.minus(a);
        compareMatrix.readCppMatrix("minus.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.sqrt(a);
        compareMatrix.readCppMatrix("sqrt.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.log(a);
        compareMatrix.readCppMatrix("log.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.exp(a);
        compareMatrix.readCppMatrix("exp.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.sign(a);
        compareMatrix.readCppMatrix("sign.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.inverse(a);
        compareMatrix.readCppMatrix("inverse.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.signedpow(a, 2);
        compareMatrix.readCppMatrix("signedpow.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.shift(a, stepSize);
        compareMatrix.readCppMatrix("shift.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.delay(a, stepSize);
        compareMatrix.readCppMatrix("delay.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.delta(a, stepSize);
        compareMatrix.readCppMatrix("delta.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.ratio(a, stepSize);
        compareMatrix.readCppMatrix("ratio.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.sum(a, stepSize);
        compareMatrix.readCppMatrix("sum.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.product(a, stepSize);
        compareMatrix.readCppMatrix("product.mat");
        res.compareMat(compareMatrix);
        res.clear();


        res = calculator.tsMax(a, stepSize);
        compareMatrix.readCppMatrix("tsmax.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.tsMin(a, stepSize);
        compareMatrix.readCppMatrix("tsmin.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.tsArgmax(a, stepSize);
        compareMatrix.readCppMatrix("tsargmax.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.tsArgmin(a, stepSize);
        compareMatrix.readCppMatrix("tsargmin.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.tsRank(a, stepSize);
        compareMatrix.readCppMatrix("tsrank.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.tsMean(a, stepSize);
        compareMatrix.readCppMatrix("tsmean.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.tsStd(a, stepSize);
        compareMatrix.readCppMatrix("tsstd.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.tsSkewness(a, stepSize);
        compareMatrix.readCppMatrix("tsskewness.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.tsKurtosis(a, stepSize);
        compareMatrix.readCppMatrix("tskurtosis.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.tsCov(a, b, stepSize);
        compareMatrix.readCppMatrix("tscov.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.tsCorr(a, b, stepSize);
        compareMatrix.readCppMatrix("tscorr.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.tsCountTrue(la, stepSize);
        compareMatrix.readCppMatrix("tscounttrue.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.tsCountConsecutiveTrue(la, stepSize);
        compareMatrix.readCppMatrix("tscountconsecutivetrue.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.decayLinear(a, stepSize);
        compareMatrix.readCppMatrix("decaylinear.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.decayExponential(a, stepSize);
        compareMatrix.readCppMatrix("decayexponential.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.neutralize(a);
        compareMatrix.readCppMatrix("neutralize.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.unify(a);
        compareMatrix.readCppMatrix("unify.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.evalValidPct(a);
        compareMatrix.readCppMatrix("evalvalidpct.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.evalAbsSum(a);
        compareMatrix.readCppMatrix("evalabssum.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.evalMean(a);
        compareMatrix.readCppMatrix("evalmean.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.evalVariance(a);
        compareMatrix.readCppMatrix("evalvariance.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.evalInnerProduction(a, b);
        compareMatrix.readCppMatrix("evalinnerproduction.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.evalCovariance(a, b);
        compareMatrix.readCppMatrix("evalcovariance.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.evalCorrelation(a, b);
        compareMatrix.readCppMatrix("evalcorrelation.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.evalBeta(a, b);
        compareMatrix.readCppMatrix("evalbeta.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.diag(a);
        compareMatrix.readCppMatrix("diag.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.cumSum(a);
        compareMatrix.readCppMatrix("cumsum.mat");
        res.compareMat(compareMatrix);
        res.clear();

    }
}
