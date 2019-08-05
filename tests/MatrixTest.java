/**
 * @CLassName MatrixTest
 * @Description This class is used for testing if Java/C++ MatrixCalculator will output the same
 * @Author BruceChin
 * @Date 2019/7/30 5:55 PM
 * @Version 1.0
 **/

package matrix;

import matrix.*;

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
        Matrix compareMatrix = new Matrix();
        LogicMatrix compareLogicMatrix = new LogicMatrix();

        //System.out.println("Det(a) : " + calculator.Det(a.value, a.getNRow()));

        System.out.println("treat(a) : " + calculator.treat(a));

        res = calculator.Normalize(a, 2.0, 10.0, 1000.0);
        compareMatrix.readMatrix("normalize.mat");
        res.compareMat(compareMatrix);
        res.clear();
        
        calculator->smoothByDecayLinear(a, 5);
        compareMatrix.readMatrix("smooth.mat");
        a.compareMat(compareMatrix);

        res = calculator.add(a, b);
        compareMatrix.readMatrix("add.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.sub(a, b);
        compareMatrix.readMatrix("sub.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.div(a, b);
        compareMatrix.readMatrix("div.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.mul(a, b);
        compareMatrix.readMatrix("mul.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.matrixMul(a, b);
        compareMatrix.readMatrix("matrixMul.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.max(a, b);
        compareMatrix.readMatrix("max.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.min(a, b);
        compareMatrix.readMatrix("min.mat");
        res.compareMat(compareMatrix);
        res.clear();

        lres = calculator.bigger(a, b);
        compareLogicMatrix.readMatrix("bigger.mat");
        lres.compareLogicMat(compareLogicMatrix);
        lres.clear();

        lres = calculator.smaller(a, b);
        compareLogicMatrix.readMatrix("smaller.mat");
        lres.compareLogicMat(compareLogicMatrix);
        lres.clear();


        lres = calculator.equal(a, b);
        compareLogicMatrix.readMatrix("equal.mat");
        lres.compareLogicMat(compareLogicMatrix);
        lres.clear();

        lres = calculator.between(a, lowerBound, upperBound);
        compareLogicMatrix.readMatrix("between.mat");
        lres.compareLogicMat(compareLogicMatrix);
        lres.clear();

        lres = calculator.and(la, lb);
        compareLogicMatrix.readMatrix("matAnd.mat");
        lres.compareLogicMat(compareLogicMatrix);
        lres.clear();

        lres = calculator.or(la, lb);
        compareLogicMatrix.readMatrix("matOr.mat");
        lres.compareLogicMat(compareLogicMatrix);
        lres.clear();

        lres = calculator.not(la);
        compareLogicMatrix.readMatrix("matNot.mat");
        lres.compareLogicMat(compareLogicMatrix);
        lres.clear();

        res = calculator.condition(la, a, b);
        compareMatrix.readMatrix("condition.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.rank(a);
        compareMatrix.readMatrix("rank.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.round(a);
        compareMatrix.readMatrix("round.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.floor(a);
        compareMatrix.readMatrix("floor.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.abs(a);
        compareMatrix.readMatrix("abs.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.minus(a);
        compareMatrix.readMatrix("minus.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.sqrt(a);
        compareMatrix.readMatrix("sqrt.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.log(a);
        compareMatrix.readMatrix("log.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.exp(a);
        compareMatrix.readMatrix("exp.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.sign(a);
        compareMatrix.readMatrix("sign.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.inverse(a);
        compareMatrix.readMatrix("inverse.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.signedpow(a, 2);
        compareMatrix.readMatrix("signedpow.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.shift(a, stepSize);
        compareMatrix.readMatrix("shift.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.delay(a, stepSize);
        compareMatrix.readMatrix("delay.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.delta(a, stepSize);
        compareMatrix.readMatrix("delta.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.ratio(a, stepSize);
        compareMatrix.readMatrix("ratio.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.sum(a, stepSize);
        compareMatrix.readMatrix("sum.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.product(a, stepSize);
        compareMatrix.readMatrix("product.mat");
        res.compareMat(compareMatrix);
        res.clear();


        res = calculator.tsmax(a, stepSize);
        compareMatrix.readMatrix("tsmax.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.tsmin(a, stepSize);
        compareMatrix.readMatrix("tsmin.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.tsargmax(a, stepSize);
        compareMatrix.readMatrix("tsargmax.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.tsargmin(a, stepSize);
        compareMatrix.readMatrix("tsargmin.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.tsrank(a, stepSize);
        compareMatrix.readMatrix("tsrank.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.tsmean(a, stepSize);
        compareMatrix.readMatrix("tsmean.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.tsstd(a, stepSize);
        compareMatrix.readMatrix("tsstd.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.tsskewness(a, stepSize);
        compareMatrix.readMatrix("tsskewness.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.tskurtosis(a, stepSize);
        compareMatrix.readMatrix("tskurtosis.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.tscov(a, b, stepSize);
        compareMatrix.readMatrix("tscov.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.tscorr(a, b, stepSize);
        compareMatrix.readMatrix("tscorr.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.tscounttrue(la, stepSize);
        compareMatrix.readMatrix("tscounttrue.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.tscountconsecutivetrue(la, stepSize);
        compareMatrix.readMatrix("tscountconsecutivetrue.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.decaylinear(a, stepSize);
        compareMatrix.readMatrix("decaylinear.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.decayexponential(a, stepSize);
        compareMatrix.readMatrix("decayexponential.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.neutralize(a);
        compareMatrix.readMatrix("neutralize.mat");
        res.compareMat(compareMatrix);
        res.clear();

        res = calculator.unify(a);
        compareMatrix.readMatrix("unify.mat");
        res.compareMat(compareMatrix);
        res.clear();

        // res = calculator.evalValidPct(a);
        // compareMatrix.readMatrix("evalvalidpct.mat");
        // res.compareMat(compareMatrix);
        // res.clear();

        // res = calculator.evalAbsSum(a);
        // compareMatrix.readMatrix("evalabssum.mat");
        // res.compareMat(compareMatrix);
        // res.clear();

        // res = calculator.evalMean(a);
        // compareMatrix.readMatrix("evalmean.mat");
        // res.compareMat(compareMatrix);
        // res.clear();

        // res = calculator.evalVariance(a);
        // compareMatrix.readMatrix("evalvariance.mat");
        // res.compareMat(compareMatrix);
        // res.clear();

        // res = calculator.evalInnerProduction(a, b);
        // compareMatrix.readMatrix("evalinnerproduction.mat");
        // res.compareMat(compareMatrix);
        // res.clear();

        // res = calculator.evalCovariance(a, b);
        // compareMatrix.readMatrix("evalcovariance.mat");
        // res.compareMat(compareMatrix);
        // res.clear();

        // res = calculator.evalCorrelation(a, b);
        // compareMatrix.readMatrix("evalcorrelation.mat");
        // res.compareMat(compareMatrix);
        // res.clear();

        // res = calculator.evalBeta(a, b);
        // compareMatrix.readMatrix("evalbeta.mat");
        // res.compareMat(compareMatrix);
        // res.clear();

        // res = calculator.diag(a);
        // compareMatrix.readMatrix("diag.mat");
        // res.compareMat(compareMatrix);
        // res.clear();

        // res = calculator.cumSum(a);
        // compareMatrix.readMatrix("cumsum.mat");
        // res.compareMat(compareMatrix);
        // res.clear();

    }
}
