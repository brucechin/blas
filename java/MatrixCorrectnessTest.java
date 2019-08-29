/**
 * @CLassName MatrixTest
 * @Description This class is used for testing if Java/C++ MatrixMatrixCalculator will output the same
 * @Author BruceChin
 * @Date 2019/7/30 5:55 PM
 * @Version 1.0
 **/

package blas.java;

import blas.java.*;

public class MatrixCorrectnessTest {
        public static void main(String[] args) {

                int rangeMin = -10000;
                int rangeMax = 10000;
                int matSize = 1000;
                int lowerBound = -100;
                int upperBound = 100;
                int stepSize = 50; // used for shift/delta/delay/ratio/sum/product
                // MatrixFactory factory = new MatrixFactory();
                // MatrixCalculator calculator = new MatrixCalculator();
                String matrix_a = "a.mat";
                String matrix_b = "b.mat";
                String logicMatrix_a = "la.mat";
                String logicMatrix_b = "lb.mat";

                Matrix a = new Matrix(matrix_a);
                Matrix b = new Matrix(matrix_b);
                LogicMatrix la = new LogicMatrix(logicMatrix_a);
                LogicMatrix lb = new LogicMatrix(logicMatrix_b);

                Matrix res;
                LogicMatrix lres;
                Matrix compareMatrix = new Matrix();
                LogicMatrix compareLogicMatrix = new LogicMatrix();
                /*
                 * System.out.println("Det(a) : " + MatrixCalculator.Det(a.value, a.getNRow() -
                 * 1));
                 * 
                 * System.out.println("treat(a) : " + MatrixCalculator.treat(a));
                 * 
                 */
                int i = 0;

                res = MatrixCalculator.add(a, b);
                compareMatrix.readMatrix("add.mat");
                res.compareMat(compareMatrix);
                res.clear();
                System.out.println(i++);

                res = MatrixCalculator.sub(a, b);
                compareMatrix.readMatrix("sub.mat");
                res.compareMat(compareMatrix);
                res.clear();
                System.out.println(i++);

                res = MatrixCalculator.div(a, b);
                compareMatrix.readMatrix("div.mat");
                res.compareMat(compareMatrix);
                res.clear();
                System.out.println(i++);

                res = MatrixCalculator.mul(a, b);
                compareMatrix.readMatrix("mul.mat");
                res.compareMat(compareMatrix);
                res.clear();
                System.out.println(i++);

                res = MatrixCalculator.matrixMul(a, b);
                compareMatrix.readMatrix("matrixMul.mat");
                res.compareMat(compareMatrix);
                res.clear();
                System.out.println(i++);

                res = MatrixCalculator.max(a, b);
                compareMatrix.readMatrix("max.mat");
                res.compareMat(compareMatrix);
                res.clear();
                System.out.println(i++);

                res = MatrixCalculator.min(a, b);
                compareMatrix.readMatrix("min.mat");
                res.compareMat(compareMatrix);
                res.clear();
                System.out.println(i++);

                lres = MatrixCalculator.bigger(a, b);
                compareLogicMatrix.readMatrix("bigger.mat");
                lres.compareLogicMat(compareLogicMatrix);
                lres.clear();
                System.out.println(i++);

                lres = MatrixCalculator.smaller(a, b);
                compareLogicMatrix.readMatrix("smaller.mat");
                lres.compareLogicMat(compareLogicMatrix);
                lres.clear();
                System.out.println(i++);

                lres = MatrixCalculator.equal(a, b);
                compareLogicMatrix.readMatrix("equal.mat");
                lres.compareLogicMat(compareLogicMatrix);
                lres.clear();
                System.out.println(i++);

                lres = MatrixCalculator.between(a, lowerBound, upperBound);
                compareLogicMatrix.readMatrix("between.mat");
                lres.compareLogicMat(compareLogicMatrix);
                lres.clear();
                System.out.println(i++);

                lres = MatrixCalculator.and(la, lb);
                compareLogicMatrix.readMatrix("matAnd.mat");
                lres.compareLogicMat(compareLogicMatrix);
                lres.clear();
                System.out.println(i++);

                lres = MatrixCalculator.or(la, lb);
                compareLogicMatrix.readMatrix("matOr.mat");
                lres.compareLogicMat(compareLogicMatrix);
                lres.clear();
                System.out.println(i++);

                lres = MatrixCalculator.not(la);
                compareLogicMatrix.readMatrix("matNot.mat");
                lres.compareLogicMat(compareLogicMatrix);
                lres.clear();
                System.out.println(i++);

                res = MatrixCalculator.condition(la, a, b);
                compareMatrix.readMatrix("condition.mat");
                res.compareMat(compareMatrix);
                res.clear();
                System.out.println(i++);

                res = MatrixCalculator.rank(a);
                compareMatrix.readMatrix("rank.mat");
                res.compareMat(compareMatrix);
                res.clear();
                System.out.println(i++);

                res = MatrixCalculator.round(a);
                compareMatrix.readMatrix("round.mat");
                res.compareMat(compareMatrix);
                res.clear();
                System.out.println(i++);

                res = MatrixCalculator.floor(a);
                compareMatrix.readMatrix("floor.mat");
                res.compareMat(compareMatrix);
                res.clear();
                System.out.println(i++);

                res = MatrixCalculator.abs(a);
                compareMatrix.readMatrix("abs.mat");
                res.compareMat(compareMatrix);
                res.clear();
                System.out.println(i++);

                res = MatrixCalculator.minus(a);
                compareMatrix.readMatrix("minus.mat");
                res.compareMat(compareMatrix);
                res.clear();
                System.out.println(i++);

                res = MatrixCalculator.sqrt(a);
                compareMatrix.readMatrix("sqrt.mat");
                res.compareMat(compareMatrix);
                res.clear();
                System.out.println(i++);

                res = MatrixCalculator.log(a);
                compareMatrix.readMatrix("log.mat");
                res.compareMat(compareMatrix);
                res.clear();
                System.out.println(i++);

                res = MatrixCalculator.exp(a);
                compareMatrix.readMatrix("exp.mat");
                res.compareMat(compareMatrix);
                res.clear();
                System.out.println(i++);

                res = MatrixCalculator.sign(a);
                compareMatrix.readMatrix("sign.mat");
                res.compareMat(compareMatrix);
                res.clear();
                System.out.println(i++);

                res = MatrixCalculator.inverse(a);
                compareMatrix.readMatrix("inverse.mat");
                res.compareMat(compareMatrix);
                res.clear();
                System.out.println(i++);

                res = MatrixCalculator.signedpow(a, 2);
                compareMatrix.readMatrix("signedpow.mat");
                res.compareMat(compareMatrix);
                res.clear();
                System.out.println(i++);

                res = MatrixCalculator.shift(a, stepSize);
                compareMatrix.readMatrix("shift.mat");
                res.compareMat(compareMatrix);
                res.clear();
                System.out.println(i++);

                res = MatrixCalculator.delay(a, stepSize);
                compareMatrix.readMatrix("delay.mat");
                res.compareMat(compareMatrix);
                res.clear();
                System.out.println(i++);

                res = MatrixCalculator.delta(a, stepSize);
                compareMatrix.readMatrix("delta.mat");
                res.compareMat(compareMatrix);
                res.clear();
                System.out.println(i++);

                res = MatrixCalculator.ratio(a, stepSize);
                compareMatrix.readMatrix("ratio.mat");
                res.compareMat(compareMatrix);
                res.clear();
                System.out.println(i++);

                res = MatrixCalculator.sum(a, stepSize);
                compareMatrix.readMatrix("sum.mat");
                res.compareMat(compareMatrix);
                res.clear();
                System.out.println(i++);

                res = MatrixCalculator.product(a, stepSize);
                compareMatrix.readMatrix("product.mat");
                res.compareMat(compareMatrix);
                res.clear();
                System.out.println(i++);

                res = MatrixCalculator.tsMax(a, stepSize);
                compareMatrix.readMatrix("tsmax.mat");
                res.compareMat(compareMatrix);
                res.clear();
                System.out.println(i++);

                res = MatrixCalculator.tsMin(a, stepSize);
                compareMatrix.readMatrix("tsmin.mat");
                res.compareMat(compareMatrix);
                res.clear();
                System.out.println(i++);

                res = MatrixCalculator.tsArgmax(a, stepSize);
                compareMatrix.readMatrix("tsargmax.mat");
                res.compareMat(compareMatrix);
                res.clear();
                System.out.println(i++);

                res = MatrixCalculator.tsArgmin(a, stepSize);
                compareMatrix.readMatrix("tsargmin.mat");
                res.compareMat(compareMatrix);
                res.clear();
                System.out.println(i++);

                res = MatrixCalculator.tsRank(a, stepSize);
                compareMatrix.readMatrix("tsrank.mat");
                res.compareMat(compareMatrix);
                res.clear();
                System.out.println(i++);

                res = MatrixCalculator.tsMean(a, stepSize);
                compareMatrix.readMatrix("tsmean.mat");
                res.compareMat(compareMatrix);
                res.clear();
                System.out.println(i++);

                res = MatrixCalculator.tsStd(a, stepSize);
                compareMatrix.readMatrix("tsstd.mat");
                res.compareMat(compareMatrix);
                res.clear();
                System.out.println(i++);

                res = MatrixCalculator.tsSkewness(a, stepSize);
                compareMatrix.readMatrix("tsskewness.mat");
                res.compareMat(compareMatrix);
                res.clear();
                System.out.println(i++);

                res = MatrixCalculator.tsKurtosis(a, stepSize);
                compareMatrix.readMatrix("tskurtosis.mat");
                res.compareMat(compareMatrix);
                res.clear();
                System.out.println(i++);

                res = MatrixCalculator.tsCov(a, b, stepSize);
                compareMatrix.readMatrix("tscov.mat");
                res.compareMat(compareMatrix);
                res.clear();
                System.out.println(i++);

                res = MatrixCalculator.tsCorr(a, b, stepSize);
                compareMatrix.readMatrix("tscorr.mat");
                res.compareMat(compareMatrix);
                res.clear();
                System.out.println(i++);

                res = MatrixCalculator.tsCountTrue(la, stepSize);
                compareMatrix.readMatrix("tscounttrue.mat");
                res.compareMat(compareMatrix);
                res.clear();
                System.out.println(i++);

                res = MatrixCalculator.tsCountConsecutiveTrue(la, stepSize);
                compareMatrix.readMatrix("tscountconsecutivetrue.mat");
                res.compareMat(compareMatrix);
                res.clear();
                System.out.println(i++);

                res = MatrixCalculator.decayLinear(a, stepSize);
                compareMatrix.readMatrix("decaylinear.mat");
                res.compareMat(compareMatrix);
                res.clear();
                System.out.println(i++);

                res = MatrixCalculator.decayExponential(a, stepSize);
                compareMatrix.readMatrix("decayexponential.mat");
                res.compareMat(compareMatrix);
                res.clear();

                res = MatrixCalculator.neutralize(a);
                compareMatrix.readMatrix("neutralize.mat");
                res.compareMat(compareMatrix);
                res.clear();

                res = MatrixCalculator.unify(a);
                compareMatrix.readMatrix("unify.mat");
                res.compareMat(compareMatrix);
                res.clear();

                res = MatrixCalculator.evalValidPct(a);
                compareMatrix.readMatrix("evalvalidpct.mat");
                res.compareMat(compareMatrix);
                res.clear();

                res = MatrixCalculator.evalAbsSum(a);
                compareMatrix.readMatrix("evalabssum.mat");
                res.compareMat(compareMatrix);
                res.clear();

                res = MatrixCalculator.evalMean(a);
                compareMatrix.readMatrix("evalmean.mat");
                res.compareMat(compareMatrix);
                res.clear();

                res = MatrixCalculator.evalVariance(a);
                compareMatrix.readMatrix("evalvariance.mat");
                res.compareMat(compareMatrix);
                res.clear();

                res = MatrixCalculator.evalInnerProduction(a, b);
                compareMatrix.readMatrix("evalinnerproduction.mat");
                res.compareMat(compareMatrix);
                res.clear();

                res = MatrixCalculator.evalCovariance(a, b);
                compareMatrix.readMatrix("evalcovariance.mat");
                res.compareMat(compareMatrix);
                res.clear();

                res = MatrixCalculator.evalCorrelation(a, b);
                compareMatrix.readMatrix("evalcorrelation.mat");
                res.compareMat(compareMatrix);
                res.clear();

                res = MatrixCalculator.evalBeta(a, b);
                compareMatrix.readMatrix("evalbeta.mat");
                res.compareMat(compareMatrix);
                res.clear();

                res = MatrixCalculator.diag(a);
                compareMatrix.readMatrix("diag.mat");
                res.compareMat(compareMatrix);
                res.clear();

                res = MatrixCalculator.cumSum(a);
                compareMatrix.readMatrix("cumsum.mat");
                res.compareMat(compareMatrix);
                res.clear();

        }
}
