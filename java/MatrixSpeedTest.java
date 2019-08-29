package blas.java;

import blas.java.*;

public class MatrixSpeedTest {
    public static void main(String[] args) {

        int rangeMin = -10000;
        int rangeMax = 10000;
        int matSize = 3000;
        int lowerBound = -100;
        int upperBound = 100;
        int stepSize = 50; // used for shift/delta/delay/ratio/sum/product and timeseries functions
        String matrix_a = "a.mat";
        String matrix_b = "b.mat";
        String logicMatrix_a = "la.mat";
        String logicMatrix_b = "lb.mat";
        /*
         * Matrix a = new Matrix(matrix_a); Matrix b = new Matrix(matrix_b); LogicMatrix
         * la = new LogicMatrix(logicMatrix_a); LogicMatrix lb = new
         * LogicMatrix(logicMatrix_b);
         */

        Matrix a = MatrixFactory.getInstanceOfRandomMatrix(matSize, matSize, rangeMin, rangeMax);
        Matrix b = MatrixFactory.getInstanceOfRandomMatrix(matSize, matSize, rangeMin, rangeMax);
        LogicMatrix la = MatrixFactory.getInstanceOfRandomLogicMatrix(matSize, matSize);
        LogicMatrix lb = MatrixFactory.getInstanceOfRandomLogicMatrix(matSize, matSize);
        Matrix res;
        LogicMatrix lres;
        long startTime;
        long endTime;
        int times = 10;
        /*
         * System.out.println("Det(a) : " + MatrixCalculator.Det(a.value, a.getNRow() -
         * 1));
         * 
         * System.out.println("treat(a) : " + MatrixCalculator.treat(a)); res =
         * MatrixCalculator.Normalize(a, 2.0, 10.0, 1000.0);
         * 
         */

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.add(a, b);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("add : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.sub(a, b);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("sub : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.div(a, b);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("div : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.mul(a, b);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("mul : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.matrixMul(a, b);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("matrixMul : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.max(a, b);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("max : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.min(a, b);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("min : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            lres = MatrixCalculator.bigger(a, b);
            lres.clear();
        }

        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("bigger : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            lres = MatrixCalculator.smaller(a, b);
            lres.clear();
        }

        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("smaller : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            lres = MatrixCalculator.equal(a, b);
            lres.clear();
        }

        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("equal : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            lres = MatrixCalculator.between(a, lowerBound, upperBound);
            lres.clear();
        }

        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("between : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            lres = MatrixCalculator.and(la, lb);
            lres.clear();
        }

        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("and : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            lres = MatrixCalculator.or(la, lb);
            lres.clear();
        }

        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("or : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            lres = MatrixCalculator.not(la);
            lres.clear();
        }

        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("not : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.condition(la, a, b);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("condition : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        // for (int i = 0; i < times; i++) {
        // res = MatrixCalculator.rank(a);
        // res.clear();
        // }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("rank : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.round(a);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("round : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.floor(a);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("floor : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.abs(a);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("abs : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.minus(a);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("minus : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.sqrt(a);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("sqrt : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.log(a);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("log : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.exp(a);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("exp : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.sign(a);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("sign : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.inverse(a);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("inverse : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.signedpow(a, 2);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("signedpow : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.shift(a, stepSize);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("shift : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.delay(a, stepSize);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("delay : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.delta(a, stepSize);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("delta : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.ratio(a, stepSize);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("ratio : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.sum(a, stepSize);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("sum : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.product(a, stepSize);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("product : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.tsMax(a, stepSize);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("tsMax : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.tsMin(a, stepSize);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("tsMin : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.tsArgmax(a, stepSize);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("tsArgmax : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.tsArgmin(a, stepSize);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("tsArgmin : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.tsRank(a, stepSize);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("tsRank : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.tsMean(a, stepSize);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("tsMean : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.tsStd(a, stepSize);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("tsStd : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.tsSkewness(a, stepSize);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("tsSkewness : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.tsKurtosis(a, stepSize);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("tsKurtosis : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.tsCov(a, b, stepSize);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("tsCov : " + (endTime - startTime) / 1000000 + "ms");

        // res = MatrixCalculator.tscorr(a, b, stepSize);
        //
        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.tsCountTrue(la, stepSize);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("tsCountTrue : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.tsCountConsecutiveTrue(la, stepSize);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("tsCountConsecutiveTrue : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.decayLinear(a, stepSize);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("decayLinear : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.decayExponential(a, stepSize);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("decayExponential : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.neutralize(a);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("neutralize : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.unify(a);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("unify : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.evalValidPct(a);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("evalValidPct : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.evalAbsSum(a);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("evalAbsSum : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.evalMean(a);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("evalMean : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.evalVariance(a);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("evalVariance : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.evalInnerProduction(a, b);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("evalInnerProduction : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.evalCovariance(a, b);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("evalCovariance : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.evalCorrelation(a, b);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("evalCorrelation : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.evalBeta(a, b);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("evalBeta : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.diag(a);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("diag : " + (endTime - startTime) / 1000000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.cumSum(a);
            res.clear();
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("cumSum : " + (endTime - startTime) / 1000000 + "ms");

    }
}
