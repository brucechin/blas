package blas.java;

import blas.java.*;

public class MatrixSpeedTest {
    public static void main(String[] args) {

        int rangeMin = -10000;
        int rangeMax = 10000;
        int matSize = 1000;
        int lowerBound = -100;
        int upperBound = 100;
        int stepSize = 50; // used for shift/delta/delay/ratio/sum/product and timeseries functions
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
        long startTime;
        long endTime;
        int times = 100;
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
            MatrixCalculator.smoothByDecayLinear(a, 5);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("smooth : " + (endTime - startTime) / 1000 + "ms");

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.add(a, b);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("add : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.sub(a, b);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("sub : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.div(a, b);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("div : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.mul(a, b);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("mul : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.matrixMul(a, b);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("matrixMul : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.max(a, b);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("max : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.min(a, b);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("min : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            lres = MatrixCalculator.bigger(a, b);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("bigger : " + (endTime - startTime) / 1000 + "ms");
        lres.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            lres = MatrixCalculator.smaller(a, b);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("smaller : " + (endTime - startTime) / 1000 + "ms");
        lres.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            lres = MatrixCalculator.equal(a, b);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("equal : " + (endTime - startTime) / 1000 + "ms");
        lres.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            lres = MatrixCalculator.between(a, lowerBound, upperBound);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("between : " + (endTime - startTime) / 1000 + "ms");
        lres.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            lres = MatrixCalculator.and(la, lb);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("and : " + (endTime - startTime) / 1000 + "ms");
        lres.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            lres = MatrixCalculator.or(la, lb);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("or : " + (endTime - startTime) / 1000 + "ms");
        lres.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            lres = MatrixCalculator.not(la);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("not : " + (endTime - startTime) / 1000 + "ms");
        lres.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.condition(la, a, b);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("condition : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.rank(a);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("rank : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.round(a);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("round : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.floor(a);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("floor : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.abs(a);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("abs : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.minus(a);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("minus : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.sqrt(a);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("sqrt : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.log(a);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("log : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.exp(a);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("exp : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.sign(a);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("sign : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.inverse(a);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("inverse : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.signedpow(a, 2);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("signedpow : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.shift(a, stepSize);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("shift : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.delay(a, stepSize);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("delay : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.delta(a, stepSize);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("delta : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.ratio(a, stepSize);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("ratio : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.sum(a, stepSize);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("sum : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.product(a, stepSize);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("product : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.tsMax(a, stepSize);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("tsMax : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.tsMin(a, stepSize);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("tsMin : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.tsArgmax(a, stepSize);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("tsArgmax : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.tsArgmin(a, stepSize);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("tsArgmin : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.tsRank(a, stepSize);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("tsRank : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.tsMean(a, stepSize);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("tsMean : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.tsStd(a, stepSize);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("tsStd : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.tsSkewness(a, stepSize);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("tsSkewness : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.tsKurtosis(a, stepSize);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("tsKurtosis : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.tsCov(a, b, stepSize);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("tsCov : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        // res = MatrixCalculator.tscorr(a, b, stepSize);
        // res.clear();
        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.tsCountTrue(la, stepSize);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("tsCountTrue : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.tsCountConsecutiveTrue(la, stepSize);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("tsCountConsecutiveTrue : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.decayLinear(a, stepSize);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("decayLinear : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.decayExponential(a, stepSize);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("decayExponential : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.neutralize(a);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("neutralize : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.unify(a);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("unify : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.evalValidPct(a);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("evalValidPct : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.evalAbsSum(a);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("evalAbsSum : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.evalMean(a);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("evalMean : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.evalVariance(a);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("evalVariance : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.evalInnerProduction(a, b);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("evalInnerProduction : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.evalCovariance(a, b);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("evalCovariance : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.evalCorrelation(a, b);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("evalCorrelation : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.evalBeta(a, b);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("evalBeta : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.diag(a);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("diag : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

        startTime = System.nanoTime(); // 获取开始时间
        for (int i = 0; i < times; i++) {
            res = MatrixCalculator.cumSum(a);
        }
        endTime = System.nanoTime(); // 获取结束时间
        System.out.println("cumSum : " + (endTime - startTime) / 1000 + "ms");
        res.clear();

    }
}
