/*************************************************************************
	> File Name: MatrixFactory.java
	> Author: Lianke QIN
	> Mail: lianke.qin@gmail.com 
	> Created Time: Mon Aug  5 16:40:20 2019
 ************************************************************************/

package blas.java;

import blas.java.*;

public class MatrixFactory {

	static {
		System.loadLibrary("MatrixFactory");
	}

	public MatrixFactory() {
	}

	public static void main(String[] args) {
	}

	public static Matrix getInstanceOfEmptyMatrix() {
		Matrix res = new Matrix();
		return res;
	}

	public static Matrix getInstanceOfNaNMatrix(int size) {
		Matrix res = new Matrix();
		res.setPtr(getInstanceOfNaNMatrixNative(size));
		return res;
	}

	public static Matrix getInstanceOfRandomMatrix(int n, int m, int min, int max) {
		Matrix res = new Matrix();
		res.setPtr(getInstanceOfRandomMatrixNative(n, m, min, max));
		return res;
	}

	public static Matrix getInstanceOfRandomMatrixWithAbnormalValue(int n, int m, int min, int max, int frequency) {
		Matrix res = new Matrix();
		res.setPtr(getInstanceOfRandomMatrixWithAbnormalValueNative(n, m, min, max, frequency));
		return res;
	}

	public static LogicMatrix getInstanceOfRandomLogicMatrix(int n, int m) {
		LogicMatrix res = new LogicMatrix();
		res.setPtr(getInstanceOfRandomLogicMatrixNative(n, m));
		return res;
	}

	public static native long getInstanceOfEmptyMatrixNative();

	public static native long getInstanceOfNaNMatrixNative(int size);

	public static native long getInstanceOfRandomMatrixNative(int n, int m, int min, int max);

	public static native long getInstanceOfRandomLogicMatrixNative(int n, int m);

	public static native long getInstanceOfRandomMatrixWithAbnormalValueNative(int n, int m, int min, int max,
			int frequency);
}
