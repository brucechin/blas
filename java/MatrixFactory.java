/*************************************************************************
	> File Name: MatrixFactory.java
	> Author: Lianke QIN
	> Mail: lianke.qin@gmail.com 
	> Created Time: Mon Aug  5 16:40:20 2019
 ************************************************************************/

package blas.java;
import blas.java.*;

public class MatrixFactory{

	static{
		System.loadLibrary("MatrixFactory");
	}

	public MatrixFactory(){}

	public static void main(String[] args){}

	public static Matrix getInstanceOfEmptyMatrix(){
		Matrix res = new Matrix();
		return res;
	}


	public static native long getInstanceOfEmptyMatrixNative(int size);
	public static native long getInstanceOfNaNMatrixNative(int size);
	public static native long getInstanceOfDiagMatrixNative(int size);
	public static native long getInstanceOfRandomMatrixNative(int n, int m, int min, int max);
	public static native long getInstanceOfRandomLogicMatrixNative(int n, int m);
	
}
