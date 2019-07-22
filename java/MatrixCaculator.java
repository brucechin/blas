package blas.java;
public class MatrixCalculator{
	
	static{
		System.loadLibrary("MatrixCalculator");
	}

	long nativeMatrixCalculator;
	
	public MatrixCalculator(){
		nativeMatrixCalculator = createMatrixCalculator();
	}

	private native long createMatrixCalculator();



	public static native long add(long mat1, long mat2);

	public static native long sub(long mat1, long mat2);



}
