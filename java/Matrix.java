package blas.java;
import java.lang.Thread;
public class Matrix{
	static{
		System.loadLibrary("Matrix");
	}
	public long ptr;//pointing to the memory allocated to this Matrix
	public int nrow;
	public int ncol;
	public Matrix(){
		ptr = ccMatrix(0, 0);
	}

	public Matrix(int n, int m){
		ptr = ccMatrix(n, m);
		System.out.println("construction done");
	}

//	public Matrix(String filename){
//		ptr = readMatrix(filename);
//	}

	//public Matrix(long p); 
	public void finalize(){
		System.out.println("deconstruct");
	}
	public native void clear();

	public native int getNRow();

	public native int getNCol();

	//public void setElement(int i, int j, double x);
	//
	//public double getElement(int i, int j);
	
	public native void print();

	//public native void saveMatrix(String filename);

	//public native long readMatrix(String filename);

	public native long ccMatrix(int n, int m);
	
	public static void main(String[] args){
		
		{
			Matrix mat = new Matrix(10000, 10000);
		}
	
		{
			Matrix mat = new Matrix(10000, 10000);
		}
		{
			Matrix mat = new Matrix(10000, 10000);
		}
		{
			Matrix mat = new Matrix(10000, 10000);
		}
		{
			Matrix mat = new Matrix(10000, 10000);
		}
		{
			Matrix mat = new Matrix(10000, 10000);
		}
		{
			Matrix mat = new Matrix(10000, 10000);
		}
		{
			Matrix mat = new Matrix(10000, 10000);
		}
		{
			Matrix mat = new Matrix(10000, 10000);
		}
		{
			Matrix mat = new Matrix(10000, 10000);
		}
		{
			Matrix mat = new Matrix(10000, 10000);
		}
		{
			Matrix mat = new Matrix(10000, 10000);
		}
		{
			Matrix mat = new Matrix(10000, 10000);
		}
		{
			Matrix mat = new Matrix(10000, 10000);
		}
		{
			Matrix mat = new Matrix(10000, 10000);
		}
		{
			Matrix mat = new Matrix(10000, 10000);
		}
		{
			Matrix mat = new Matrix(10000, 10000);
		}
		{
			Matrix mat = new Matrix(10000, 10000);
		}
		{
			Matrix mat = new Matrix(10000, 10000);
		}
		{
			Matrix mat = new Matrix(10000, 10000);
		}
		{
			Matrix mat = new Matrix(10000, 10000);
		}
		{
			Matrix mat = new Matrix(10000, 10000);
		}
		{
			Matrix mat = new Matrix(10000, 10000);
		}
		{
			Matrix mat = new Matrix(10000, 10000);
		}
		{
			Matrix mat = new Matrix(10000, 10000);
		}
		{
			Matrix mat = new Matrix(10000, 10000);
		}
		{
			Matrix mat = new Matrix(10000, 10000);
		}
		{
			Matrix mat = new Matrix(10000, 10000);
		}
		{
			Matrix mat = new Matrix(10000, 10000);
		}
		{
			Matrix mat = new Matrix(10000, 10000);
		}
		{
			Matrix mat = new Matrix(10000, 10000);
		}
		{
			Matrix mat = new Matrix(10000, 10000);
		}
		{
			Matrix mat = new Matrix(10000, 10000);
		}
		{
			Matrix mat = new Matrix(10000, 10000);
		}
		{
			Matrix mat = new Matrix(10000, 10000);
		}
		{
			Matrix mat = new Matrix(10000, 10000);
		}
		{
			Matrix mat = new Matrix(10000, 10000);
		}
	}
}
