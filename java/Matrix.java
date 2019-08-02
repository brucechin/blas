package blas.java;
import java.lang.Thread;
import java.io.*;
public class Matrix implements AutoCloseable{
	static{
		System.loadLibrary("Matrix");
	}
	public long ptr;//pointing to the memory allocated to this Matrix
	public int nrow;
	public int ncol;
	public Matrix(){
		ptr = ccMatrixNative(0, 0);
	}

	public Matrix(int n, int m){
		nrow = n;
		ncol = m;
		ptr = ccMatrixNative(n, m);
	}

	public Matrix(String filename){
		ptr = ccMatrixNative(0, 0);
		readMatrix(filename);
	}

	public Matrix(Matrix other){
	
		//TODO : deep copy
	
	}

	public void finalize(){
		System.out.println("deconstruct");
	}

	public void close() throws IOException{
		//REMINDER : all instances initialized in try block will be auto-released when try block ends;
		//
		//example : 
		//
		//
		//	Matrix a = new Matrix(10, 10);
		//	try(Matrix b = new Matrix(20, 20)){
		//		do something
		//	}
		//
		//Matrix b will be deconstructed here, but Matrix a will not.
		System.out.println("closing " + nrow);
		nrow = 0;
		ncol = 0;
		clearNative(ptr);
	}

	public long getPtr(){
		return ptr;
	}

	public int getNRow(){
		int res = getNRowNative(ptr);
		return res;
	}

	public int getNCol(){
		int res = getNColNative(ptr);
		return res;
	}

	public void setElement(int i, int j, double x){
		setElementNative(i, j, x, ptr);
	}
	
	public double getElement(int i, int j){
		return getElementNative(i, j, ptr);
	}
	
	public Matrix getRowVector(int i){
		Matrix res = new Matrix(0, 0);
		long resPtr = res.getPtr();
		getRowVectorNative(i, ptr, resPtr);
		return res;
	}

	public Matrix getColVector(int j){
		Matrix res = new Matrix(0, 0);
		long resPtr = res.getPtr();
		getColVectorNative(j, ptr, resPtr);
		return res;
	}

	public void print(){printNative(ptr);}

	public void saveMatrix(String filename){saveMatrixNative(filename, ptr);}
	
	public void readMatrix(String filename){readMatrixNative(filename, ptr);}
	
	public native double setElementNative(int i, int j, double x, long ptr);

	public native double getElementNative(int i, int j, long ptr);

	public native void getRowVectorNative(int i, long ptr, long resPtr);

	public native void getColVectorNative(int j, long ptr, long resPtr);
	
	public native int getNRowNative(long ptr);

	public native int getNColNative(long ptr);

	public native void printNative(long ptr);

	public native void clearNative(long ptr);
	
	public native void saveMatrixNative(String filename, long ptr);

	public native void readMatrixNative(String filename, long ptr);

	public native long ccMatrixNative(int n, int m);
	
	public static void main(String[] args){
	
		Matrix a = new Matrix(5, 5);
		try(Matrix b = new Matrix(10, 10);
			Matrix c = new Matrix(15, 15)){

			a.print();
				b.saveMatrix("test.mat");
				a.readMatrix("test.mat");
				a.print();
				System.out.println("running");
			}
		catch(IOException e){
			System.out.println("exception");
		}
	
	
	}
}
