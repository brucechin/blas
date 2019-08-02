package blas.java;
import java.lang.Thread;
import java.io.*;
public class LogicMatrix implements AutoCloseable{
	static{
		System.loadLibrary("Matrix");
	}
	public long ptr;//pointing to the memory allocated to this Matrix
	//public int nrow;
	//public int ncol;
	public LogicMatrix(){
		ptr = ccLogicMatrixNative(0, 0);
	}

	public LogicMatrix(int n, int m){
		ptr = ccLogicMatrixNative(n, m);
	}

	public LogicMatrix(String filename){
		ptr = ccLogicMatrixNative(0, 0);
		readMatrix(filename);
	}

	public LogicMatrix(LogicMatrix other){
	
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
		System.out.println("closing ");
		clearNative(ptr);
	}

	public long getPtr(){
		return ptr;
	}

	public void setPtr(long p){
		ptr = p;
	}

	public int getNRow(){
		int res = getNRowNative(ptr);
		return res;
	}

	public int getNCol(){
		int res = getNColNative(ptr);
		return res;
	}

	public void setElement(int i, int j, boolean x){
		setElementNative(i, j, x, ptr);
	}
	
	public boolean getElement(int i, int j){
		return getElementNative(i, j, ptr);
	}
	
	public LogicMatrix getRowVector(int i){
		LogicMatrix res = new LogicMatrix(0, 0);
		long resPtr = getRowVectorNative(i, ptr);
		res.setPtr(resPtr);
		return res;
	}

	public LogicMatrix getColVector(int j){
		LogicMatrix res = new LogicMatrix(0, 0);
		long resPtr = getColVectorNative(j, ptr);
		res.setPtr(resPtr);
		return res;
	}

	public void print(){printNative(ptr);}

	public void saveMatrix(String filename){saveMatrixNative(filename, ptr);}
	
	public void readMatrix(String filename){readMatrixNative(filename, ptr);}
	
	public native void setElementNative(int i, int j, boolean x, long ptr);

	public native boolean getElementNative(int i, int j, long ptr);

	public native long getRowVectorNative(int i, long ptr);

	public native long getColVectorNative(int j, long ptr);
	
	public native int getNRowNative(long ptr);

	public native int getNColNative(long ptr);

	public native void printNative(long ptr);

	public native void clearNative(long ptr);
	
	public native void saveMatrixNative(String filename, long ptr);

	public native void readMatrixNative(String filename, long ptr);

	public native long ccLogicMatrixNative(int n, int m);
	
	public static void main(String[] args){
	
		Matrix a = new Matrix(5, 5);
		try(Matrix b = new Matrix(10, 10);
			Matrix c = new Matrix(15, 15)){

			a.print();
			b.setElement(0, 0, 0);
			b.saveMatrix("test.mat");
			a.readMatrix("test.mat");
			a.print();
			double t = a.getElement(1, 1);
			System.out.println(t);
			
			Matrix arow = a.getRowVector(1);
			System.out.println("\n");
			arow.print();
			Matrix acol = a.getColVector(2);
			acol.print();
			
			}
		catch(IOException e){
			System.out.println("exception");
		}
	
	
	}
}
