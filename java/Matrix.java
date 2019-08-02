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
		ptr = ccMatrix(0, 0);
	}

	public Matrix(int n, int m){
		nrow = n;
		ncol = m;
		ptr = ccMatrix(n, m);
	}

//	public Matrix(String filename){
//		ptr = readMatrix(filename);
//	}

	//public Matrix(long p); 
	public void finalize(){
		System.out.println("deconstruct");
	}

	public void close() throws IOException{
	
		System.out.println("closing" + nrow);
		nrow = 0;
		ncol = 0;
		clear();
	}
	public native void clear();

	public int getNRow(){return this.nrow;}

	public int getNCol(){return this.ncol;}

	//public void setElement(int i, int j, double x);
	
	//public double getElement(int i, int j);
	
	//public native void print();

	public native void saveMatrix(String filename);

	public native long readMatrix(String filename);

	public native long ccMatrix(int n, int m);
	
	public static void main(String[] args){
	
		Matrix a = new Matrix(100, 100);
		try(Matrix b = new Matrix(200, 200);
			Matrix c = new Matrix(300, 300)){
				System.out.println("running");
			}
		catch(IOException e){
			System.out.println("exception");
		}
	
	
	}
}
