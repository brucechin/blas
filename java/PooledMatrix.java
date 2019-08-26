package blas.java;

import blas.java.*;
import java.util.concurrent.*;

public class PooledMatrix extends Matrix {

    public PooledMatrix() {
        super();
    }

    public PooledMatrix(int n, int m) {
        super(n, m);
    }

    public static void main(String[] args) {
        PooledMatrix test = new PooledMatrix(10, 10);
        for (int i = 0; i < test.getNRow(); i++) {
            for (int j = 0; j < test.getNCol(); j++) {
                test.setElement(i, j, j);
            }
        }
        test.print();
    }
}