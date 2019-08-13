package blas.java;

import blas.java.*;
import java.lang.*;

public class Consumer implements Runnable {
    private MatrixPoolManager _manager;
    private int _ncol;
    private int _nrow;

    public Consumer(int nrow, int ncol) {
        _manager = MatrixPoolManager.getMatrixPoolManager();
        _nrow = nrow;
        _ncol = ncol;
    }

    public void run() {
        while (true) {

            try {
                Thread.sleep(501 + (int) (Math.random() * 500));
                Matrix mat = _manager.get(_nrow, _ncol);
                // mat.readMatrix("a.mat");
                Thread.sleep(2001 + (int) (Math.random() * 2000));
                _manager.release(mat);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }

        }
    }
}