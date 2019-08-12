package blas.java;

import blas.java.*;
import java.lang.*;

public class Consumer implements Runnable {
    private MatrixPool _pool;
    private int _ncol;
    private int _nrow;

    public Consumer(MatrixPool pool, int nrow, int ncol) {
        _pool = pool;
        _nrow = nrow;
        _ncol = ncol;
    }

    public void run() {
        while (true) {

            try {
                Thread.sleep(501 + (int) (Math.random() * 500));
                PooledMatrix mat = _pool.get(_nrow, _ncol);
                mat.readMatrix("a.mat");
                Thread.sleep(2001 + (int) (Math.random() * 2000));
                _pool.release(mat);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }

        }
    }
}