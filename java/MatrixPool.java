package blas.java;

import java.util.*;
import blas.java.*;
import java.lang.*;

public class MatrixPool {
    private int _size;
    private int _nrow;
    private int _ncol;
    private boolean _isOn;
    private LinkedList<PooledMatrix> _freeMatrixList;
    private LinkedList<PooledMatrix> _usedMatrixList;

    public static void main(String[] args) throws InterruptedException {
        // MatrixPool pool = new MatrixPool(50, 2000, 2000);
        // for (int i = 0; i < 200; i++) {
        // Thread.sleep(500);
        // PooledMatrix mat1 = pool.get(2000, 2000);
        // mat1.readMatrix("a.mat");
        // PooledMatrix mat2 = pool.get(2000, 2000);
        // mat2.readMatrix("b.mat");
        // pool.release(mat1);
        // pool.release(mat2);
        // System.out.println("releasing");
        // if (i == 100) {
        // pool.shutDown();
        // }

        // }
    }

    public MatrixPool(int poolsize, int nrow, int ncol) {
        _size = poolsize;
        _freeMatrixList = new LinkedList<PooledMatrix>();
        _usedMatrixList = new LinkedList<PooledMatrix>();
        _ncol = ncol;
        _nrow = nrow;
        _isOn = true;
        for (int i = 0; i < _size; i++) {
            PooledMatrix mat = new PooledMatrix(nrow, ncol);
            _freeMatrixList.add(mat);
        }

    }

    public synchronized PooledMatrix get(int nrow, int ncol) {
        PooledMatrix res;

        if (_isOn) {
            if (_freeMatrixList.size() > 0 && _ncol == ncol && _nrow == nrow) {
                res = _freeMatrixList.getFirst();
                _freeMatrixList.remove();
                _usedMatrixList.add(res);
            } else {
                System.out.println("allocating new space");
                res = new PooledMatrix(nrow, ncol);
            }
            System.out.println(_freeMatrixList.size());
            return res;
        } else {
            res = new PooledMatrix();
            return res;
        }

    }

    public synchronized boolean release(PooledMatrix mat) {
        if (_usedMatrixList.contains(mat)) {
            _usedMatrixList.remove(mat);
            _freeMatrixList.add(mat);
            return true;
        } else if (!_freeMatrixList.contains(mat) && mat.getNCol() == _ncol && mat.getNRow() == _nrow) {
            _freeMatrixList.add(mat);
            _size++;
            return true;
        }

        return false;
    }

    // no need to implement this API for now
    // public synchronized boolean setPoolSize(int newSize) {
    //
    // return true;
    // }

    public synchronized boolean shutDown() throws InterruptedException {
        // TODO implement it later
        _isOn = false;
        System.out.println("Start shutting down the matrix pool");
        while (_usedMatrixList.size() != 0) {
            Thread.sleep(1000);
        }
        System.out.println("Matrix pool shuts down completely");
        return true;
    }

}
