package blas.java;

import java.util.*;
import java.util.concurrent.locks.*;

import blas.java.*;
import java.lang.*;

public class MatrixPool {
    private int _minPoolSize;
    private int _maxPoolSize;
    private int _nrow;
    private int _ncol;
    private boolean _isOn;
    private LinkedList<Matrix> _freeMatrixList;
    private LinkedList<Matrix> _usedMatrixList;
    final Lock lock = new ReentrantLock();

    public static void main(String[] args) throws InterruptedException {
        int nrow = 2000;
        int ncol = 2000;
        int minsize = 10;
        int maxsize = 30;
        MatrixPool pool = new MatrixPool(minsize, maxsize, nrow, ncol);
        for (int i = 0; i < maxsize; i++) {
            new Thread(new Consumer(pool, nrow, ncol), "Thread " + String.valueOf(i)).start();
        }

    }

    public MatrixPool(int minsize, int maxsize, int nrow, int ncol) {
        _minPoolSize = minsize;
        _maxPoolSize = maxsize;
        _freeMatrixList = new LinkedList<Matrix>();
        _usedMatrixList = new LinkedList<Matrix>();
        _ncol = ncol;
        _nrow = nrow;
        _isOn = true;
        for (int i = 0; i < _minPoolSize; i++) {
            PooledMatrix mat = new PooledMatrix(nrow, ncol);
            _freeMatrixList.add(mat);
        }

    }

    public Matrix get(int nrow, int ncol) {
        Matrix res = null;
        lock.lock();
        try {
            if (_isOn && _ncol == ncol && _nrow == nrow) {
                if (_freeMatrixList.size() > 0) {
                    res = _freeMatrixList.getFirst();
                    _freeMatrixList.remove();
                    _usedMatrixList.add(res);

                    System.out.println(
                            Thread.currentThread().getName() + " --get-- free matrix num : " + _freeMatrixList.size());
                } else if (_freeMatrixList.size() + _usedMatrixList.size() < _maxPoolSize) {
                    res = new PooledMatrix(nrow, ncol);
                    _usedMatrixList.add(res);

                    System.out.println(Thread.currentThread().getName() + " --get-- increase pool size ");
                }

            }

        } finally {
            lock.unlock();
        }
        if (res == null) {
            res = new Matrix(nrow, ncol);
            System.out.println(
                    Thread.currentThread().getName() + " --get-- empty free PooledMatrixList, allocate Matrix instead");
        }

        return res;

    }

    public boolean release(Matrix mat) {
        lock.lock();
        try {
            if (_ncol == mat.getNCol() && _nrow == mat.getNRow()) {
                _usedMatrixList.remove(mat);
                _freeMatrixList.add(mat);
                System.out.println(
                        Thread.currentThread().getName() + " --release-- free matrix num : " + _freeMatrixList.size());
                return true;
            }

            // if (_usedMatrixList.contains(mat)) {
            // _usedMatrixList.remove(mat);
            // _freeMatrixList.add(mat);
            // System.out.println(
            // Thread.currentThread().getName() + " --release-- free matrix num : " +
            // _freeMatrixList.size());
            // return true;
            // } else if (!_freeMatrixList.contains(mat) && mat.getNCol() == _ncol &&
            // mat.getNRow() == _nrow) {
            // _freeMatrixList.add(mat);
            // System.out.println(
            // Thread.currentThread().getName() + " --release-- free matrix num : " +
            // _freeMatrixList.size());
            // return true;
            // }

        } finally {
            lock.unlock();
        }

        return false;
    }

    // no need to implement this API for now
    // public synchronized boolean setPoolSize(int newSize) {
    //
    // return true;
    // }

    public boolean clear() throws InterruptedException {
        _isOn = false;
        _freeMatrixList.clear();
        _usedMatrixList.clear();
        _minPoolSize = 0;
        _maxPoolSize = 0;
        System.out.println("Matrix pool shuts down completely");
        return true;
    }

}
