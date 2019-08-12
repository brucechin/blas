package blas.java;

import java.util.*;
import java.util.concurrent.locks.*;

import blas.java.*;
import java.lang.*;

public class MatrixPool {
    private int _size;
    private int _nrow;
    private int _ncol;
    private boolean _isOn;
    private LinkedList<PooledMatrix> _freeMatrixList;
    private LinkedList<PooledMatrix> _usedMatrixList;
    final Lock lock = new ReentrantLock();

    public static void main(String[] args) throws InterruptedException {
        int nrow = 2000;
        int ncol = 2000;
        int poolsize = 50;
        MatrixPool pool = new MatrixPool(poolsize, nrow, ncol);
        for (int i = 0; i < poolsize; i++) {
            new Thread(new Consumer(pool, nrow, ncol), "Thread " + String.valueOf(i)).start();
        }

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

    public PooledMatrix get(int nrow, int ncol) {
        PooledMatrix res;
        lock.lock();
        try {
            if (_isOn) {
                if (_freeMatrixList.size() > 0 && _ncol == ncol && _nrow == nrow) {
                    res = _freeMatrixList.getFirst();
                    _freeMatrixList.remove();
                    _usedMatrixList.add(res);
                } else {
                    System.out.println("allocating extra memory");
                    res = new PooledMatrix(nrow, ncol);
                }
                System.out.println(
                        Thread.currentThread().getName() + " --get-- free matrix num : " + _freeMatrixList.size());
                return res;
            } else {
                // what shall we return in get() after the pool is down
                res = new PooledMatrix();
                return res;
            }

        } finally {
            lock.unlock();
        }

    }

    public boolean release(PooledMatrix mat) {
        lock.lock();
        try {
            if (_usedMatrixList.contains(mat)) {
                _usedMatrixList.remove(mat);
                _freeMatrixList.add(mat);
                System.out.println(
                        Thread.currentThread().getName() + " --release-- free matrix num : " + _freeMatrixList.size());
                return true;
            } else if (!_freeMatrixList.contains(mat) && mat.getNCol() == _ncol && mat.getNRow() == _nrow) {
                _freeMatrixList.add(mat);
                _size++;
                System.out.println(
                        Thread.currentThread().getName() + " --release-- free matrix num : " + _freeMatrixList.size());
                return true;
            }

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

    public boolean shutDown() throws InterruptedException {
        _isOn = false;
        _freeMatrixList.clear();
        _usedMatrixList.clear();
        _size = 0;
        System.out.println("Matrix pool shuts down completely");
        return true;
    }

}
