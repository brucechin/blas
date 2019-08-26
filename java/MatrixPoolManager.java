package blas.java;

import java.util.*;
import java.util.concurrent.locks.*;
import java.util.concurrent.*;
import blas.java.*;
import java.lang.*;

public class MatrixPoolManager {

    // fixed min/max pool size
    private int _minPoolSize = 5;
    private int _maxPoolSize = 30;
    private ConcurrentHashMap<String, MatrixPool> _poolMap;// key is a string "nrow+ncol"

    private static MatrixPoolManager manager;

    public static void main(String[] args) {
        int threadNum = 100;
        for (int i = 0; i < threadNum; i++) {
            new Thread(new Consumer((i % 5 + 1) * 100, (i % 5 + 1) * 100), "Thread " + String.valueOf(i)).start();
        }
    }

    public Matrix get(int nrow, int ncol) {
        Matrix res = null;
        String key = String.valueOf(nrow) + "+" + String.valueOf(ncol);

        if (!_poolMap.containsKey(key)) {
            // allocate new matrix pool for this (nrow, ncol)
            _poolMap.put(key, new MatrixPool(_minPoolSize, _maxPoolSize, nrow, ncol));
        }

        res = _poolMap.get(key).get(nrow, ncol);

        return res;
    }

    public boolean release(Matrix mat) {
        String key = String.valueOf(mat.getNRow()) + "+" + String.valueOf(mat.getNCol());
        _poolMap.get(key).release(mat);
        return false;
    }

    private MatrixPoolManager() {
        _poolMap = new ConcurrentHashMap<String, MatrixPool>();
    }

    public static synchronized MatrixPoolManager getMatrixPoolManager() {
        if (manager == null) {
            manager = new MatrixPoolManager();
        }
        return manager;
    }

}