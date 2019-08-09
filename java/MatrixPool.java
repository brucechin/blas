package blas.java;

import blas.java.*;
import java.util.*;

import com.sun.istack.internal.Pool;

public class MatrixPool {
    private int _size;
    private int _nrow;
    private int _ncol;
    private LinkedList _freeMatrixList;
    private LinkedList<E> _usedMatrixList;

    public MatrixPool(int poolsize, int nrow, int ncol) {
        _size = poolsize;
        _freeMatrixList = new LinkedList<PooledMatrix>();
        _usedMatrixList = new LinkedList<PooledMatrix>();
        _ncol = ncol;
        _nrow = nrow;
        for (int i = 0; i < _size; i++) {
            PooledMatrix mat = new PooledMatrix(nrow, ncol);
            _freeMatrixList.add(mat);
        }

    }

    public synchronized PooledMatrix get(int nrow, int ncol) {
        PooledMatrix res;

        if (_freeMatrixList.size() > 0 && _ncol == ncol && _nrow == nrow) {
            res = _freeMatrixList.getFirst();
            _freeMatrixList.remove();
            _usedMatrixList.add(res);
        } else {
            res = new PooledMatrix(nrow, ncol);
        }

        return res;
    }

    public synchronized boolean release(PooledMatrix mat) {
        if (_usedMatrixList.contains(mat)) {
            _usedMatrixList.remove(mat);
            _freeMatrixList.add(mat);
            return true;
        }
        // if mat is the same size with other PooledMatrix in the pool, shall we add it
        // into the free list directly?
        return false;

    }

}