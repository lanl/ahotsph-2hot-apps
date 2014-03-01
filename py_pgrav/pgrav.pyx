"""
pgrav.pyx

pgrav defines an m x n interaction

"""

import numpy as np
cimport numpy as np

cdef extern void _pgravm(float *p, float *accp, int n, int stride, 
                         float *f, int source_n)

cdef extern int pgrav_vec_init()

cdef extern void pMinteract(float *p, float *accp, int n, int stride, 
                            float *f, int source_n)

#@cython.boundscheck(False)
#@cython.wraparound(False)
def monopole(np.ndarray[float, ndim=1, mode="c"] p not None, \
             np.ndarray[float, ndim=1, mode="c"] accp not None, \
             int stride, \
             np.ndarray[float, ndim=1, mode="c"] f not None):

    m = p.shape[0]/4
    n = f.shape[0]/4
    
    _pgravm(&p[0], &accp[0], m, stride, &f[0], n)
    
    return None

#@cython.boundscheck(False)
#@cython.wraparound(False)
def vec_monopole(np.ndarray[float, ndim=1, mode="c"] p not None, \
                 np.ndarray[float, ndim=1, mode="c"] accp not None, \
                 int stride, \
                 np.ndarray[float, ndim=1, mode="c"] f not None):
    
    nsse = pgrav_vec_init()

    m = p.shape[0]/4
    n = f.shape[0]/(4*nsse)
    
    pMinteract(&p[0], &accp[0], m, stride, &f[0], n)
    
    return None
