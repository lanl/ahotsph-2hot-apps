"""
pgrav_cu.pyx

pgrav_cu defines an m x n interaction

"""

import numpy as np
cimport numpy as np

cdef extern void CUDA_Init()
cdef extern void pinteractCUDA(float *p, float *accp, int n, int stride, 
                               float *f, int source_n, int sz)

def init():
    CUDA_Init()
    return None

#@cython.boundscheck(False)
#@cython.wraparound(False)
def monopole(np.ndarray[float, ndim=1, mode="c"] p not None, \
             np.ndarray[float, ndim=1, mode="c"] accp not None, \
             int stride, \
             np.ndarray[float, ndim=1, mode="c"] f not None):

    m = p.shape[0]/4
    n = f.shape[0]/4
    
    pinteractCUDA(&p[0], &accp[0], m, stride, &f[0], n, 4)
    
    return None
