#!/usr/bin/env python

import numpy as np
import pgrav

m = 1
n = 16
np.random.seed(123)

source = np.empty((n,4), dtype=np.float32)
source[:,0] = 1.0 / n
if 1:
    source[:,1] = np.random.uniform( 0.0,1.0,n)
    source[:,2] = np.random.uniform(-0.5,0.5,n)
    source[:,3] = np.random.uniform(-0.5,0.5,n)
else:
    source[:,1] = np.arange(n)
    source[:,2] = n+np.arange(n)
    source[:,3] = 2*n+np.arange(n)

sink = np.array([1.0, 0.0, 0.0, 0.0], dtype=np.float32)
accp = np.zeros((4), dtype=np.float32)

sink_stride = 4
accp = np.zeros((4), dtype=np.float32)
pgrav.monopole(sink, accp, sink_stride, source.flatten())
print accp

nsse = 8
source_swiz = source.reshape(-1,nsse,sink_stride).transpose(0,2,1)

accp = np.zeros((4), dtype=np.float32)
pgrav.vec_monopole(sink, accp, sink_stride, source_swiz.flatten())
print accp
