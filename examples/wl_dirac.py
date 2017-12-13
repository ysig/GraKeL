"""
====================================
Calculating a WL-dirac kernel matrix
====================================

An example plot of :class:`grakel.graph_kernels`
"""
import numpy as np
from grakel.graph_kernels import GraphKernel

X = np.array([[0,1,2,1,0,0.5,0],
              [1,0,0,0,1,2,0.5],
              [2,0,0,3,0,0,2],
              [1,0,3,0,0,0,0],
              [0,1,0,0,0,3,1],
              [0.5,2,0,0,3,0,0],
              [0,0.5,2,0,1,0,0]])

L = {0:'banana', 1:'cherry', 2:'banana', 3:'cherry', 4:'peach',5:'cherry',6:'lime'}

k = 10 
XX = list(zip(k*[X],k*[L]))
base_kernel = dict()
    
gk = GraphKernel(kernel=[{"name":"weisfeiler_lehman","niter":5}, {"name":"dirac"}])
gkf = gk.fit(XX)
print(gkf.transform())
