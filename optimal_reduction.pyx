# -*- coding: utf-8 -*-
"""
@author: Kacper Pluta kacper.pluta@bdslabs.com.br
Licensed under the terms of the GPL 3 License 
"""

import numpy as np
cimport numpy as np
cimport cython

from libc.math cimport fmin
from libc.math cimport fmax

def computeNLinksCostsReduced(np.ndarray[dtype=np.float32_t, ndim=3, negative_indices=False] costs,
                              np.ndarray[dtype=np.int32_t, ndim=3, negative_indices=False] indices, 
                              np.ndarray[dtype=np.int32_t, ndim=3, negative_indices=False] n_index):
    
    cdef np.ndarray[np.int32_t] index_i = np.empty(0, np.int32)
    cdef np.ndarray[np.int32_t] index_j = np.empty(0, np.int32)
    cdef np.ndarray[np.float32_t] n_costs = np.empty(0, np.float32)

    cdef unsigned int H = n_index.shape[0]
    cdef unsigned int W = n_index.shape[1]
    cdef unsigned int D = n_index.shape[2]
    cdef unsigned int size = H * W * D
    cdef unsigned int i = 0
    cdef unsigned int j = 0
    cdef unsigned int k = 0
    cdef unsigned int c = 0
    cdef unsigned int start0 = 0
    cdef unsigned int start1 = 0
    cdef unsigned int end0 = 0
    cdef unsigned int end1 = 0
    
    while i < H:
        while j < W - 1:
            while k < D:
                if k < D - 1:
                    start0 = n_index[i, j, k]
                    end0 = n_index[i, j, k + 1]
                else: # for terminal group
                    start0 = n_index[i, j, k]
                    end0 = D + 2
                while c < D: # TODO remove this loop and move C to loop one level above
                    if c < D - 1:
                        start1 = n_index[i, j + 1, c]
                        end1 = n_index[i, j + 1, c + 1]
                    else: # for terminal group
                        start1 = n_index[i, j, c]
                        end1 = D + 2
                    if start1 > end0:
                        break
                    intersect = np.arange(fmax(start0,start1),fmin(end0,end1), dtype=np.int32)
                    if intersect.size != 0:
                        np.append(index_i, indices[i, j, k])
                        np.append(index_j, indices[i,j + 1, c])
                        np.append(n_costs,np.sum(costs[i,j, intersect]))
                    c += 1
                k += 1
            j += 1
        i += 1
                        
    return index_i, index_j, n_costs

def hoge(np.ndarray[dtype=np.int32_t,ndim=1],dat):

