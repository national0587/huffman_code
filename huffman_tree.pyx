# -*- coding: utf-8 -*-
# distutils: language=c++
# distutils: sources=huffman_tree_.cpp
"""
@author: Kacper Pluta kacper.pluta@bdslabs.com.br
Licensed under the terms of the GPL 3 License
https://bl.ocks.org/benjaminirving/436262a58f9da5a68532
"""



#import numpy as np
#cimport numpy as np
cimport cython
from libcpp.vector cimport vector
#from libc.math cimport fmin
#from libc.math cimport fmax


cdef extern from "huffman_tree_.h" namespace "huffman_test":
    cdef cppclass INode:
        int f
        INode()

    cdef cppclass HuffmanCpp:
        HuffmanCpp() except +
        void makeCodeBook(vector[int])
        #getCodeMap()
        void Vec2CodeVec(vector[int], bitstream)
    cdef cppclass bitstream:
        bitstream()
        void reset()
        void pushbit(char c)
        void flushbit()
        char readbit()
        vector[char] data
        int data_size
        int last_bit_offset



cdef class PyHuffman:
    cdef HuffmanCpp *thisptr
    cdef bitstream bs
    def __cinit__(self):
        self.thisptr = new HuffmanCpp()

    def __dealloc__(self):
        del self.thisptr

    # def out(self, input_list):
    #     cdef vector[int] vect = input_list
    #     self.thisptr.makeCodeBook(vect)
    def makeCodeBook(self, input_vec):
        cdef vector[int] vec = input_vec
        self.thisptr.makeCodeBook(vec)
    def V2CV(self, llist):
        cdef vector[int] hoge = llist
        self.thisptr.Vec2CodeVec(hoge,self.bs)
        return self.bs.data, self.bs.data_size


