# -*- coding: utf-8 -*-
"""
@author: Kacper Pluta kacper.pluta@bdslabs.com.br
Licensed under the terms of the GPL 3 License 
"""

from sys import argv
import pymaxflow
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
from skimage.io import imread
import heapq as h
from skimage import color
import optimal_reduction as opr

import unittest

"""
For stereo matching we are interested only with images of same size
"""
def checkInput(imLeft, imRight):
    assert(imLeft.shape[0] == imRight.shape[0])
    assert(imLeft.shape[1] == imRight.shape[1])
    

def smoothGaussian(im,sigma):
    N=3*sigma
    x=np.arange(-N,N+1).astype(float).reshape(-1,1)
    weights=np.exp(-(x**2)/(2*sigma**2))
    weights=weights/np.sum(weights)
    im2 = ndimage.convolve(im.astype(float), weights)
    im_smooth = ndimage.convolve(im2.astype(float), weights.transpose())  
    return im_smooth    

    
def movingAVG(im, N=3):
    W = N*N
    x= np.ones(N*N).reshape(N,N) / W
    im_avg = ndimage.convolve(im.astype(float), x.astype(float))  
    return im_avg
    
    
def localSum(im, N=3):
    x= np.ones(N*N).reshape(N,N)
    im_avg = ndimage.convolve(im.astype(float), x.astype(float))  
    return im_avg


"""
Local method which return a sum of absolute differences in the given window.
Implemented for comparisions with graph cut like methods. To avoid loops we use
Gaussian smoothing to calculate scores.
"""
def computeSAD(imLeft, imRight, sigma, minDis, maxDis):
    W, H = imLeft.shape[:]

    left = imLeft[:,: -maxDis]
    disMap = np.zeros(left.size).reshape(left.shape[:]) 
    prevScores = np.ones(left.size).reshape(left.shape[:]) * 65532
    
    for i in range( minDis, maxDis + 1 ):
        tmp = np.abs(left - imRight[:,i:left.shape[1] + i])
        gaussian = smoothGaussian(tmp,sigma)
        maskedData = ( ( np.ma.masked_where(prevScores < gaussian, prevScores) * 0 ) + gaussian )
        disMap = ( ( np.ma.masked_array(data=disMap,mask=maskedData.mask) * 0 ) + i ).data
        prevScores = maskedData.data
    return disMap
    
"""
Local method which return a sum of squered differences in the given window.
Implemented for comparisions with graph cut like methods. To avoid loops we use
Gaussian smoothing to calculate scores.
"""
def computeSSD(imLeft, imRight, sigma, minDis, maxDis):
    W, H = imLeft.shape[:]

    
    left = imLeft[:,: -maxDis]
    disMap = np.zeros(left.size).reshape(left.shape[:]) 
    prevScores = np.ones(left.size).reshape(left.shape[:]) * 65532
    
    for i in range( minDis, maxDis + 1 ):
        tmp = (left - imRight[:,i:left.shape[1] + i])**2
        gaussian = smoothGaussian(tmp,sigma)
        maskedData = ( ( np.ma.masked_where(prevScores < gaussian, prevScores) * 0 ) + gaussian )
        disMap = ( ( np.ma.masked_array(data=disMap,mask=maskedData.mask) * 0 ) + i ).data
        prevScores = maskedData.data
    return disMap
    

"""
Local method which return a normalized cross-corelation in the given window.
This is a version proposed by de La Gorce which is a modification of version proposed by Yehu Shen.
Implemented for comparisions with graph cut like methods. To avoid loops we use
convolution with average and local sum like kernel. 
"""
def computeNCC(imLeft, imRight, win, minDis, maxDis):
    
    left = imLeft[:,: -maxDis]
    F4 = localSum(left**2,win)
    F5ii = movingAVG(left**2,win)
    Da = F4 - F5ii    
    
    disMap = np.zeros(left.size).reshape(left.shape[:]) 
    prevScores = np.zeros(left.size).reshape(left.shape[:])
    
    for i in range( minDis, maxDis + 1 ):
        right =imRight[:,i:left.shape[1] + i]
        
        F1 = localSum(left * right,win)
        Nd = F1  - movingAVG(left * right,win)
        
        F2 = localSum(right**2,win)
        F3i = movingAVG(right**2,win)
        
        Db = F2 - F3i

        Dab = np.sqrt(Da * Db)
        tmp = Nd/(np.ma.masked_where(Dab != 0, Dab) + 1).data
        tmp = smoothGaussian(tmp, win)
        maskedData = ( ( np.ma.masked_where(prevScores > tmp, prevScores) * 0 ) + tmp )
        disMap = ( ( np.ma.masked_array(data=disMap,mask=maskedData.mask) * 0 ) + i ).data
        prevScores = maskedData.data
        
    return disMap
      
      
"""
This function returns scorse for graph cut method computed by SAD local method.
Input: imLeft - left image
Input: win - a one side of window
Input: imRight - right image
Input: minDis - start disparity in the range
Input: maxDis - end disparity in the range
"""
def computeCostSAD(imLeft,win, imRight, minDis, maxDis):
    left = imLeft[:,: -maxDis]
    disRange = maxDis - minDis + 1
    disCostsRightToLeft = np.zeros(left.size * disRange).reshape(left.shape[0], left.shape[1], disRange)
    
    for i in range (minDis, maxDis + 1):
        right = imRight[:,i:left.shape[1] + i]
        tmp = np.abs(left - right)
        disCostsRightToLeft[:,:,i - minDis] = localSum(tmp,win)
        
    return disCostsRightToLeft
    
"""
This function returns scorse for graph cut method computed by SSD local method.
Input: imLeft - left image
Input: win - a one side of window
Input: imRight - right image
Input: minDis - start disparity in the range
Input: maxDis - end disparity in the range
"""
def computeCostSSD(imLeft,win, imRight, minDis, maxDis):
    left = imLeft[:,: -maxDis]
    disRange = maxDis - minDis + 1
    disCostsRightToLeft = np.zeros(left.size * disRange).reshape(left.shape[0], left.shape[1], disRange)
    
    for i in range (minDis, maxDis + 1):
        right = imRight[:,i:left.shape[1] + i]
        tmp = (left - right)**2
        disCostsRightToLeft[:,:,i - minDis] = localSum(tmp,win)
        
    return disCostsRightToLeft
    

"""
This function returns scorse for graph cut method computed by NCC local method.
This is a version proposed by de La Gorce which is a modification of version proposed by Yehu Shen.
Input: imLeft - left image
Input: win - a one side of window
Input: imRight - right image
Input: minDis - start disparity in the range
Input: maxDis - end disparity in the range
"""  
def computeCostNCC(imLeft,win, imRight, minDis, maxDis):
    left = imLeft[:,: -maxDis]
    disRange = maxDis - minDis + 1
    F4 = localSum(left**2,win)
    F5ii = movingAVG(left**2,win)
    Da = F4 - F5ii
    
    disCostsRightToLeft = np.empty(left.size * disRange).reshape(left.shape[0], left.shape[1], disRange)
    
    for i in range (minDis, maxDis + 1):
        right = imRight[:,i:left.shape[1] + i]
        
        F1 = localSum(left * right,win)
        Nd = F1  - movingAVG(left * right,win)
        
        F2 = localSum(right**2,win)
        F3i = movingAVG(right**2,win)
        
        Db = F2 - F3i

        Dab = np.sqrt(Da * Db)
        disCostsRightToLeft[:,:,i - minDis] = Nd/(np.ma.masked_where(Dab != 0, Dab) + 1).data
        
    return disCostsRightToLeft
    
"""
This function returns N highest scores and their indices
Input: N - number of scores to keep
Input: costs - costs to keep
Input: costs1 - weights of costs - this function looks only on them
"""
def leaveOnly(costs, costs1, N):
    H = costs.shape[0]
    W = costs.shape[1]    
    n_costs = np.ones(W * H * N).reshape(H, W, N)
    n_index = np.ones(W * H * N).reshape(H, W, N)
    
    for i in range(0,H):
        for j in range(0,W):
            indices = [t[0] for t in h.nlargest(N, enumerate(costs1[i,j]), lambda t: t[1])]
            indices.sort()
            n_index[i,j,:] = indices;
            n_costs[i,j,range(0,N)] = costs[i,j,indices]
            
    
    return n_costs, n_index

"""
This function returns N highest scores and their indices
Input: N - number of scores to keep
Input: costs - costs to keep
Input: costs1 - weights of costs - this function looks only on them
Input: inhibitor - scores which we want to keep in output N scores
"""    
def leaveOnlyAndPreserve(costs, costs1, N, inhibitor):
    H = costs.shape[0]
    W = costs.shape[1]    
    n_costs = np.ones(W * H * N).reshape(H, W, N)
    n_index = np.ones(W * H * N).reshape(H, W, N)
    
    for i in range(0,H):
        for j in range(0,W):
            indices = [t[0] for t in h.nlargest(N, enumerate(costs1[i,j]), lambda t: t[1])]
            if np.size(np.where(indices != inhibitor[i,j] )) == 0:
                indices[N-1] = inhibitor[i,j] # Change last one to the value from inhibitor
                print "Preserved disparity:" + str(inhibitor[i,j])
            indices.sort()
            n_index[i,j,:] = indices;
            n_costs[i,j,range(0,N)] = costs[i,j,indices]
            
    
    return n_costs, n_index
    

"""
Graph reduction scheme proposed by de La Gorce
!!!Cython version available!!!
Input: costs - whole set of scores for whole range of disparities
Input: indices of vertices like for full graph scheme
Input: n_index - indices of N best disparities
""" 
def computeNLinksCostsReduced(costs, indices, n_index):
    index_i = []
    index_j = []
    n_costs = []
    
    H, W, D = n_index.shape[:]
    for i in range(0, H):
        for j in range(0, W - 1):
            for k in range(0, D):
                if k < D - 1:
                    start0 = n_index[i, j, k]
                    end0 = n_index[i, j, k + 1]
                else: # for terminal group
                    start0 = n_index[i, j, k]
                    end0 = D + 2
                for c in range(0, D):
                    if c < D - 1:
                        start1 = n_index[i, j + 1, c]
                        end1 = n_index[i, j + 1, c + 1]
                    else: # for terminal group
                        start1 = n_index[i, j, c]
                        end1 = D + 2
                    if start1 > end0:
                        break
                    intersect = np.arange(max(start0,start1),min(end0,end1),dtype=np.int32)
                    if intersect.size != 0:
                        index_i.append(indices[i, j, k])
                        index_j.append(indices[i,j + 1, c])
                        n_costs.append(np.sum(costs[i,j, intersect]))
                        
    return np.array(index_i).astype(np.int32), np.array(index_j).astype(np.int32), np.array(n_costs).astype(np.float32)
    

"""
Implementation of stereo matching with a full graph scheme

@inproceedings{zureiki2007SteMatusiRedCut,
	author = {Zureiki, Ayman and Devy, Michel and Chatila, Raja},
	booktitle = {{ICIP (1)}},
	pages = {237–240},
	publisher = {IEEE},
	title = {{Stereo Matching using Reduced-Graph Cuts.}},
	year = 2007
}

Small modification also decribed in this paper but proposed by (Ishikawa, 2000)

Input: minDis - start disparity in the range
Input: maxDis - end disparity in the range
Input: costs0 - scores computed by local method
Input: imLeft - left image
Input: alpha - varible which allows to change importance of n-links
"""
def computeDispMapByGraph(disMin, disMax, costs0, imLeft, alpha):

    W = costs0.shape[0]
    H = costs0.shape[1]
    disRange = disMax - disMin + 1    
    nbpixels = W * H
    V = nbpixels * disRange
    # Thus we are working with oriented graph we need more edges, thus note V * 6 but V * 12
    E = 12 * V - 4 * disRange * (W + H)
    
    g = pymaxflow.PyGraph(V,E)
    g.add_node(V)
    indices = np.arange(V).reshape(W, H, disRange)   
    
   #nodes, source, sink
    """ is just calling add_tweights(index[i],source[i],sink[i]) -- create terminal
    edges which connect s and t with non terminal nodes. To avoid cut near to sink we just add a big weights
    """
    g.add_tweights_vectorized(indices[:,:,0].flatten().astype(np.int32), costs0[:,:,0].flatten(), 0 * costs0[:,:,-1].flatten())
    g.add_tweights_vectorized(indices[:,:,-1].flatten().astype(np.int32), 0 * costs0[:,:,0].flatten(), costs0[:,:,-1].flatten() + 1000)
        
    # Add t-links - To eliminate K term we use oriented graph with a big cost for a one egde
    indices_node_i = indices[:,:,1:].flatten().astype(np.int32)
    indices_node_j = indices[:,:,:-1].flatten().astype(np.int32)
    cost_diff0 = costs0[:,:,1:].flatten()
    g.add_edge_vectorized(indices_node_i, indices_node_j, cost_diff0.astype(np.float32) + 1000, cost_diff0.astype(np.float32))    
    
    # add n-links - horizontal
    indices_node_i = indices[:-1,:].ravel().astype(np.int32)
    indices_node_j = indices[1:,:].ravel().astype(np.int32)
    
    # U term represents a gradient of intencity between neighbor pixels
    U = 1-np.abs(imLeft[:-1,:] - imLeft[1:,:]).T #U term has to be low for big difference of intensity
    # V term represents a gradient of scores neighbor pixels
    V = np.abs(costs0[:-1,:] - costs0[1:,:])
    # Multiply costs for each disparity by U term
    cost_diff = np.array([U * v for v in V.T]).T.flatten().astype(np.float32) * alpha
    
    g.add_edge_vectorized(indices_node_i, indices_node_j, cost_diff ,cost_diff)
    
    # add n-links - vertical
    indices_node_i = indices[:,:-1].ravel().astype(np.int32)
    indices_node_j = indices[:,1:].ravel().astype(np.int32)
    
    # U term represents a gradient of intencity between neighbor pixels
    U = 1-np.abs(imLeft[:,:-1] - imLeft[:,1:]).T #U term has to be low for big difference of intensity
    # V term represents a gradient of scores neighbor pixels
    V = np.abs(costs0[:,:-1] - costs0[:,1:])
    # Multiply costs for each disparity by U term
    cost_diff = np.array([U * v for v in V.T]).T.flatten().astype(np.float32) * alpha
    
    g.add_edge_vectorized(indices_node_i, indices_node_j, cost_diff ,cost_diff)    

    g.maxflow()

    out = g.what_segment_vectorized()
    out = out.reshape(W, H, disRange)
    im = (disMax - disMin) - np.sum(out,axis=2) # calculate the final result by suming scores and moving them to proper range
        
    return im, out
    
"""
Implementation of stereo matching with a reduced graph scheme

@inproceedings{zureiki2007SteMatusiRedCut,
	author = {Zureiki, Ayman and Devy, Michel and Chatila, Raja},
	booktitle = {{ICIP (1)}},
	pages = {237–240},
	publisher = {IEEE},
	title = {{Stereo Matching using Reduced-Graph Cuts.}},
	year = 2007
}

Small modification also decribed in this paper but proposed by (Ishikawa, 2000)
Input: costs0 - N best scores computed by local method
Input: imLeft - left image
Input: n_index - indices of N best scores in whole matrix of indices
Input: alpha - varible which allows to change importance of n-links
"""
def computeDispMapByReducedGraph(costs0, imLeft, n_index, alpha):

    W = costs0.shape[0]
    H = costs0.shape[1]
    V = n_index.size
    rangeDisp = n_index.shape[2]
    # Thus we are working with oriented graph we need more edges, thus note V * 6 but V * 12
    E = 12 * V - 4 * rangeDisp * (W + H)
    
    g = pymaxflow.PyGraph(V,E)
    g.add_node(V)
    indices = np.arange(V).reshape(W, H, rangeDisp)
    
    #nodes, source, sink
    """ is just calling add_tweights(index[i],source[i],sink[i]) -- create terminal
    edges which connect s and t with non terminal nodes. To avoid cut near to sink we just add a big weights
    """
    g.add_tweights_vectorized(indices[:,:,0].flatten().astype(np.int32), costs0[:,:,0].flatten(),0 * costs0[:,:,-1].flatten())
    g.add_tweights_vectorized(indices[:,:,-1].flatten().astype(np.int32), 0 * costs0[:,:,0].flatten(), costs0[:,:,-1].flatten() + 1000)
        
    # Add t-links - To eliminate C term we use oriented graph with a big cost for a one egde
    indices_node_i = indices[:,:,1:].flatten().astype(np.int32)
    indices_node_j = indices[:,:,:-1].flatten().astype(np.int32)
    cost_diff0 = costs0[:,:,1:].flatten()
    g.add_edge_vectorized(indices_node_i, indices_node_j, cost_diff0.astype(np.float32) + 1000, cost_diff0.astype(np.float32))   

    # add n-links - horizontal
    indices_node_i = indices[:-1,:].ravel().astype(np.int32)
    indices_node_j = indices[1:,:].ravel().astype(np.int32)

    # U term represents a gradient of intencity between neighbor pixels
    U = 1-np.abs(imLeft[:-1,:] - imLeft[1:,:]).T #U term has to be low for big difference of intensity
    # V term represents a difference between disparities of neighbor pixels
    V = np.abs(n_index[:-1,:] - n_index[1:,:]) + 1
    # Multiply costs for each disparity by U term
    cost_diff = np.array([U * v for v in V.T]).T.flatten().astype(np.float32) * alpha
    
    g.add_edge_vectorized(indices_node_i, indices_node_j, cost_diff ,cost_diff)
    
    # add n-links - vertical
    indices_node_i = indices[:,:-1].ravel().astype(np.int32)
    indices_node_j = indices[:,1:].ravel().astype(np.int32)
    
    # U term represents a gradient of intencity between neighbor pixels
    U = 1-np.abs(imLeft[:,:-1] - imLeft[:,1:]).T #U term has to be low for big difference of intensity
    # V term represents a difference between disparities of neighbor pixels
    V = np.abs(n_index[:,:-1] - n_index[:,1:]) + 1
    # Multiply costs for each disparity by U term
    cost_diff = np.array([U * v for v in V.T]).T.flatten().astype(np.float32) * alpha            

    g.maxflow()
    out = g.what_segment_vectorized()
    out = out.reshape(n_index.shape[:])
    im = np.zeros(W * H).reshape(W, H)
    
    for i in range(0,W):
        for j in range(0,H):
               if np.size(np.where(out[i,j,:] !=0)) != 0:
                   im[i,j] = n_index[i,j, np.where(out[i,j,:] !=0)[0][0]] # look for good disparity
        
    return im
    
    
"""
Implementation of stereo matching with a reduced graph scheme

Reduced graph scheme proposed by de La Gorce

Input: costs0 - N best scores computed by local method
Input: Leftimage - left image
Input: n_index - indices of N best scores in whole matrix of indices
Input: alpha - varible which allows to change importance of n-links
Input: range of whole disparities
"""    
def computeDispMapByReducedGraph2(costs0, Leftimage, n_index, alpha, dispRang):

    W = costs0.shape[0]
    H = costs0.shape[1]
    V = n_index.size
    rangeDisp = n_index.shape[2]
    
    Leftimage = Leftimage.astype(np.float32)    
    n_index = n_index.astype(np.int32)      
    
    indices = np.arange(V).reshape(W, H, rangeDisp).astype(np.int32)
    mask = np.ones(W * H * dispRang).reshape(W, H, dispRang)
    # Calculate absolute differences in left image and propagate them by multiplication with ones
    U = 1-np.abs(Leftimage[:,:-1] - Leftimage[:,1:])
    cost_diff = np.array([U.T * v[:-1,:] for v in mask.T]).T.astype(np.float32)
    # Calculate horizontal n-links
    indices_node_i_h, indices_node_j_h, cost_diff_h = opr.computeNLinksCostsReduced(cost_diff,indices,n_index)
    
    # Calculate absolute differences in left image and propagate them by multiplication with ones
    U = 1-np.abs(Leftimage[:-1,:] - Leftimage[1:,:])
    cost_diff = np.array([U.T * v[:,:-1] for v in mask.T]).T.astype(np.float32)
    # Calculate vertical n-links
    indices_node_i_v, indices_node_j_v, cost_diff_v = opr.computeNLinksCostsReduced(cost_diff.transpose(1,0,2),indices.transpose(1,0,2),n_index.transpose(1,0,2))
    
    # Calculate edges directly from output data
    E = indices_node_i_h.size * 2 + indices_node_i_v.size * 2 + V * 2
    
    g = pymaxflow.PyGraph(V,E)
    g.add_node(V)
    
    #nodes, source, sink
    """ is just calling add_tweights(index[i],source[i],sink[i]) -- create terminal
    edges which connect s and t with non terminal nodes. To avoid cut near to sink we just add a big weights
    """
    g.add_tweights_vectorized(indices[:,:,0].flatten().astype(np.int32), costs0[:,:,0].flatten(),0 * costs0[:,:,-1].flatten())
    g.add_tweights_vectorized(indices[:,:,-1].flatten().astype(np.int32), 0 * costs0[:,:,0].flatten(), costs0[:,:,-1].flatten() + 1000)
        
    # Add t-links - To eliminate C term we use oriented graph with a big cost for a one egde
    indices_node_i = indices[:,:,1:].flatten().astype(np.int32)
    indices_node_j = indices[:,:,:-1].flatten().astype(np.int32)
    cost_diff0 = costs0[:,:,1:].flatten()
    g.add_edge_vectorized(indices_node_i, indices_node_j, cost_diff0.astype(np.float32) + 1000, cost_diff0.astype(np.float32))   

    # add n-links - horizontal
    cost_diff_h *= alpha
    g.add_edge_vectorized(indices_node_i_h, indices_node_j_h, cost_diff_h ,cost_diff_h)
    # add n-links - vertical
    cost_diff_v *= alpha
    g.add_edge_vectorized(indices_node_i_v, indices_node_j_v, cost_diff_v, cost_diff_v)
    
    g.maxflow()
    out = g.what_segment_vectorized()
    out = out.reshape(n_index.shape[:])
    im = np.zeros(W * H).reshape(W, H)
    
    for i in range(0,W):
        for j in range(0,H):
               if np.size(np.where(out[i,j,:] !=0)) != 0:
                   im[i,j] = n_index[i,j, np.where(out[i,j,:] !=0)[0][0]] # look for good disparity
        
    return im

"""
Unit tests for de La Gorce scheme
"""
class PowerReducedGraphTestCase2D(unittest.TestCase):
    index_i = None
    index_j = None
    n_costs = None    
    def setUp(self):
        n_costs2D = np.ones(24).reshape(1,4,6)
        n_index2D = np.array([[0,1,4,5],
                        [0,2,4,5],
                        [0,3,4,5],
                        [0,2,3,5]]).reshape(1,4,4)
        indices2D = np.arange(n_costs2D.size).reshape(n_costs2D.shape[:])
        self.index_i, self.index_j, self.n_costs = computeNLinksCostsReduced(n_costs2D, indices2D, n_index2D)
    def CheckIndices_I(self):
        assert( (self.index_i == np.array([0, 1, 1, 2, 3, 6, 7, 7, 8, 9, 12, 12, 13, 14, 15])).all())
    def CheckIndices_J(self):
        assert( (self.index_j == np.array([6, 6, 7, 8, 9, 12, 12, 13, 14, 15, 18, 19, 20, 20, 21])).all())
    def CheckCosts(self):
        assert( (self.n_costs == np.array([1.0, 1.0, 2.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 1.0])).all())
        
"""
Unit tests for de La Gorce scheme
"""        
class PowerReducedGraphTestCase3D(unittest.TestCase):
    index_i = None
    index_j = None
    n_costs = None    
    def setUp(self):
        n_costs3D = np.ones(24 * 2).reshape(2,4,6)
        n_index3D = np.empty(4 * 4 * 2).reshape(2,4,4)
    
        for i in range(0, n_index3D.shape[0]):
            n_index3D[i,:,:] = [[0,1,4,5], [0,2,4,5], [0,3,4,5], [0,2,3,5]]
                        
        indices3D = np.arange(n_costs3D.size).reshape(n_costs3D.shape[:])
        self.index_i, self.index_j, self.n_costs = computeNLinksCostsReduced(n_costs3D, indices3D, n_index3D.astype(int))
    def CheckIndices_I(self):
        assert( (self.index_i == np.array([0, 1, 1, 2, 3, 6, 7, 7, 8, 9, 12, 12, 13, 14, 15, 24, 25, 25, 26, 27, 30, 31, 31, 32, 33, 36, 36, 37, 38, 39])).all())
    def CheckIndices_J(self):
        assert( (self.index_j == np.array([6, 6, 7, 8, 9, 12, 12, 13, 14, 15, 18, 19, 20, 20, 21, 30, 30, 31, 32, 33, 36, 36, 37, 38, 39, 42, 43, 44, 44, 45])).all())
    def CheckCosts(self):
        assert( (self.n_costs == np.array([1.0, 1.0, 2.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 1.0])).all())
       
"""
Unit tests for de La Gorce scheme
"""       
class PowerReducedGraphTestCase3DVertical(unittest.TestCase):
    index_i = None
    index_j = None
    n_costs = None    
    def setUp(self):
        n_costs3D = np.ones(24 * 2).reshape(2,4,6)
        n_index3D = np.empty(4 * 4 * 2).reshape(2,4,4)
    
        for i in range(0, n_index3D.shape[0]):
            n_index3D[i,:,:] = [[0,1,4,5], [0,2,4,5], [0,3,4,5], [0,2,3,5]]
                        
        indices3D = np.arange(n_costs3D.size).reshape(n_costs3D.shape[:])
        self.index_i, self.index_j, self.n_costs = computeNLinksCostsReduced(n_costs3D.transpose(1,0,2), indices3D.transpose(1,0,2), n_index3D.transpose(1,0,2).astype(int))
    def CheckIndices_I(self):
        assert( (self.index_i == np.array([0, 1, 2, 3, 6, 7, 8, 9, 12, 13, 14, 15, 18, 19, 20, 21])).all())
    def CheckIndices_J(self):
        assert( (self.index_j == np.array([24, 25, 26, 27, 30, 31, 32, 33, 36, 37, 38, 39, 42, 43, 44, 45])).all())
    def CheckCosts(self):
        assert( (self.n_costs == np.array([1.0, 3.0, 1.0, 1.0, 2.0, 2.0, 1.0, 1.0, 3.0, 1.0, 1.0, 1.0, 2.0, 1.0, 2.0, 1.0])).all())
  

"""
Buld unit tests
"""
def suite():
    suite = unittest.TestSuite()
    suite.addTest(PowerReducedGraphTestCase2D("CheckIndices_I"))
    suite.addTest(PowerReducedGraphTestCase2D("CheckIndices_J"))
    suite.addTest(PowerReducedGraphTestCase2D("CheckCosts"))
    suite.addTest(PowerReducedGraphTestCase3D("CheckIndices_I"))
    suite.addTest(PowerReducedGraphTestCase3D("CheckIndices_J"))
    suite.addTest(PowerReducedGraphTestCase3D("CheckCosts"))
    suite.addTest(PowerReducedGraphTestCase3DVertical("CheckIndices_I"))
    suite.addTest(PowerReducedGraphTestCase3DVertical("CheckIndices_J"))
    suite.addTest(PowerReducedGraphTestCase3DVertical("CheckCosts"))
    return suite
    
"""
Comapre outputs TODO remove in the future
"""
def isInFullGraph(fullGraphOut,imFromReduce):
    H = fullGraphOut.shape[0]
    W = fullGraphOut.shape[1]
    truth = np.empty(H * W).reshape(H, W)
    for i in range(0,fullGraphOut.shape[0]):
        for j in range(0,fullGraphOut.shape[1]):
            if fullGraphOut[i,j, imFromReduce[i,j]] == 1:
                truth[i,j] = True
            else:
                truth[i,j] = False
            
    print truth.all()
    return truth


def main():
#    script, disMin, disMax, alpha, K, leftImagePath, rightImagePath = argv
    
#    disMin = np.int(disMin)
#    disMax = np.int(disMax)
    
#    alpha = np.float(alpha)
    
    disMin = 0
    disMax = 30
    alpha = 0.1
    
    leftImagePath = "./art/view1.png"
    rightImagePath = "./art/view0.png"
    
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
    imLeft=np.array(imread(leftImagePath))
    imRight=np.array(imread(rightImagePath))
    
    checkInput(imLeft, imRight)
    
#    We will work only with a grayscale images
    imLeftGray = color.rgb2gray(imLeft)
    imRightGray = color.rgb2gray(imRight)
    imLeft = imLeft.astype(np.float)
    imRight = imRight.astype(np.float)
    
    plt.figure()
    plt.imshow(imLeftGray,cmap=plt.cm.Greys_r) 
    plt.title('Left image')
    plt.show()
    plt.figure()
    plt.imshow(imRightGray,cmap=plt.cm.Greys_r) 
    plt.title('Right image')
    plt.show()
    
    # Compute disparities by local methods
    disMap = computeSAD(imLeftGray, imRightGray, 3, disMin, disMax)
    
    plt.figure()
    plt.imshow(disMap,cmap=plt.cm.Greys_r) 
    plt.title('Disparity map - SAD')
    plt.show()
    
    disMap = computeNCC(imLeftGray, imRightGray,3, disMin, disMax)
    
    plt.figure()
    plt.imshow(disMap,cmap=plt.cm.Greys_r) 
    plt.title('Disparity map - NCC')
    plt.show()
    
    disMap = computeSSD(imLeftGray, imRightGray,3, disMin, disMax)
    
    plt.figure()
    plt.imshow(disMap,cmap=plt.cm.Greys_r) 
    plt.title('Disparity map - SSD')
    plt.show()
    
    # Compute scores for the full graph    
    
#    costs0 = 1-computeCostNCC(imLeftGray, 3, imRightGray, disMin, disMax)
    costs0 = computeCostSAD(imLeftGray, 7, imRightGray, disMin, disMax)
#    costs0 = computeCostSSD(imLeftGray, 3, imRightGray, disMin, disMax)
    im, out = computeDispMapByGraph(disMin,disMax, costs0.astype(np.float32),imLeftGray[:,: -disMax].astype(np.float32), alpha)
    
    plt.figure()
    plt.imshow(im,cmap=plt.cm.Greys_r)
    plt.title('Disparity map - Graph cut - NCC FULL')
    plt.show()
    
    # Compute scores for reduced graphs
    costs0 = computeCostSAD(imLeftGray, 7, imRightGray, disMin, disMax)
    ncc = computeCostNCC(imLeftGray, 7, imRightGray, disMin, disMax)    
    
    n_costs, n_index = leaveOnly(costs0, ncc, 4)
#    n_costs, n_index = leaveOnlyAndPreserve(costs0, ncc, 4, im)
    n_costs = n_costs.astype(np.float32)
    # Reduction method -- from paper, look in documentation
    im = computeDispMapByReducedGraph(n_costs, imLeftGray[:,: -disMax].astype(np.float32), n_index, alpha)
        
    plt.figure()
    plt.imshow(im,cmap=plt.cm.Greys_r)
    plt.title('Disparity map - Graph cut - SAD/NCC REDUCED')
    plt.show()
    # Reduction method by de La Gorce scheme
    im = computeDispMapByReducedGraph2(n_costs, imLeftGray[:,: -disMax].astype(np.float32), n_index, alpha, disMax - disMin)

#    isInFullGraph(out,im)
        
    plt.figure()
    plt.imshow(im,cmap=plt.cm.Greys_r)
    plt.title('Disparity map - Graph cut - SAD/NCC REDUCED2')
    plt.show()
    
if __name__ == "__main__":
    main()