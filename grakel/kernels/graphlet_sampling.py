""" This file contains the shortest graphlet as defined in :cite:`Shervashidze2009EfficientGK`
    
"""

import math
import itertools
import random
import warnings

import numpy as np

import pynauty

from ..graph import graph
from ..tools import matrix_to_dict

random.seed(15487469)
np.random.seed(15487103)

def graphlet_sampling(X, Y, k=5, delta=0.05, epsilon=0.05, a=-1):
    """ Applies the graphlet sampling kernel as proposed in :cite:`Shervashidze2009EfficientGK`.
        
    
    Parameters
    ----------
    X,Y : valid graph format
        The pair of graphs on which the kernel is applied.
        
    k : int, default=5
        The dimension of the given graphlets.
        
    delta : float, default=0.05
        Confidence level (typically 0.05 or 0.1).
    
    epsilon : float, default=0.05
        Precision level (typically 0.05 or 0.1).
        
    a : int
        Number of isomorphism classes of graphlets.
        If -1 the number is the maximum possible, from a database 1 until 9
        or else predicted through interpolation.
    
    Returns
    -------
    kernel : number
        The kernel value.
    """
    Gx = graph(X)
    Gy = graph(Y)
    
    return graphlet_sampling_inner(Gx, Gy, k=k, delta=delta, epsilon=epsilon, a=a)
    
def graphlet_sampling_inner(Gx, Gy, k=5, delta=0.05, epsilon=0.05, a=-1):
    """ Applies the graphlet sampling kernel, as proposed in :cite:`Shervashidze2009EfficientGK`.

    Parameters
    ----------
    G{x,y} : graph
        The pair of graphs on which the kernel is applied.
        
    k : int, default=5
        The dimension of the given graphlets.
        
    delta : float, default=0.05
        Confidence level (typically 0.05 or 0.1).
    
    epsilon : float, default=0.05
        Precision level (typically 0.05 or 0.1).
        
    a : int
        Number of isomorphism classes of graphlets.
        If -1 the number is the maximum possible, from a database 1 until 9
        or else predicted through interpolation.
    
    Returns
    -------
    kernel : number
        The kernel value.

    """
    # Steps:
    # Feature Space
    nsamples, graphlets, P, graph_bins, nbins = sample_graphlets(k, delta, epsilon, a)

    # Calculate Features
    return graphlet_sampling_core(Gx, Gy, nsamples, graphlets, P, graph_bins, nbins, k)

def graphlet_sampling_core(Gx, Gy, nsamples, graphlets, P, graph_bins, nbins, k=5):
    """ Applies the graphlet sampling kernel given a graphlet feature space, as proposed in :cite:`Shervashidze2009EfficientGK`.
        
    Parameters
    ----------
    G{x,y} : graph
        The pair of graphs on which the kernel is applied.
    nsamples : int
        The number of samples in the graphlet feature space.
    graphlets : dict
        A dictionary of graphlets with keys from 0 ... nsamples-1.
    P : np.array
        A matrix that has zeros where the pair of indexes corresponds to isomorphic graphlets and 1 otherwise.
    graph_bins : dict
        A dictionary of bined graphlets, based on isomorphism classes.
    nbins : int
        The number of isomorphism classes.
    k : int
        Graphlet dimension.
    
    Returns
    -------
    kernel : number
        The kernel value.
    """
    if k<1:
        raise ValueError('k must be bigger than 1')
    
    
    Gx.desired_format("adjacency")
    Gy.desired_format("adjacency")
    X, Y = Gx.adjacency_matrix, Gy.adjacency_matrix    

    # Calculate frequencies for each graph
    # Check that if each matrix is a principal
    # minor of the adjoint and how many times
    i = 0

    # Frequency vector for x
    fx = np.zeros(nsamples)
    
    # Frequency vector for y
    fy = np.zeros(nsamples)

    # To transform the adjacency to edge dictionary
    # needed for nauty graph initialise a small lambda
    to_edge_dict_real = lambda x : matrix_to_dict(x, '>', .0, k, False)

    # For all kminors
    for c in itertools.combinations(list(range(k)),k):
        for s in range(0,nbins):
            idxs = list(c)
            # Extract k minors
            X_m = X[idxs,:][:,idxs]
            Y_m = Y[idxs,:][:,idxs]
            graphlet_idx = graph_bins[s][0]
            # Test isomorphism with each graphlet
            if(pynauty.isomorphic(pynauty.Graph(k, True, to_edge_dict_real(X_m)), graphlets[graphlet_idx])):
                fx[graphlet_idx]+=1
            if(pynauty.isomorphic(pynauty.Graph(k, True, to_edge_dict_real(Y_m)), graphlets[graphlet_idx])):
                fy[graphlet_idx]+=1
                
    for s in range(0,nbins):
        for i in range(1,len(graph_bins[s])):
            idx = graph_bins[s][i]
            fx[idx]=fx[graph_bins[s][0]]
            fy[idx]=fy[graph_bins[s][0]]

    # normalize fx
    sfx = np.sum(fx,axis=0)
    if(sfx != 0):
        fx = np.divide(fx,sfx)
    
    # normalize fy
    sfy = np.sum(fy,axis=0)
    if(sfy != 0):
        fy = np.divide(fy,sfy)

    # Calculate the kernel
    kernel = np.dot(fx.T,np.dot(P,fy))

    return kernel
 
def sample_graphlets(k=5, delta=0.05, epsilon=0.05, a=-1):
    """ A function that samples graphlets, based on statistic parameters as proposed in :cite:`Shervashidze2009EfficientGK`.
    
    Parameters
    ----------
    k : int
        Graphlet dimension.

    delta : float, default=0.05
        Confidence level (typically 0.05 or 0.1).
    
    epsilon : float, default=0.05
        Precision level (typically 0.05 or 0.1).
        
    a : int
        Number of isomorphism classes of graphlets.
        If -1 the number is the maximum possible, from a database 1 until 9
        or else predicted through interpolation.
        
    Returns
    -------
    nsamples : int
        The number of samples in the graphlet feature space.
    graphlets : dict
        A dictionary of graphlets with keys from 0 ... nsamples-1.
    P : np.array
        A matrix that has zeros where the pair of indexes corresponds to isomorphic graphlets and 1 otherwise.
    graph_bins : dict
        A dictionary of bined graphlets, based on isomorphism classes.
    nbins : int
        The number of isomorphism classes.
    """
    # apply checks
    if k>10:
        warnings.warn('graphlets are too big - computation may be slow')
    elif k==1:
        warnings.warn('nonce graphlets')
    elif k<1:
        raise ValueError('k must be bigger than 1')
    
    if delta>1 or delta<0:
        raise ValueError('delta must be in the range (0,1)')
    
    if epsilon>1 or epsilon<0:
        raise ValueError('epsilon must be in the range (0,1)')
    
    
    fallback_map = {1: 1, 2: 2, 3: 4, 4: 8, 5: 19, 6: 53, 7: 209, 8: 1253, 9: 13599}
    
    if type(a) is not int:
        raise ValueError('a must be an integer')
    elif a==0:
        raise ValueError('a cannot be zero')
    elif a<-1:
        raise ValueError('negative a smaller than -1 have no meaning')
    
    if(a==-1):
        if(k>9):
            warnings.warn('warning for such size number of isomorphisms is not known - interpolation on know values will be used')
            # Use interpolations
            isomorphism_prediction = interp1d(list(fallback_map.keys()), list(fallback_map.values()), kind='cubic')
            a = isomorphism_prediction(k)
        else:
            a = fallback_map[k]
    
    # Calculate number of samples
    nsamples = math.ceil(2*(a*np.log10(2) + np.log10(1/delta))/(epsilon**2))

    # stores graphlets as pynauty graphs
    graphlets = dict()
    # stores the graphlet certificates
    graphlet_set = set()
    
    # transforms a matrix to a dictionary 
    to_edge_dict_binary = lambda x : matrix_to_dict(x, '==', 1, k, False)
    # different isomorphism bins
    graph_bins = dict()
    nbins = 0
    # produce all the binary sequences up to 2^(k-1)
    max_n = 2**(k-1)
    all_bin = {i: np.array(list(j)) for (i,j) in zip(range(0,max_n),itertools.product([0, 1], repeat=k-1))}
    for i in range(0,nsamples):
        gr = np.empty(shape = (k,k))
        cert = list()
        f = True        
        while f:
            # Calculates a symmetric random matrix
            # is calculated characteristics 
            # with at least line and row.
            for j in range(0,k-1):
                gr[j,j] = .0
                # pick randomly a non zero line
                line = all_bin[random.randrange(1,max_n)]
                # assign
                gr[0:j,j], gr[j,(j+1):k] = line[0:j], line[j:k-1]
                # Apply also for symmetric
                gr[(j+1):k,j], gr[j,0:j] = line[j:k-1],line[0:j]
                # calculate the certificate of the graph
                cert += list(line[j:k-1])
            certificate = str(cert)
            # checks if this graph has already existed
            f = certificate in graphlet_set
        graphlet_set.add(str(cert))
        graphlets[i] = pynauty.Graph(k, True, to_edge_dict_binary(gr))
        
        # add the graph to an isomorphism class
        if i==0:
            graph_bins[0] = [0]
            nbins+=1
        else:
            newbin = True
            for j in range(nbins):
                if pynauty.isomorphic(graphlets[graph_bins[j][0]],graphlets[i]):
                    newbin = False
                    graph_bins[j].append(i)
                    break
            if newbin:
                graph_bins[nbins] = [i]
                nbins+=1
                
    # Produce Pij Matrix:
    # Based on the idea that 
    # if Pij = 1 and Pjk = 1 then Pik=1
    P = np.zeros(shape = (nsamples,nsamples))
    for i in range(nbins):
        pair = list(graph_bins[i])
        for (i,j) in itertools.combinations(pair,2):
            P[i,j]=1

    # also Pij = Pji
    for i in range(0,nsamples):
        for j in range(i+1,nsamples):
            P[j,i]=P[j,i]

    return nsamples, graphlets, P, graph_bins, nbins
