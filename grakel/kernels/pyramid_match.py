""" This file contains the pyramid match kernel
    as defined in [1]:
    
    @inproceedings{
        Nikolentzos2017MatchingNE,
        title={Matching Node Embeddings for Graph Similarity},
        author={Giannis Nikolentzos and Polykarpos Meladianos and Michalis Vazirgiannis},
        booktitle={AAAI},
        year={2017}
    }
"""
import itertools

import numpy as np

from scipy.sparse.linalg import eigs
from scipy.sparse import csr_matrix

from ..graph import graph

def pyramid_match(X, Y, Lx, Ly, L=4, d=6) :
    """ The neighborhood hashing kernel as proposed
        in [Hido, Kashima, 2009]

        X,Y: Valid graph formats to be compared
        L{x,y}: labels for nodes for graphs X, Y
        
        --- Both labels should correspond to the 
        graph format given on X, Y.
        
        d: the dimension of the hypercube
        L: pyramid histogram level
    """
    Gx = graph(X,Lx)
    Gy = graph(Y,Ly)
    
    # If no labels are given the kernel is calculated without labels
    if Lx != {} and Ly != {}:
        with_labels = True
    else:
        with_labels = False
        
    return float(pyramid_match_matrix({0: Gx}, {0: Gy}, with_labels=with_labels, L=L, d=d)[0,0])

def pyramid_match_matrix(Graphs_x, Graphs_y=None, with_labels=True, L=4, d=6):
    """ The calculate_similarity matrix function as defined
        in [1], p. 4

        Graphs_{x,y}: A dictionary o graph type objects that are
                      going to be compared with keys from 0 to
                      the number of values.
        d: the dimension of the hypercube
        L: pyramid histogram level
    """

    # Input checks
    if type(with_labels) != bool:
        raise ValueError('with labels must be a boolean variable')
        
    if L<0:
        raise ValueError('L: the number of levels must be bigger equal to 0')

    if d<1:
        raise ValueError('d: hypercube dimension must be bigger than 1')
        
    # Parse input
    h_x = len(Graphs_x.keys())
    Gs = dict()
    
    if Graphs_y==None:
        h_y = h_x
        g_iter = list(Graphs_x.values())
        ng = h_x
        pairs = [(i,j) for i in range(0,h_x) for j in range(i, h_x)]
        offset = 0
    else:
        h_y = len(Graphs_y.keys())
        g_iter = list(itertools.chain(Graphs_x.values(), Graphs_y.values()))
        ng = h_x+h_y
        pairs = list(itertools.product(range(h_x, ng), range(0,h_x)))
        offset = h_x
    
    # hash labels
    num_labels = 0
    if with_labels:
        labels = set()
        for g in g_iter:
            g.desired_format('adjacency')
            labels |= set(g.get_labels(purpose='adjacency').values())
        num_labels = len(labels)
        labels = {l:i for (i,l) in enumerate(labels)}
    
    # calculate eigenvalues
    Us = []
    for g in g_iter:
        n = g.nv()
        if n==0:
            Us.append(np.zeros((1,d)))
        else:
            A = g.get_adjacency_matrix()
            if n > d+1:
                Lambda, U = eigs(csr_matrix(A, dtype=np.float), k=d, ncv=10*d)
                idx = Lambda.argsort()[::-1]
                U = U[:, idx]
            else:
                Lambda,U = np.linalg.eig(A)
                idx = Lambda.argsort()[::-1]   
                U = U[:,idx]
                U = U[:,:d]
            U = np.absolute(U)
            Us.append(U)
    
    # calculate histograms
    Hs = {}
    for (i,g) in enumerate(g_iter):
        n = g.nv()
        if n > 0:
            Hs[i] = []
            for j in range(L):
                l = 2**j
                if with_labels:
                    D = np.zeros((d*num_labels, l))
                else:
                    D = np.zeros((d, l))
                T = np.floor(Us[i]*l)
                T[np.where(T==l)] = l-1
                for p in range(Us[i].shape[0]):
                    if p >= n:
                        continue
                    for q in range(Us[i].shape[1]):
                        if with_labels:
                            D[labels[g.label(p, purpose='adjacency')]*d + q,int(T[p,q])] = D[labels[g.label(p, purpose='adjacency')]*d + q,int(T[p,q])] + 1
                        else:
                            D[q,int(T[p,q])] = D[q,int(T[p,q])] + 1

                Hs[i].append(D)
        
    # calculate kernel
    K = np.zeros(shape=(h_x,h_y))
    for (i,j) in pairs:
        if i in Hs and j in Hs:
            k = 0
            intersec = np.zeros(L)
            for p in range(L):
                # calculate histogram intersection
                intersec[p] = np.sum(np.minimum(Hs[i][p], Hs[j][p]))
                k += intersec[L-1]
                for p in range(L-1):
                    k += (1.0/(2**(L-p-1)))*(intersec[p]-intersec[p+1])
            K[i-offset,j] = k
        
    if Graphs_y is None:
        K = np.triu(K) + np.triu(K, 1).T
    
    return K
