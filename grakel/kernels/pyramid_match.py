""" This file contains the pyramid match kernel as defined in :cite:`Nikolentzos2017MatchingNE`.
"""
import itertools

import numpy as np

from scipy.sparse.linalg import eigs
from scipy.sparse import csr_matrix

from ..graph import graph

def pyramid_match(X, Y, Lx, Ly, L=4, d=6) :
    """ The pyramid match kernel as defined in :cite:`Nikolentzos2017MatchingNE`.

    Parameters
    ----------
    X,Y : *valid-graph-format*
        The pair of graphs on which the kernel is applied.
        
    L{x,y} : dict
        Corresponding graph labels for vertices.
        
    L : int, default=4
        Pyramid histogram level.
        
    d : int, default=6
        The dimension of the hypercube.

    Returns
    -------
    kernel : number
        The kernel value.
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
    """ The pyramid match kernel function as defined in :cite:`Nikolentzos2017MatchingNE`.

    Parameters
    ----------
    Graphs_{x,y} : dict, default_y=None
        Enumerative dictionary of graph type objects with keys from 0 to the number of values.
        If value of Graphs_y is None the kernel matrix is computed between all pairs of Graphs_x
        where in another case the kernel_matrix rows correspond to elements of Graphs_y, and columns
        to the elements of Graphs_x.
        
    with_labels : bool, default=True
        A flag that determines if the kernel computation will consider labels.
        
    L : int
        The pyramid histogram level.
        
    d : int
        The dimension of the hypercube.
        
    Returns
    -------
    kernel_mat : np.array
        The kernel matrix.
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
        pairs = [(i,j) for i in range(0, h_x) for j in range(i, h_x)]
        offset = 0
    else:
        h_y = len(Graphs_y.keys())
        g_iter = list(itertools.chain(Graphs_x.values(), Graphs_y.values()))
        ng = h_x+h_y
        pairs = list(itertools.product(range(h_x, ng), range(0, h_x)))
        offset = h_x
    
    # Map labels to values between 0 and L-1 where L is the number of distinct labels
    num_labels = 0
    if with_labels:
        labels = set()
        for g in g_iter:
            g.desired_format('adjacency')
            labels |= set(g.get_labels(purpose='adjacency').values())
        num_labels = len(labels)
        labels = {l:i for (i,l) in enumerate(labels)}
    
    # Embed vertices into the d-dimensional space
    Us = []
    for g in g_iter:
        n = g.nv()
        if n==0:
            Us.append(np.zeros((1,d)))
        else:
            # Perform eigenvalue decomposition. Rows of matrix U correspond to vertex representations
            A = g.get_adjacency_matrix()
            if n > d+1:
                # If size of graph smaller than d, pad with zeros
                Lambda, U = eigs(csr_matrix(A, dtype=np.float), k=d, ncv=10*d)
                idx = Lambda.argsort()[::-1]
                U = U[:, idx]
            else:
                Lambda,U = np.linalg.eig(A)
                idx = Lambda.argsort()[::-1]   
                U = U[:,idx]
                U = U[:,:d]
            # Replace all components by their absolute values
            U = np.absolute(U)
            Us.append(U)
    
    # Calculate histograms
    Hs = {}
    for (i,g) in enumerate(g_iter):
        n = g.nv()
        if n > 0:
            Hs[i] = []
            for j in range(L):
                # Number of cells along each dimension at level j
                l = 2**j
                if with_labels:
                    # To store the number of vertices that are assigned a
                    # specific label and lie in each of the 2^j cells of
                    # each dimension at level j
                    D = np.zeros((d*num_labels, l))
                else:
                    # Determines the cells in which each vertex lies
                    # along each dimension since nodes lie in the unit
                    # hypercube in R^d 
                    D = np.zeros((d, l))
                T = np.floor(Us[i]*l)
                T[np.where(T==l)] = l-1
                for p in range(Us[i].shape[0]):
                    if p >= n:
                        continue
                    for q in range(Us[i].shape[1]):
                        # Identify the cell into which the i-th vertex 
                        # lies and increase its value by 1
                        if with_labels:
                            D[labels[g.label(p, purpose='adjacency')]*d + q,int(T[p,q])] = D[labels[g.label(p, purpose='adjacency')]*d + q,int(T[p,q])] + 1
                        else:
                            D[q,int(T[p,q])] = D[q,int(T[p,q])] + 1

                Hs[i].append(D)
        
    # Calculate kernel
    K = np.zeros(shape=(h_y,h_x))
    
    for (i,j) in pairs:
        if i in Hs and j in Hs:
            k = 0
            intersec = np.zeros(L)
            for p in range(L):
                # Calculate histogram intersection (equation 6 in :cite:`Nikolentzos2017MatchingNE`)
                intersec[p] = np.sum(np.minimum(Hs[i][p], Hs[j][p]))
                k += intersec[L-1]
                for p in range(L-1):
                    # Computes the new matches that occur at level p.
                    # These matches weight less than those that occur at
                    # higher levels (e.g. p+1 level)
                    k += (1.0/(2**(L-p-1)))*(intersec[p]-intersec[p+1])
            K[i-offset,j] = k

    if Graphs_y is None:
        K = np.triu(K) + np.triu(K, 1).T
    
    return K
