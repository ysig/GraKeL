""" The propagation kernel as defined in :cite:`Neumann2015PropagationKE`.
"""
import itertools

import numpy as np

from ..graph import graph

np.random.seed(1235476565)

def propagation(X, Y, Lx, Ly, Tx=None, Ty=None, t_max=5, w=10, M="TV", base_kernel=(lambda x, y: sum(a*b for (a,b) in zip(x,y)))):
    """ The propagation kernel for fully labeled graphs :cite:`Neumann2015PropagationKE`: Algorithms 1, 3, p. 216, 221.

    Parameters
    ----------
    X,Y : *valid-graph-format*
        The pair of graphs on which the kernel is applied.
        
    L{x,y} : dict
        Corresponding graph labels for vertices.
        
    T{x,y} : np-array, default=None
        Transition matrices corresponding to Graphs X, Y (*if None default value is adjacency matrix*).
        
    t_max : int, default=5
        Maximum number of iterations.
        
    w : int, default=10
        Bin width.
        
    M : str, default="TV"
        The preserved distance metric (on local sensitive hashing):
            - "H": hellinger
            - "L1": l1-norm
            - "L2": l2-norm
            - "TV": total-variation
            
    base_kernel : function (list, list -> number), default=:math:`f(x,y)=\sum_{i} x_{i}*y_{i}`
        A base_kernel between two lists of numbers that outputs a number.
    
    Returns
    -------
    kernel : number
        The kernel value.
    """
    Gx = graph(X,Lx)
    Gy = graph(Y,Ly)
    return float(propagation_matrix({0: Gx}, {0: Gy}, t_max=t_max, T={0:Tx, 1:Ty}, w=w, M=M, base_kernel=base_kernel)[0,0])
    
        
def propagation_matrix(Graphs_x, Graphs_y=None, t_max=5, T=None, w=10, M="TV", base_kernel=(lambda x, y: sum(a*b for (a,b) in zip(x,y)))):
    """ The propagation kernel for fully labeled graphs :cite:`Neumann2015PropagationKE`: p. 216, 221.

    Parameters
    ----------
    Graphs_{x,y} : dict, default_y=None
        Enumerative dictionary of graph type objects with keys from 0 to the number of values.
        If value of Graphs_y is None the kernel matrix is computed between all pairs of Graphs_x
        where in another case the kernel_matrix rows correspond to elements of Graphs_y, and columns
        to the elements of Graphs_x.
        
    L{x,y} : dict
        Corresponding graph labels for vertices.
        
    t_max : int, default=5
        Maximum number of iterations.
        
    T{x,y} : np.array,  default=None
        Transition matrices corresponding to Graphs X, Y (*if None default value is adjacency matrix*).
        
    M : str, default="TV"
        The preserved distance metric (on local sensitive hashing):
            + "H": hellinger
            + "L1": l1-norm
            + "L2": l2-norm
            + "TV": total-variation
            
    w : int
        The bin width.
        
    base_kernel : function (list, list -> number), default=:math:`f(x,y)=\sum_{i} x_{i}*y_{i}`
        A base_kernel between two lists of numbers that outputs a number.
        
    Returns
    -------
    kernel_mat : np.array
        The kernel matrix.
    """
    if M not in ["H", "L1", "L2", "TV"]:
        raise ValueError('Invalid metric type')
    
    if t_max<0:
        raise ValueError('The number of iterations must be positive')
    
    if w<=0:
        raise ValueError('The bin width must be bigger than 0')
    
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
    
    # semantic type-checking to the transition matrices
    t_none = {i:True for i in range(0,ng)}
    if T is not None:
        if type(T) is not dict:
            raise ValueError('wrong type for T')
        elif sorted(list(T.keys())) != list(range(0,ng)):
            raise ValueError('wrong indexing for T')
        else:
            for i in range(0,ng):
                if T[i] is not None:
                    t_none[i] = False
                    if T[i].shape[0] != T.shape[1]:
                        raise ValueError('transition matrix must be a square matrix')
                    if T[i].shape[0] != g_iter[i].nv():
                        raise ValueError('propagation matrix must have the same dimension as number of vertices')
    else:
        T = dict()

    # initialize kernel matrix
    K = np.zeros(shape=(h_x, h_y))
    
    # calculate labels
    labels = set()
    for i in range(0,ng):
        labels |= set(g_iter[i].get_labels(purpose='any').values())
        
    # enumerate labels
    enum_labels = {l:i for (i,l) in enumerate(sorted(list(labels)))}
    
    # number of features
    D = len(enum_labels.keys())
    
    # make a matrix for all graphs that contains label vectors
    P = dict()
    for (i,g) in enumerate(g_iter):
        # calculate feature vectors
        p = np.zeros(shape=(g.nv(), D))
        for (j,v) in enumerate(sorted(list(g.get_vertices(purpose='any')))):
            p[j,enum_labels[g.label(v, purpose='any')]]=1
        P[i] = p
        
        # calculate adjacency matrices
        if t_none[i]:
            T[i] = g.get_adjacency_matrix()
        # Cast to binary array    
        T[i] = (T[i] > 0).astype(int)
        
    # feature vectors
    phi = dict()
    for t in range(0, t_max):
        # for each graph hash P and produce the feature vectors
        for i in range(0,ng):
            phi[i] = np.bincount(calculate_LSH(P[i], w, M)).tolist()
            
        # calculate the base kernel between all pairs of feature vectors
        for (i,j) in pairs:
            K[i-offset,j] += base_kernel(phi[i],phi[j])
        
        # calculate the Propagation matrix if needed    
        if t<t_max-1:
            for i in range(0,ng):
                P[i] = np.dot(T[i],P[i])
    
    if Graphs_y is None:
        K = np.triu(K) + np.triu(K, 1).T
        
    return K
    
def calculate_LSH(X, w, M):
    """ A function for calculating Local Sensitive Hashing needed for propagation kernels defined in :cite:`Neumann2015PropagationKE`, p.12.
        
    Parameters
    ----------
    X : np.array
        A float array of shape (N, D) with N vertices and D features.
        
    w : int
        Bin width.
        
    M : str
        The preserved distance metric:
           - "H": hellinger
           - "L1": l1-norm
           - "L2": l2-norm
           - "TV": total-variation

    Returns
    -------
    lsh : np.array.
        The local sensitive hash coresponding to each vertex.
    """
    if len(X.shape)!=2:
        raise ValueError('X must be an N x D array')
    D = X.shape[1]
    
    if M not in ["H", "L1", "L2", "TV"]:
        raise ValueError('Invalid metric type')
    
    if M is "H":
        X = np.sqrt(X)
    
    # simple normal
    u = np.random.randn(D)
    if M in ["TV", "L1"]:
        # cauchy
        u = np.divide(u, np.random.randn(D))

    # random offset
    b = w*np.random.rand()
    
    # hash
    h = np.floor((np.dot(X,u)+b)/w)
    return np.unique(h, return_inverse=True)[1]

