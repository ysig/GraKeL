""" The propagation kernel as defined in [1]:
' Neumann M.; Garnett R., P.A..: Propagation kernels: efficient graph kernels from propagated information. In: Springer. Mach Learn (2016) p.209–245'
"""
import itertools

import numpy as np

from ..graph import graph

np.random.seed(1235476565)

def calculate_LSH(X, w, M):
    """ A function for calculating
        Local Sensitive Hashing needed
        for propagation kernels defined on [1], p.12
        
        X: A float numpy array of shape (N, D)
        w: Bin width
        M: The preserved distance metric
           "H": hellinger
           "L1": l1-norm
           "L2": l2-norm
           "TV": total-variation
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

def propagation(X, Y, Lx, Ly, Tx=None, Ty=None, t_max=5, w=10, M="TV", base_kernel=(lambda x, y: sum(a*b for (a,b) in zip(x,y)))):
    """ The propagation kernel for fully labeled graphs
        [1]: Algorithms 1, 3, p. 216, 221

        X,Y: Valid graph formats to be compared
        L{x,y}: labels for nodes for graphs X, Y
        T{x,y}: np-arrays corresponding to Graphs X, Y
                (if None default value is adjacency matrix)
        
        --- Both labels should correspond to the 
        graph format given on X, Y.
        
        t_max: maximum number of iterations
                   
        M: The preserved distance metric (on local sensitive hashing)
           "H": hellinger
           "L1": l1-norm
           "L2": l2-norm
           "TV": total-variation
           
        w: bin width
        
        base_kernel: a base_kernel between two lists of numbers that outputs a number
    """
    Gx = graph(X,Lx)
    Gy = graph(Y,Ly)
    return float(propagation_matrix({0: Gx}, {0: Gy}, t_max=t_max, T={0:Tx, 1:Ty}, w=w, M=M, base_kernel=base_kernel)[0,0])
    
        
def propagation_matrix(Graphs_x, Graphs_y=None, t_max=5, T=None, w=10, M="TV", base_kernel=(lambda x, y: sum(a*b for (a,b) in zip(x,y)))):
    """ The propagation kernel for fully labeled graphs
        [1]: Algorithms 1, 3 p. 216, 221

        Graphs_{x,y}: A dictionary o graph type objects that are
                      going to be compared with keys from 0 to
                      the number of values.
        
        t_max: maximum number of iterations
        
        T: a dictionary of np.arrays corresponding to Graphs_x, Grpahs_y
           in sequence keys from 0 ... len(Graphs_x) ... len(Graphs_x)+len(Graphs_y)
           (if None default value is adjacency matrix)
           
        M: The preserved distance metric (on local sensitive hashing)
           "H": hellinger
           "L1": l1-norm
           "L2": l2-norm
           "TV": total-variation
           
        w: bin width
        
        base_kernel: a base_kernel between two lists of numbers that outputs a number
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
