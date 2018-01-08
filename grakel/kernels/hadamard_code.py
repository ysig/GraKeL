""" This file contains the hadamard code kernel as defined in :cite:`icpram16`.
"""
import itertools
import os
import math

import numpy as np

from scipy.linalg import hadamard

from ..graph import graph
from ..tools import rotl, rotr

def hadamard_code(X, Y, Lx, Ly, base_kernel, niter=5, hc_type='simple', **kwargs) :
    """ The hadamard code kernel as proposed in :cite:`icpram16`.

    Parameters
    ----------
    X,Y : *valid-graph-format*
        The pair of graphs on which the kernel is applied.
        
    L{x,y} : dict
        Corresponding graph labels for vertices.
    
    base_kernel : function (graph, graph -> number)
        A valid pairwise base_kernel for graphs.
    
    hc_type : str, valid_inputs={"simple", "shortened"}, default="simple"
        The hadamard code kernel type as defined in :cite:`icpram16`.
    
    rho : int, condition_of_appearance: hc_type=="shortened", default=-1
        The size of each single bit arrays. If -1 is chosen r is calculated as
        the biggest possible that satisfies an equal division.

    L : int, condition_of_appearance: hc_type=="shortened", default=4
        The number of bytes to store the bitarray of each label.

    niter : int, default=5
        The number of iterations.
    
    Returns
    -------
    kernel : number
        The kernel value.
    """
    Gx = graph(X,Lx)
    Gy = graph(Y,Ly)
    bkm = lambda x, y: np.array([[base_kernel(x[0], y[0])]])
    return float(hadamard_code_matrix({0: Gx}, bkm, {0: Gy}, niter=5, hc_type=hc_type, **kwargs)[0,0])

def hadamard_code_matrix(Graphs_x, base_kernel_matrix, Graphs_y=None, niter=5, hc_type='simple', **kwargs):
    """ The hadamard code kernel for matrices as derived from :cite:`icpram16`.

    Parameters
    ----------
    Graphs_{x,y} : dict, default_y=None
        Enumerative dictionary of graph type objects with keys from 0 to the number of values.
        If value of Graphs_y is None the kernel matrix is computed between all pairs of Graphs_x
        where in another case the kernel_matrix rows correspond to elements of Graphs_y, and columns66ertfgv 
        to the elements of Graphs_x.

    base_kernel_matrix : function (dict(graph), dict(graph) -> np.array)
        Rows of the np.array should correspond to the second dictionary of graphs and cols to the first.
        If the second graph is None the kernel should be computed upon itself.

    hc_type : str, valid_inputs={"simple", "shortened"}, default="simple"
        The hadamard code kernel type as defined in :cite:`icpram16`.
    
    rho : int, condition_of_appearance: hc_type=="shortened", default=-1
        The size of each single bit arrays. If -1 is chosen r is calculated as
        the biggest possible that satisfies an equal division.

    L : int, condition_of_appearance: hc_type=="shortened", default=4
        The number of bytes to store the bitarray of each label

    niter : int, default=5
        The number of iterations.
    
    Returns
    -------
    kernel_matrix : np.array
        The kernel matrix. If the Graphs_y is not None the rows correspond to Graphs_y
        and the cols to the Graphs_x (based on the given order).
    """

    # parmatrise and extract more arguments based on the input of hc_type
    if hc_type == 'simple':
        shortened = False
        add = np.add
        get = lambda A, i: A[i,:]
    elif hc_type == 'shortened':
        # extract rho
        rho = kwargs.get('rho', -1)
        if rho <= 0 and rho != -1:
            raise ValueError('rho must be bigger than zero or equal to 1')

        # extract L
        L = kwargs.get('L', 4)
        if L <= 0:
            raise ValueError('L must be a positive integer as it corresponds to the number of bits')
        L = L*8

        # initialise addition labels and get for an internal use
        shortened = True
        q = int('0b'+(L*'1'), 2)
        add = lambda x, y: (x + y) & q
        get = lambda A, i: A[i]
    else:
        raise ValueError('unrecognised hadamard code kernel type')

    if niter<1:
        raise ValueError('niter must be an integer bigger than zero')
     
    ng_a = len(Graphs_x.keys())
    Gs = dict()
    
    if Graphs_x == None:
        g_iter = [Graphs_x[i] for i in range(0, ng_a)]
        ng = ng_a
        ng_b = ng_a
    else:
        ng_b = len(Graphs_y.keys())
        g_iter = itertools.chain([Graphs_x[i] for i in range(0, ng_a)], [Graphs_y[i] for i in range(0, ng_b)])
        ng = ng_a + ng_b
        
    # Enumerate all labels and count their quantity
    nl = 0
    labels_enum = dict()
    for g in g_iter:
        for l in g.get_labels(purpose='any'):
            if l not in labels_enum:
                labels_enum[l] = nl
                nl+=1

    # Calculate the hadamard matrix
    ord_H = 2**(math.ceil(math.log2(nl)))
    H = hadamard(ord_H)
    
    if shortened:
        if rho == -1:
            rho = int(L/(ord_H-1))
            if rho < 0:
                raise ValueError('The default calculated rho is too small, raise more L ~ L must approach the number of labels ('+str(ord_H))

        elif rho*(ord_H-1)> L:
            raise ValueError('rho*(ord_H-1)='+str(rho)+'*('+str(ord_H)+'-1) > L='+str(L))
        bit_array = dict()
        bool_H = (H > 0).astype(bool)
        for i in range(ord_H):
            bit_array[i] = (L - rho*(ord_H-1))*str(int(bool_H[i,0]))
            for j in range(1, ord_H):
                bit_array[i] = bit_array[i] + rho*str(int(bool_H[i,j]))
            bit_array[i] = int('0b'+bit_array[i], 2)
        Labeling = bit_array
    else:
        Labeling = H

    # Intial labeling of vertices based on their
    # corresponding Hadamard code (i-th row of the Hadamard matrix)
    # where i is the i-th label on enumeration
    orig_labels = dict()
    for (i,g) in enumerate(g_iter):
        orig_labels[i] = g.get_labels(purpose='any')
        new_labels = dict()
        for v in orig_labels[i].keys():
            new_labels[v] = get(Labeling, labels_enum[ol[i][v]])
        g.relabel(nl, purpose='any')
    
    # Initialize kernel matrix to define shape
    kernel = np.zeros(shape=(ng_b, ng_a))
    
    # Add the zero iteration element
    kernel += base_kernel_matrix(Graphs_x, Graphs_y)
    
    # Main
    for i in range(niter):
        for g in g_iter:
            # Find unique labels and sort them for both graphs
            # Keep for each node the temporary
            old_labels = g.get_labels(purpose='any')
            new_labels = dict()
            for v in g.get_vertices(purpose='any'):
                new_labels[v] = old_labels[v]
                for q in g.get_neighbors(v):
                    new_labels[v] = add(new_labels[v], old_labels[v])
            g.relabel(new_labels)
            
        # calculate kernel 
        kernel += base_kernel_matrix(Graphs_x, Graphs_y)
        
    # Restore original labels  
    for g in g_iter:
        Gs[i].relabel(orig_labels[i])
    
    if Graphs_y is None:
        kernel = np.triu(kernel) + np.triu(kernel, 1).T

    return kernel
