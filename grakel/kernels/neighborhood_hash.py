""" This file contains the neighborhood hashing kernel as defined in :cite:`Hido2009ALG`
"""
import itertools
import os

import numpy as np

from ..graph import graph
from ..tools import rotl, rotr

def neighborhood_hash(X, Y, Lx, Ly, nh_type='simple', R=3, bytes=2) :
    """ The neighborhood hashing kernel as proposed in :cite:`Hido2009ALG`

    arguments:
        - X,Y (valid graph format): the pair of graphs on which the kernel is applied
        - L{x,y} (dict): coresponding graph labels for nodes
        - R (int): the maximum number of neighborhood hash
        - nh_type (str): 'simple' or 'count-sensitive'
        - bytes (int): number of bytes for hashing of original labels
    returns:
        number. The kernel value
    """
    Gx = graph(X,Lx)
    Gy = graph(Y,Ly)
    return float(neighborhood_hash_matrix({0: Gx}, {0: Gy}, nh_type=nh_type, R=R, bytes=bytes)[0,0])

def neighborhood_hash_matrix(Graphs_x, Graphs_y=None, nh_type='simple', R=3, bytes=2):
    """ The calculate_similarity matrix function as defined :cite:`Hido2009ALG`

    arguments:
        - Graphs_{x,y} (dict): A dictionary of graph type objects that are going to be compared with keys from 0 ... to the dictionary length.
        - R (int): the maximum number of neighborhood hash
        - nh_type (str): 'simple' or 'count-sensitive'
        - bytes (int): number of bytes for hashing of original labels
    returns:
        np.array. The kernel matrix        
    """
    if R<=0:
        raise ValueError('R must be bigger than zero')
    
    if nh_type == 'simple':
        noc_f = False
        NH = lambda G: neighborhood_hash_simple(G)
    elif nh_type == 'count-sensitive':
        noc_f = True
        NH = lambda G: neighborhood_hash_count_sensitive(G)
    else:
        raise ValueError('unrecognised neighborhood hashing type')
        
    bytes = int(bytes)
    if bytes <= 0:
        raise ValueError('illegal number of bytes for hashing')
        
    h_x = len(Graphs_x.keys())
    Gs = dict()
    
    if Graphs_y==None:
        h_y = h_x
        g_iter = Graphs_x.values()
        ng = h_x
        pairs = [(i,j) for i in range(0,h_x) for j in range(i, h_x)]
        offset = 0
    else:
        h_y = len(Graphs_y.keys())
        g_iter = itertools.chain(Graphs_x.values(), Graphs_y.values())
        ng = h_x+h_y
        pairs = list(itertools.product(range(h_x, ng), range(0,h_x)))
        offset = h_x
    labels_hash_dict, labels_hash_set = dict(), set()
    for (i, g) in enumerate(g_iter):
        g.desired_format('adjacency')
        vertices = list(g.get_vertices())
        labels = hash_labels(g.get_labels(),labels_hash_dict, labels_hash_set, bytes)
        Gs[i] = tuple(radix_sort(vertices,labels))+({v: g.neighbors(v, purpose="adjacency") for v in vertices},)
        if noc_f:
            noc = dict()
            for k in labels.keys():
                if labels[k] not in noc:
                    noc[labels[k]] = 1
                else:
                    noc[labels[k]] += 1
            Gs[i] += (noc,)
        
    K = dict()
    S = np.zeros(shape=(h_x,h_y))
    for r in range(0,R):
        K = np.eye(h_x,h_y)
        for i in range(0, ng):
            Gs[i] = NH(Gs[i])
        for (i,j) in pairs:
            # targets - y - graph take the row
            K[i-offset,j] = nh_compare_labels(Gs[i], Gs[j])
        S = np.add(K,S)
        
    kernel_mat = np.divide(S,R)
    
    if Graphs_y is None:
        kernel_mat = np.triu(kernel_mat) + np.triu(kernel_mat, 1).T
        
    return kernel_mat
  
def hash_labels(labels, labels_hash_dict, labels_hash_set, bytes=2):
    """ A function that hashes existing labels to 16-bit integers without collisions and with consistency in same labels

    arguments:
        - Graphs_{x,y} (dict): A dictionary of graph type objects that are going to be compared with keys from 0 ... to the dictionary length.
        - labels (dict): labels for nodes
        - labels_hash_dict (dict): a hash table for labels
        - labels_hash_set (set): a set of all possible labels
        - nh_type (str): 'simple' or 'count-sensitive'
        - bytes (int): a positive integer denoting the size of hashes
    returns:
        dict. The new hashed labels for nodes
    """
    bytes = int(bytes)
    if bytes <= 0:
        raise ValueError('illegal number of bytes for hashing')
    
    new_labels = dict()
    for k in labels.keys():
        if labels[k] not in labels_hash_dict:
            f = True
            while f:
                r = int.from_bytes(os.urandom(bytes), 'little')
                f = r in labels_hash_set
            labels_hash_set.add(r)
            labels_hash_dict[labels[k]] = r
            new_labels[k] = r
        else:
            new_labels[k] = labels_hash_dict[labels[k]]
    return new_labels
            
def radix_sort(vertices, labels):
    """ sorts vertices based on labels
    
        arguments:
            - vertices (listable): listable of vertices
            - labels (dictionary): dictionary of labels for nodes
        returns:
            tuple. The sorted vertices based on labels, labels for nodes
    """
    
    return (sorted(list(vertices), key = lambda x: labels[x]),labels)
 
def neighborhood_hash_simple(G):
    """ Produces the (simple) neighborhood hash as defined in :cite:`Hido2009ALG`
    
    arguments:
        G (tuple): a tuple of three elements consisting of vertices sorted by labels, node label dictionary and edge dictionary.
    returns:
        tuple. a tuple of nodes, new_labels dictionary and edges
    """
    nodes, labels, edges = G
    new_labels = dict()
    for u in nodes:
        l = ROT(labels[u],1)
        for n in edges[u]:
            l ^= labels[n]
        new_labels[u] = l
    return (nodes, new_labels, edges)
    
def neighborhood_hash_count_sensitive(G):
    """ Produce the count sensitive neighborhood hash
        as defined in [Hido, Kashima et al., 2009]
    
    arguments:
        G (tuple): a tuple of three elements consisting of vertices sorted by labels, node label dictionary, edge dictionary and number of occurencies dictionary for labels
    returns:
        tuple. A tuple of three elements consisting of vertices sorted by labels, node label dictionary, edge dictionary and number of occurencies dictionary
    """
    nodes, labels, edges, noc= G
    new_labels = dict()
    new_noc = dict()
    for u in nodes:
        l = ROT(labels[u],1)
        for n in edges[u]:
            o = noc[labels[n]]
            l ^= ROT(labels[n] ^ o, o)
        if l not in new_noc:
            new_noc[l] = 1
        else:
            new_noc[l] += 1
        new_labels[u] = l
    return (nodes, new_labels, edges, new_noc)
    
def ROT(l,n):
    """ The rot operation for binary numbers
        
    arguments:
        - n (int): an integer
        - l (int): a valid number
    returns:
        (int). The result of a rot operation
    """
    if n < 0:
        return rotr(l,-n)
    elif n > 0:
        return rotl(l,n)
    else:
        return l
    
def nh_compare_labels(Gx, Gy):
    """ The compared labels function as defined in :cite:`Hido2009ALG`

    arguments:
        G_{x,y} (tuple): graph tuples of two elements, consisting of vertices sorted by labels and node, edge label dictionary

    returns:
        number. The kernel value
    """
    
    # get vertices
    vx, vy = iter(Gx[0]), iter(Gy[0])
    
    # get size of nodes
    nv_x, nv_y = len(Gx[0]), len(Gy[0])
    
    # get labels for nodes
    Lx, Ly = Gx[1], Gy[1]
       
    c = 0
    ui, uj = next(vx, None), next(vy, None)
    while (ui is not None) and (uj is not None):
        if Lx[ui]==Ly[uj]:
            c+=1
            ui = next(vx, None)
            uj = next(vy, None)
        elif Lx[ui]<Ly[uj]:
            ui = next(vx, None)
        else:
            uj = next(vy, None)
   
    return c/float(nv_x+nv_y-c)
   
   
