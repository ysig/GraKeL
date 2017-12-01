""" This file contains the neighborhood hashing kernel
    as defined in [Hido, Kashima, 2009]
"""
import itertools
import os

import numpy as np

from ..graph import graph
from ..tools import rotl, rotr

def neighborhood_hash(X, Y, Lx, Ly, nh_type='simple', R=3, bytes=2) :
    """ The neighborhood hashing kernel as proposed
        in [Hido, Kashima, 2009]

        X,Y: Valid graph formats to be compared
        L{x,y}: labels for nodes for graphs X, Y
        
        --- Both labels should correspond to the 
        graph format given on X, Y.
        
        R: the maximum number of neighborhood hash
        nh_type: 'simple' or 'count-sensitive'
        bytes: number of bytes for hashing of original labels
    """
    Gx = graph(X,Lx)
    Gy = graph(Y,Ly)
    return float(neighborhood_hash_similarity_matrix([Gx], [Gy], nh_type=nh_type, R=R, bytes=bytes)[0,0])

def neighborhood_hash_similarity_matrix(Graph_list_x, Graph_list_y=None, nh_type='simple', R=3, bytes = 2):
    """ The calculate_similarity matrix function as defined
        in [Hido, Kashima, 2009]
        
        Graph_list_{x,y}: A list o graph type objects that
                          are going to be compared.
        R: the maximum number of neighborhood hash
        nh_type: 'simple' or 'count-sensitive'
        bytes: number of bytes for hashing of original labels
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
        
    h_x = len(Graph_list_x)
    Gs = dict()
    
    if Graph_list_y==None:
        h_y = h_x
        g_iter = Graph_list_x
        pairs = list(itertools.product(range(0,h_x), range(0, h_x)))
        ng = h_x
    else:
        h_y = len(Graph_list_y)
        g_iter = itertools.chain(Graph_list_x, Graph_list_y)
        pairs = list(itertools.product(range(0,h_x), range(h_x, h_y)))
        ng = h_x+h_y
        
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
            K[i,j] = nh_compare_labels(Gs[i], Gs[j])
        S = np.add(K,S)
    return np.divide(S,R)
  
def hash_labels(labels, labels_hash_dict, labels_hash_set, bytes=2):
    """ A function that hashes existing labels to 16-bit
        integers without collisions and with consistency 
        in same labels
        
        bytes: a positive integer denoting the hashing size
        labels: dictionary of labels
        labels_hash_dict: labels_hash_dictionary
        labels_hash_set: labels_hash_set
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
    return (sorted(list(vertices), key = lambda x: labels[x]),labels)
 
def neighborhood_hash_simple(G):
    """ Produces the (simple) neighborhood hash
        as defined in [Hido, Kashima et al., 2009]
        
        G: a tuple of three elements consisting of
           vertices sorted by labels, label dictionary
           for nodes and edge dictionary.
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
    
        G: a tuple of three elements consisting of
           vertices sorted by labels, label dictionary
           for nodes, edge dictionary and number of 
           occurencies for labels
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
        n: an integer
        l: a valid number
    """
    if n < 0:
        return rotr(l,-n)
    elif n > 0:
        return rotl(l,n)
    else:
        return l
    
def nh_compare_labels(Gx, Gy):
    """ The compared labels function as defined
        in [Hido, Kashima et al., 2009]

        G_{x,y}: Graph  tuples of three elements 
           consisting of vertices sorted by labels,
           label dictionary for nodes and edges.
    """
    # get vertices
    vx = iter(Gx[0])
    vy = iter(Gy[0])
    
    # get size of nodes
    nv_x = len(Gx[0])
    nv_y = len(Gy[0])
    
    # get labels for nodes
    Lx = Gx[1]
    Ly = Gy[1]
       
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
   
   
