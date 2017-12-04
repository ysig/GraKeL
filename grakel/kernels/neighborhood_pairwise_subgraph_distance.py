""" This file contains the neighborhood subgraph pairwise distance kernel
    as defined in [Costa et al., 2010]
"""
import itertools

from ..graph import graph

def neighborhood_pairwise_subgraph_distance(X, Y, Lx, Ly, LEx, LEy, r=3, d=4):
    """ The neighborhood subgraph pairwise distance kernel
        as proposed in [Costa et al., 2010]

        X,Y: Valid graph formats to be compared
        L{x,y}: labels for nodes for graphs X, Y
        LE{x,y}: labels for edges for graphs X, Y
        
        --- Both labels should correspond to the 
        graph format given on X, Y.
    """
    Gx = graph(X,Lx,LEx)
    Gy = graph(Y,Ly,LEy)
    return neighborhood_pairwise_subgraph_distance_inner(Gx, Gy, r=r, d=d)

def neighborhood_pairwise_subgraph_distance_inner(Gx, Gy, r=3, d=4):
    """ The neighborhood subgraph pairwise distance kernel
        as proposed in [Costa et al., 2010]
        
        Gx, Gy: Graph type objects
        r: radius
        d: depth
    """
    if r<0:
        raise ValueError('r must be a positive integer')
    
    if d<0:
        raise ValueError('d must be a positive integer')
    
    Gx.desired_format('adjacency')
    Gy.desired_format('adjacency')
    
    Nx, Dx, Dx_pair = Gx.produce_neighborhoods(r, with_distances=True, d=d)
    Ny, Dy, Dy_pair = Gy.produce_neighborhoods(r, with_distances=True, d=d)
    
    Hx = hash_neighborhoods(Gx, Nx, Dx_pair, r)
    Hy = hash_neighborhoods(Gy, Ny, Dy_pair, r)
    
    kernel = 0
    
    for distance in range(0,d+1):
        if distance in Dx and distance in Dy:
            pairs = list(itertools.product(Dx[distance],Dy[distance]))
            for radius in range(0,r+1):
                krd = 0
                for ((A,B),(Ap,Bp)) in pairs:
                    if Hx[radius][A] == Hy[radius][B]:
                        krd+= int(Hx[radius][Ap] == Hy[radius][Bp])
                if len(pairs)>0:
                    # normalization by the number of pairs
                    kernel += float(krd)/len(pairs)
                    
                
                              
    return kernel
    
def hash_neighborhoods(G, N, D_pair, r, purpose = "adjacency"):
    """ A function that calculates the hash for all
        neighborhoods and all root nodes.
        
        G: graph type object of the original graph
        N: Dictionary with int keys as levels and dictionaries as values
           with symbol keys as root nodes and values a list of symbols for
           that correspond to this neighborhood
        D_pairs: a dictionary with keys as pairs of symbols corresponding to 
                 and values corresponding to element distances
    """
    H = {ra: dict() for ra in range(0,r+1)}
    for v in G.get_vertices(purpose):
        for radius in range(r, -1, -1):
            H[radius][v] = hash_graph(G.get_subgraph(N[radius][v]), D_pair)
    return H

def hash_graph(G, D, purpose='adjacency'):
    """ Make labels for hashing according to the proposed method
        and produce the graph hash needed for fast comparison
        
        G: graph type object of the original graph
        D: a dictionary with keys as pairs of symbols corresponding to 
           and values corresponding to element distances
    """
    encoding = ""
    s = str()
    
    # Make labels for vertices
    Lv = dict()
    glv = G.get_labels(label_type="vertex", purpose=purpose)
    for i in sorted(G.get_vertices(purpose=purpose)):
        l = sorted([(str(D[(i,j)])+str(',')+str(glv[j])) for j in G.get_vertices() if j!=i])
        encoding += str(l)+"."
        Lv[i] = l
    encoding = encoding[:-1]+":"
    
    # Expand to labels for edges
    Le = dict()
    gle = G.get_labels(label_type="edge", purpose=purpose)
    for (i,j) in sorted(list(G.get_edges(purpose=purpose))):
        Le[(i,j)] = str(Lv[i])+str(',')+str(Lv[j])+str(',')+str(gle[(i,j)])
        encoding += str(l)+"_"
    
    return APHash(encoding)

    
global bit_ones
bit_ones = dict()

def as_n_bit(num, n):
    n = int(n)
    
    if n<=0:
        raise ValueError('n must be an integer bigger than zero')
        
    if n not in bit_ones:
        x = int('0b'+n*'1',2)
        bit_ones[n] = x
    else:
        x = bit_ones[n]
        
    return x & num

def APHash(string):
    """ Arash Partov hashing as implemented in the original
        implementation of NPSDK.
    """
    hash_num = (int('0xAAAAAAAA',16))
    n = 32
    for s in list(string):
        if hash_num & 1 == 0:
            hash_num ^= as_n_bit((hash_num << 7),n) ^ as_n_bit(as_n_bit(ord(s),n) * as_n_bit((hash_num >> 3),n),n)
        else:
            hash_num ^= ~as_n_bit((as_n_bit((hash_num << 11),n) + (as_n_bit(ord(s),n) ^ as_n_bit((hash_num >> 5),2))),n)
            
    return hash_num
