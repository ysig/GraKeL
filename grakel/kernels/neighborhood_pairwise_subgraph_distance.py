""" This file contains the neighborhood subgraph pairwise distance kernel as defined in :cite:`Costa2010FastNS`.
"""
import itertools

from ..graph import graph

def neighborhood_pairwise_subgraph_distance(X, Y, Lx, Ly, LEx, LEy, r=3, d=4):
    """ The neighborhood subgraph pairwise distance kernel as proposed in :cite:`Costa2010FastNS`.

    Parameters
    ----------
    X,Y : valid graph format.
        The pair of graphs on which the kernel is applied.
        
    L{x,y} : dict
        Corresponding graph labels for nodes.
        
    LE{x,y} : dict
        Corresponding graph labels for edges.
        
    r : int, default=3
        The maximum considered radius between vertices.
        
    d : int, default=4
        Neighborhood depth.
        
    Returns
    -------
    kernel : number
        The kernel value.
    """
    Gx = graph(X,Lx,LEx)
    Gy = graph(Y,Ly,LEy)
    return neighborhood_pairwise_subgraph_distance_inner(Gx, Gy, r=r, d=d)

def neighborhood_pairwise_subgraph_distance_inner(Gx, Gy, r=3, d=4):
    """ The neighborhood subgraph pairwise distance kernel as proposed in :cite:`Costa2010FastNS`.
        
    Parameters
    ----------
    G{x,y} : graph
        The pair of graphs on which the kernel is applied.
        
    r : int 
        The maximum considered radius between vertices.
        
    d : int
        Neighborhood depth.
        
    Returns
    -------
    kernel : number
        The kernel value.
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
    """ A function that calculates the hash for all neighborhoods and all root nodes.
    
    Parameters
    ----------
    G : graph
        The original graph.
        
    N : dict
        Neighborhoods that map levels (int) to dictionaries of root node symbols (keys) to list of vertex symbols, which correspond to the neighbors, that belong to this neighborhood.
        
    D_pairs : dict
        A dictionary that maps edges (tuple pairs of vertex symbols) to element distances (int - as produced from a BFS traversal).
    
    Returns
    -------
    H : dict
        The hashed neighborhoods as a 2-level dict from radious, vertex to the hashed values.
    """
    H = {ra: dict() for ra in range(0,r+1)}
    for v in G.get_vertices(purpose):
        for radius in range(r, -1, -1):
            H[radius][v] = hash_graph(G.get_subgraph(N[radius][v]), D_pair)
    return H

def hash_graph(G, D, purpose='adjacency'):
    """ Make labels for hashing according to the proposed method and produce the graph hash needed for fast comparison.
        
    Parameters
    ----------
    G : graph
        Original graph.
        
    D_pairs : dict
        A dictionary that maps edges (tuple pairs of vertex symbols) to element distances (int - as produced from a BFS traversal).
    
    Returns
    -------
    hash : int.
        The hash value for the given graph.
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
    """ A function that conserves the n-first digits of a given number.
        
    Parameters
    ----------
    num : int
        The binary number.

    n : int
        The number of bits.
        
    Returns
    -------
    num_n : int
        The n first numbers of the integer num.
    """
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
    """ Arash Partov hashing as implemented in the original implementation of NPSDK.

    Parameters
    ----------
    string : str
        The string that is going to be hashed.
        
    Returns
    -------
    AP_hash : int
        The AP hash value for the given text.
    """
    hash_num = (int('0xAAAAAAAA',16))
    n = 32
    for s in list(string):
        if hash_num & 1 == 0:
            hash_num ^= as_n_bit((hash_num << 7),n) ^ as_n_bit(as_n_bit(ord(s),n) * as_n_bit((hash_num >> 3),n),n)
        else:
            hash_num ^= ~as_n_bit((as_n_bit((hash_num << 11),n) + (as_n_bit(ord(s),n) ^ as_n_bit((hash_num >> 5),2))),n)
            
    return hash_num
