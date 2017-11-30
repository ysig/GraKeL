""" This file contains the neighborhood hashing kernel
    as defined in [Hido, Kashima, 2009]
"""
import itertools

from ..graph import graph

def neighborhood_hash(X, Y, Lx, Ly, :
    """ The neighborhood hashing kernel as proposed
        in [Hido, Kashima, 2009]

        X,Y: Valid graph formats to be compared
        L{x,y}: labels for nodes for graphs X, Y
        
        --- Both labels should correspond to the 
        graph format given on X, Y.
    """
    Gx = graph(X,Lx)
    Gy = graph(Y,Ly)
    return neighborhood_pairwise_subgraph_distance_inner(Gx, Gy

def nh_calculate_similarity_matrix(Graph_list_x, Graph_list_y=None, R=3, nh_type='simple'):
    """ The calculate_similarity matrix function as defined
        in [Hido, Kashima, 2009]
        
        Graph_list_{x,y}: A list o graph type objects that
                          are going to be compared.
        R: the maximum number of neighborhood hash
        nh_type: 'simple' or 'count-sensitive'
    """
    if R<=0:
        raise ValueError('R must be bigger than zero')
    
    if nh_type is 'simple':
        NH = lambda G: neighbourhood_hash_simple(G)
    elif nh_type is 'count-sensitive':
        NH = lambda G: neighbourhood_hash_count_sensitive(G)
        
    h_x = len(Graph_list_x)
    Gs = dict()
    
    if Graph_list_y==None:
        h_y = h_x
        g_iter = Graph_list_x
        pairs = list(itertools.product(range(0,h_x), range(0, h_x))
        ng = h_x
    else:
        h_y = len(Graph_list_y)
        g_iter = itertools.chain(Graph_list_x, Graph_list_y)
        pairs = list(itertools.product(range(0,h_x), range(h_x, h_y))
        ng = h_x+h_y
        
    for (i, g) in enumerate(g_iter):
        g.desired_format('adjacency')
        Gs[i]= tuple(radix_sort(list(g.get_vertices()),g.get_labels()))+tuple(g.get_edges())
        
        
    K = dict()
    
    S = np.zeros(shape=(h_x,h_y))
    for r in range(0,R):
        K = np.eye(shape=(h_x,h_y))
        for i in range(0, ng):
            G[i] = NH(G[i])
        for (i,j) in pairs
            K[i,j] = nh_compare_labels(G[i], G[j])
        S = np.add(K,S)
    return np.divide(S,R)
            
def radix_sort(vertices, labels):
    return sorted(list(vertices), key = lambda x: labels[x])
    
def nh_compare_labels(Gx, Gy):
    """ The compared labels function as defined
        in [Hido, Kashima et al., 2009]

        G_{x,y}: Graph tuples as follows
    """
    # get size of nodes
    nv_x = Gx.nv()
    nv_y = Gy.nv()
    
    # get labels for nodes
    Lx = Gx.get_labels()
    Ly = Gy.get_labels()
    
    c, i, j = 0, 1, 1
    while (i<=nv_x) and (j<=nv_y):
        if Lx[i]==Ly[j]:
            c+=1
            i+=1
            j+=1
        elif Lx[i]<Ly[j]:
            i+=1
        else
            j+=1
   
   return c/float(nv_x+nv_y-c)
   
   
