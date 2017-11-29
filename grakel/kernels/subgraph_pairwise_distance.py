""" This file contains the subgraph pairwise distance kernel
    as defined in [Costa et al., 2010]
"""

from ..graph import graph

def subgraph_pairwise_distance(X, Y, Lx, Ly, LEx, LEy):
    """ The svm theta kernel as proposed
        in [Costa et al., 2010]

        X,Y: Valid graph formats to be compared
        L{x,y}: labels for nodes for graphs X, Y
        LE{x,y}: labels for edges for graphs X, Y
        
        --- Both labels should correspond to the 
        graph format given on X, Y.
    """
    Gx = graph(X,Lx,LEx)
    Gy = graph(Y,Ly,LEy)
    return svm_theta_inner(Gx, Gy)

def subgraph_pairwise_distance_inner(Gx, Gy, r, d):
    """ The lovasz theta kernel as proposed
        in [Costa et al., 2010]

        Gx, Gy: Graph type objects
    """
    Gx.desired_format('adjacency')
    Gy.desired_format('adjacency')
    
    Nx, Dx, Dx_pair = Gx.produce_neighbourhoods(r)
    Ny, Dy, Dy_pair = Gy.produce_neighbpurhoods(r)
    
    Sx = neighbourhoods_to_strings(Gx, Nx, Dx_pair, r):
    Sy = neighbourhoods_to_strings(Gy, Ny, Dy_pair, r):
    
    kernel = 0
    kratas gia ola ta zeugaria
    
    
    for distance in range(0,d+1):
        if distance in Dx[d] and distance in Dy[d]:
            for radius in range(0,r+1):
                for ((A,B),(Ap,Bp)) in itertools.product(Dx[d].keys(),Dy[d].keys()):
                    if bool(dirac_hash(Sx[r][A], Sy[r][B])):
                        kernel+= dirac_hash(Sx[r][Ap], Sy[r][Bp])
                        
   return kernel
   
def dirac_hash(Lx,Ly):
    return int(merkle_damgard(Lx)==merkle_damgard(Ly))
    
def neighbourhoods_to_labels(G, N, D_pair, r, purpose = "adjacency"):
    """ 
    """
    S = {ra: dict() for ra in range(0,r+1)}
    for v in G.get_vertex(purpose):
        subgraph = G
        for radius in range(r, -1, -1):
            subgraph.get_subgraph(N[radius][i])
            S[r][v] = make_labels_for_hashing(subgraph)
            

def make_labels_for_hashing(G):
    """ Make labels for hashing according to the proposed method
        G: a graph type object
    """
    s = str()
    # Make labels for nodes
    Ln = dict()
    for i in sorted(G.get_nodes()):
        l = sorted([(str(D[(i,j)])+str(',')+str(ln[j])) for j in range(0,n) if j!=i])
    Ln[i] = l
    
    # Expand to labels for edges
    Le = dict()
    for (i,j) in sorted(list(G.get_edges())):
        Le[(i,j)] = str(Ln[i])+str(',')+str(Ln[j])+str(',')+str(le[(i,j)])
    
    return (Ln, Le)
