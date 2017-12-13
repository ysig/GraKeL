""" This file contains the sugraph mathing kernel as defined by :cite:`Kriege2012SubgraphMK`
"""

from ..graph import graph
from ..tools import inv_dict

global kv_default, ke_default, lw_default


# Define defaul vertex, edge and lambda weight functions

kv_default = lambda x, y, Lx, Ly: 1 if Lx[x]==Ly[y] else 0

def ke_default(x, y, Ex, Ey, Lex, Ley):
    cond_a, cond_b = x in Ex, x in Ey
    if (cond_a and cond_b and Lex[x]==Ley[y]) or ((not cond_a) and (not cond_b)):
        return 1
    else:
        return 0

lw_default = lambda x: int(bool(x))

def subgraph_matching(X, Y, Lx, Ly, Lex, Ley, kv=kv_default, ke=ke_default, lw=lw_default):
    """ The subgraph matching kernel as proposed in :cite:`Kriege2012SubgraphMK`
        
    attributes:
        - X,Y (valid graph format): the pair of graphs on which the kernel is applied
        - L{x,y} (dict): coresponding graph labels for nodes
        - kv (function: (symbol, symbol, dict, dict) -> number): the kernel function for two nodes based on label dictionaries
        - ke (function: (tuple, tuple, dict, dict, dict, dict) -> number): the kernel function for edges (tuples) of two graphs based on their edges (dict 1,2) and edge labels (dict 3,4)          
        - lw (function: set -> number) : a lambda weight function for cliques
    returns:
        number. The kernel value 
    """
    Gx = graph(X, Lx, Ley)
    Gy = graph(Y, Ly, Ley)
    return subgraph_matching_inner(Gx, Gy, kv=kv, ke=ke, lw=lw)

def subgraph_matching_inner(Gx, Gy, kv=kv_default, ke=ke_default, lw=lw_default):
    """ The subgraph matching kernel as proposed in :cite:`Kriege2012SubgraphMK`
    
    attributes:
        - G{x,y} (graph): the pair graphs on which the kernel is applied
        - kv (function: (symbol, symbol, dict, dict) -> number): the kernel function for two nodes based on label dictionaries
        - ke (function: (tuple, tuple, dict, dict, dict, dict) -> number): the kernel function for edges (tuples) of two graphs based on their edges (dict 1,2) and edge labels (dict 3,4)          
        - lw (function: set -> number) : a lambda weight function for cliques 
    return:
        number. the kernel value
    """
    # Calculate product graph
    Vp, Ep, c = weighted_product_graph(Gx, Gy, kv=kv, ke=ke)
    
    # Initialize values
    C = set()
    value = [.0]
    
    subgraph_matching_core(1, C, Vp, Ep, c, lw, value)
    return value[0]
    
def subgraph_matching_core(w, C, P, Ep, c, lw, value):
    """ The core algorithm of the subgraph matching kernel as proposed in :cite:`Kriege2012SubgraphMK` (p.5, Algorithm 1)

    attributes:
        - G{x,y} (graph): the pair graphs on which the kernel is applied
        - w (number): weight
        - C (set): clique
        - P (set): current vertices
        - c (dict): costs for vertices, edges
        - lw (function: set -> number) : a lambda weight function for cliques 
        - value (list): one element mutable list used as global variable for the recursion of the algorithm
    returns:
        None
    """
    while len(P)>0:
        v = P.pop()
        w = w*c[v]
        
        for u in list(C):
            w = w*c[(u,v)]
        C.add(v)
        value[0] += w*lw(C)
        
        # Create the new P
        Pp = P.copy()
        Pp.add(v)
        Pp &= Ep[v]
        
        # apply subgraph matching for the new clique
        subgraph_matching_core(w, C, Pp, Ep, c, lw, value)
        
def weighted_product_graph(Gx, Gy, kv, ke):
    """ Calculates the weighted product graph as defined in :cite:`Kriege2012SubgraphMK` (p.5, Definition 5)
        
    attributes:
        - G{x,y} (graph): the pair graphs on which the kernel is applied
        - kv (function: (symbol, symbol, dict, dict) -> number): the kernel function for two nodes based on label dictionaries
        - ke (function: (tuple, tuple, dict, dict, dict, dict) -> number): the kernel function for edges (tuples) of two graphs based on their edges (dict 1,2) and edge labels (dict 3,4)
    returns:
        - set. Vertices
        - dict. Edges
        - dict. A dictionary of cost for vertices and edges
    """
    Gx.desired_format("all")
    Gy.desired_format("all")
    
    Lx = Gx.get_labels()
    Ly = Gy.get_labels()
    Lex = Gx.get_labels("edge")
    Ley = Gy.get_labels("edge")
    Ex = Gx.edge_dictionary
    Ey = Gy.edge_dictionary
        
    Kv = lambda x, y: kv(x,y,Lx,Ly)
    Ke = lambda x, y: ke(x,y,Ex,Ey,Lex,Ley)

    # initialise cost function
    c = dict()
    # Calculate valid vertices
    Vp = set()
    Ep = dict()
    
    # calculate product graph vertex set
    for i in range(0, Gx.n):
        for j in range(0, Gy.n):
            if(i!=j):
                value = Kv(i,j)
                if(value>0):
                    # add to vertex set
                    Vp.add((i,j))
                    # add cost
                    c[(i,j)]= value
                    # initialise an empty set for neighbors
                    if (i,j) not in Ep:
                        Ep[(i,j)] = set()
                else:
                    c[(i,j)] = 0
            else:
                c[(i,j)] = 0
            
                    
    
    np = len(Vp)
    # calculate product graph valid edges
    for v in list(Vp):
        for w in list(Vp):
            if(v[0]!=w[0] and v[1]!=w[1]):
                value = Ke(v,w)
                if(value>0):
                    # add edge
                    Ep[v].add(w)
                    # store value
                    c[(v,w)] = value
                else:
                    c[(v,w)] = 0
            else:
                c[(v,w)] = 0
    return Vp, Ep, c
