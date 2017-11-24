""" This file contains the standard
    and simple dirac kernel.
"""

from ..graph import graph
from ..tools import inv_dict

def dirac(X, Y, Lx, Ly):
    """ The simple dirac kernel for labelled graphs
        k(X,Y) = \sum_{v \in V_{1}}\sum_{v' \in V_{2}}dirac(l(v),l(v'))

        X,Y: Valid graph formats to be compared
        L{x,y}: The corresponding label dictionaries.
    """
    Gx = graph(X, Lx)
    Gy = graph(Y, Ly)
    return dirac_inner(Gx,Gy)

def dirac_inner(Gx,Gy):
    """ The simple dirac kernel for labelled graphs
        k(X,Y) = \sum_{v \in V_{1}}\sum_{v' \in V_{2}}dirac(l(v),l(v'))

        G_{x,y}: corresponding graph structures
    """
    Gx.desired_format("dictionary")
    Gy.desired_format("dictionary")
    
    # Calculate kernel
    linv_x = inv_dict(Gx.get_labels(purpose="dictionary"))
    linv_y = inv_dict(Gy.get_labels(purpose="dictionary"))
    
    kernel = 0
    for lx in linv_x:
        if lx in linv_y:
            kernel += len(linv_x[lx])*len(linv_y[lx])   
    
    return kernel
