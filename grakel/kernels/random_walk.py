""" This file contains the random walk kernel
    as defined by [Kashima et al., 2003; Gartner et al., 2003]
"""

import numpy as np

from scipy.linalg import solve_sylvester
from numpy.linalg import inv

from ..graph import graph
from ..tools import inv_dict

def random_walk(X, Y, lamda=0.1, method_type="sylvester"):
    """ This is a function that implements the random walk kernel.

        X, Y: to valid graph types i.e. edge dictionary or adjacency_mat
        algorithm_type: "dijkstra" or "floyd_warshall" or "auto" 
        (for the best concerning current format)
    """
    g_x = graph(X)
    g_y = graph(Y)

    return random_walk_inner(g_x, g_y,lamda, method_type)


def random_walk_inner(Gx, Gy,lamda=0.1, method_type="sylvester"):
    """ This is a function that implements the random walk kernel.

        Gx, Gy: are two graph type objects representing relevant graphs to 
        be compared. Their format is a s    uare array.
        lamda: the factor concerning summation
        method_type: "simple" [O(|V|^6)] or "sylvester" [O(|V|^3)]
    """
    Gx.desired_format("adjacency")
    Gy.desired_format("adjacency")
    X = Gx.adjacency_matrix
    Y = Gy.adjacency_matrix
    
    # calculate the product graph
    XY = np.kron(X,Y)

    if(method_type == "simple"):
        # algorithm presented in [Kashima et al., 2003; Gartner et al., 2003]
        # complexity of O(|V|^6)

        # XY is a square matrix
        s = XY.shape[0]
        I = np.identity(s)
        k = np.dot(np.dot(np.ones(s),inv(I - lamda*XY)).T,np.ones(shape=(s)))
    elif(method_type == "sylvester"):
        # algorithm presented in [Vishwanathan et al., 2006]
        # complexity of O(|V|^3)
        
        # Doing some initialization
        X_dimension = X.shape[0]
        Y_dimension = Y.shape[0]

        # For efficiency reasons multiply lambda
        # with the smallest e_{x,y} in dimension
        e_x = np.ones(shape=(X_dimension,1))
        e_y = np.ones(shape=(1,Y_dimension))

        # Prepare parameters for sylvester equation
        A = Y
        B = np.divide(inv(X.T), -lamda)
        C = -np.dot(e_x ,np.dot(e_y,B))
        
        R = solve_sylvester(A, B, C)
        
        # calculate kernel
        k = - np.sum(np.sum(R,axis=1),axis=0)
    else:
        pass
        # raise exception?
        # such method does not exist

    return k

