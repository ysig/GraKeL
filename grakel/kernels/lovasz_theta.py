""" This file contains a set of functions implementing
    various graph kernel metrics.

"""
import itertools
import math
import random

import numpy as np

import pynauty

from numpy.linalg import inv
from ..graph import graph
from ..tools import inv_dict, matrix_to_dict, nested_dict_get, nested_dict_add

np.random.seed(238537)
random.seed(3456789)

def lovasz_theta(X, Y, Lx, Ly,metric=(lambda x, y:x*y)):
    """ The simple dirac kernel for labelled graphs
        k(X,Y) = \sum_{v \in V_{1}}\sum_{v' \in V_{2}}dirac(l(v),l(v'))

        X,Y: Valid graph formats to be compared
        L{x,y}: The corresponding label dictionaries.
    """
    Gx = graph(X, Lx)
    Gy = graph(Y, Ly)
    return lovasz_theta_inner(Gx,Gy,metric)

def lovasz_theta_inner(Gx,Gy,metric=(lambda x, y:x*y)):
    """ The simple dirac kernel for labelled graphs
        k(X,Y) = \sum_{v \in V_{1}}\sum_{v' \in V_{2}}dirac(l(v),l(v'))

        G_{x,y}: corresponding graph structures
    """
    Gx.desired_format("dictionary")
    Gy.desired_format("dictionary")
    
    # Calculate orthonormal matrix
    # Calculate lovasz theta value
    # Calculate kernel
    # For B subset of V
    # For C subset of V', |B|=|C|
    # k += 1/((n |B|) * (n' |B|)) * k(theta_B, theta_C)
    # Make these calculation accesible for subsets V, V'
    
