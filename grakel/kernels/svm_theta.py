""" This file contains the svm theta kernel
    as defined in [Johansson et al., 2014]
"""

import itertools

import numpy as np

from ..graph import graph

def svm_theta(X, Y, n_samples=50, subsets_size_range=(2,8), metric=(lambda x, y:x*y)):
    """ The svm theta kernel as proposed
        in [Johansson et al., 2014]

        X,Y: Valid graph formats to be compared
        metric: the applied metric between the lovasz numbers
        of the two graphs: number, number -> number
        n_samples: number of samples
        subsets_size_range: a touple of min, max size of the vertex
                            set of sampled subgraphs
    """
    Gx = graph(X)
    Gy = graph(Y)
    return svm_theta_inner(Gx, Gy, n_samples, subsets_size_range, metric)

def svm_theta_inner(Gx, Gy, n_samples=50, subsets_size_range=(2,8), metric=(lambda x, y:x*y)):
    """ The lovasz theta kernel as proposed
        in [Johansson et al., 2014]

        Gx, Gy: Graph type objects
        metric: the applied metric between the svm_theta numbers
        of the two graphs: number, number -> number
        n_samples: number of samples
        subsets_size_range: a touple of min, max size of the vertex
                            set of sampled subgraphs
    """
    
    Ldx = Gx.calculate_subgraph_samples_metric_dictionary("svm", n_samples=n_samples, subsets_size_range=subsets_size_range)
    Ldy = Gy.calculate_subgraph_samples_metric_dictionary("svm", n_samples=n_samples, subsets_size_range=subsets_size_range)
    
    kernel = 0
    for level in Ldx.keys():
        if level in Ldy:
            if bool(Ldx[level]) and bool(Ldy[level]):
                Z = len(Ldx[level])*len(Ldy[level])
                kernel += sum(metric(x,y) for (x,y) in itertools.product(Ldx[level],Ldy[level]))

    return kernel
