""" This file contains the svm theta kernel as defined in :cite:`Johansson2015LearningWS`
"""

import itertools

import numpy as np

from ..graph import graph

def svm_theta(X, Y, n_samples=50, subsets_size_range=(2,8), metric=(lambda x, y:x*y)):
    """ The svm theta kernel as proposed in :cite:`Johansson2015LearningWS`

    arguments:
        - X,Y (valid graph format): the pair of graphs on which the kernel is applied
        - metric (function: number, number -> number): the applied metric between the svm_theta numbers of the two graphs
        - n_samples (int): number of samples
        - subsets_size_range (tuple): (min, max) size of the vertex set of sampled subgraphs
    returns:
        number. The kernel value
    """
    Gx = graph(X)
    Gy = graph(Y)
    return svm_theta_inner(Gx, Gy, n_samples, subsets_size_range, metric)

def svm_theta_inner(Gx, Gy, n_samples=50, subsets_size_range=(2,8), metric=(lambda x, y:x*y)):
    """ The svm theta kernel as proposed in :cite:`Johansson2015LearningWS`

    arguments:
        - G_{x,y} (graph): the pair graphs on which the kernel is applied
        - metric (function: number, number -> number): the applied metric between the svm_theta numbers of the two graphs
        - n_samples (int): number of samples
        - subsets_size_range (tuple): (min, max) size of the vertex set of sampled subgraphs
    returns:
        number. The kernel value
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
