""" This file contains the svm theta kernel as defined in :cite:`Johansson2015LearningWS`.
"""

import itertools

import numpy as np

from ..graph import graph

def svm_theta(X, Y, n_samples=50, subsets_size_range=(2,8), metric=(lambda x, y:x*y)):
    """ The svm theta kernel as proposed in :cite:`Johansson2015LearningWS`.

    Parameters
    ----------
    X,Y : *valid-graph-format*
        The pair of graphs on which the kernel is applied.
        
    n_samples : int, default=50
        Number of samples.
        
    subsets_size_range : tuple, len=2, default=(2,8)
        (min, max) size of the vertex set of sampled subgraphs.
        
    metric : function (number, number -> number), default=:math:`f(x,y)=x*y`
        The applied metric between the svm_theta numbers of the two graphs.
    
    Returns
    -------
    kernel : number
        The kernel value.
    """
    Gx = graph(X)
    Gy = graph(Y)
    return svm_theta_inner(Gx, Gy, n_samples, subsets_size_range, metric)

def svm_theta_inner(Gx, Gy, n_samples=50, subsets_size_range=(2,8), metric=(lambda x, y:x*y)):
    """ The svm theta kernel as proposed in :cite:`Johansson2015LearningWS`.

    Parameters
    ----------
    G_{x,y} : graph
        The pair of graphs on which the kernel is applied.
        
    n_samples : int, default=50
        Number of samples.
        
    subsets_size_range : tuple, len=2, default=(2,8)
        (min, max) size of the vertex set of sampled subgraphs.
        
    metric : function (number, number -> number), default=:math:`f(x,y)=x*y`
        The applied metric between the svm_theta numbers of the two graphs.
    
    Returns
    -------
    kernel : number
        The kernel value.
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
