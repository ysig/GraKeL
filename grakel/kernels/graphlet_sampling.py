"""The graphlet sampling kernel :cite:`Shervashidze2009EfficientGK`."""

import itertools
import math
import random
import warnings

import numpy as np

import pynauty

from scipy.interpolate import interp1d

from grakel.graph import graph
from grakel.tools import matrix_to_dict

random.seed(15487103)


def graphlet_sampling(X, Y, k=5, **kargs):
    """Graphlet sampling kernel :cite:`Shervashidze2009EfficientGK`.

    Parameters
    ----------
    X,Y : valid graph format
        The pair of graphs on which the kernel is applied.

    k : int, default=5
        The dimension of the given graphlets.

    delta : float, default=0.05
        Confidence level (typically 0.05 or 0.1).
        For calculation of the number of samples achieving the certain bound.
        n_samples argument must not be provided and for initialising the
        default value either "epsilon" or "a" must be set.

    epsilon : float, default=0.05
        Precision level (typically 0.05 or 0.1).
        For calculation of the number of samples achieving the certain bound.
        n_samples argument must not be provided and for initialising the
        default value either "delta" or "a" must be set.

    a : int
        Number of isomorphism classes of graphlets.
        If -1 the number is the maximum possible, from a database 1 until 9
        or else predicted through interpolation.
        For calculation of the number of samples achieving the certain bound.
        n_samples argument must not be provided and for initialising the
        default value either "delta" or "epsilon" must be set.

    n_samples : int
        Sets the value of randomly drawn random samples,
        from sizes between 3..k

    Returns
    -------
    kernel : number
        If either "delta", "epsilon", "a" or "n_samples" is given calculates
        the kernel value for the given (or derived) random picked n_samples, by
        randomly sampling from k from 3 to 5.
        Otherwise calculates the kernel value drawing all possible connected
        samples of size k.
        The kernel value.

    """
    Gx = graph(X)
    Gy = graph(Y)

    return graphlet_sampling_matrix({0: Gx}, {0: Gy}, k=k, **kargs)[0, 0]


def graphlet_sampling_matrix(Graphs_x, Graphs_y, k=5, **kargs):
    """Graphlet sampling kernel :cite:`Shervashidze2009EfficientGK`.

    Parameters
    ----------
    Graphs_{x,y} : dict, default_y=None
        Enumerative dictionary of graph type objects with keys from 0 to the
        number of values. If value of Graphs_y is None the kernel matrix is
        computed between all pairs of Graphs_x, where in another case the
        kernel_matrix rows correspond to elements of Graphs_y, and columns
        to the elements of Graphs_x.

    k : int, default=5
        The dimension of the given graphlets.

    delta : float, default=0.05
        Confidence level (typically 0.05 or 0.1).
        For calculation of the number of samples achieving the certain bound.
        n_samples argument must not be provided and for initialising the
        default value either "epsilon" or "a" must be set.

    epsilon : float, default=0.05
        Precision level (typically 0.05 or 0.1).
        For calculation of the number of samples achieving the certain bound.
        n_samples argument must not be provided and for initialising the
        default value either "delta" or "a" must be set.

    a : int
        Number of isomorphism classes of graphlets.
        If -1 the number is the maximum possible, from a database 1 until 9
        or else predicted through interpolation.
        For calculation of the number of samples achieving the certain bound.
        n_samples argument must not be provided and for initialising the
        default value either "delta" or "epsilon" must be set.

    n_samples : int
        Sets the value of randomly drawn random samples,
        from sizes between 3..k.

    Returns
    -------
    kernel_matrix : np.array
        If either "delta", "epsilon", "a" or "n_samples" is given calculates
        the kernel matrix for the given (or derived) random picked n_samples,
        by randomly sampling from k from 3 to 5.
        Otherwise calculates the kernel value drawing all possible
        connected samples of size k.
        The kernel matrix. If the Graphs_y is not None the rows correspond
        to Graphs_y and the cols to the Graphs_x (based on the given order).

    """
    if k > 10:
        warnings.warn('graphlets are too big - computation may be slow')
    elif k < 3:
        raise ValueError('k must be bigger than 3')

    if "n_samples" in kargs:
        n_samples = kargs["n_samples"]
        l = []
        if "delta" in kargs:
            l.append('"delta"')
        if "epsilon" in kargs:
            l.append('"epsilon"')
        if "a" in kargs:
            l.append('"a"')
        if len(l):
            warnings.warn('Number of samples defined as input, \
            ignoring arguments:', ', '.join(l))
        sample_graphlets = sample_graphlets_probabilistic
    elif "delta" in kargs or "epsilon" in kargs or "a" in kargs:
        sample_graphlets = sample_graphlets_probabilistic
        delta = kargs.get("delta", 0.05)
        epsilon = kargs.get("epsilon", 0.05)
        a = kargs.get("a", -1)

        if delta > 1 or delta < 0:
            raise ValueError('delta must be in the range (0,1)')

        if epsilon > 1 or epsilon < 0:
            raise ValueError('epsilon must be in the range (0,1)')

        fallback_map = {1: 1, 2: 2, 3: 4, 4: 8, 5: 19, 6: 53, 7:
                        209, 8: 1253, 9: 13599}

        if type(a) is not int:
            raise ValueError('a must be an integer')
        elif a == 0:
            raise ValueError('a cannot be zero')
        elif a < -1:
            raise ValueError('negative a smaller than -1 have no meaning')

        if(a == -1):
            if(k > 9):
                warnings.warn('warning for such size number of isomorphisms is \
                not known - interpolation on know values will be used')
                # Use interpolations
                isomorphism_prediction = \
                    interp1d(list(fallback_map.keys()),
                             list(fallback_map.values()), kind='cubic')
                a = isomorphism_prediction(k)
            else:
                a = fallback_map[k]
        # Calculate number of samples
        n_samples = math.ceil(2*(a*np.log10(2) +
                              np.log10(1/delta))/(epsilon**2))
    else:
        sample_graphlets = sample_graphlets_all_connected

    nx = len(Graphs_x.keys())

    if Graphs_y is None:
        ng = nx
        Gs = Graphs_x
    else:
        ng = nx + len(Graphs_y.keys())
        Gs = {i: g for (i, g) in
              enumerate(itertools.chain(Graphs_x.values(), Graphs_y.values()))}

    graph_bins = dict()
    local_values = dict()
    for i in range(ng):
        # take n_samples from graphs Gs[i] of sizes between 3..k
        samples = sample_graphlets(
            (Gs[i].get_adjacency_matrix() > 0).astype(int),
            k=k,
            n_samples=n_samples
        )

        for sg in samples:
            # add the graph to an isomorphism class
            if len(graph_bins) == 0:
                graph_bins[0] = sg
                local_values[(i, 0)] = 1
            else:
                newbin = True
                for j in range(len(graph_bins)):
                    if pynauty.isomorphic(graph_bins[j], sg):
                        newbin = False
                        if (i, j) not in local_values:
                            local_values[(i, j)] = 1
                        local_values[(i, j)] += 1
                        break
                if newbin:
                    local_values[(i, len(graph_bins))] = 1
                    graph_bins[len(graph_bins)] = sg

    phi_x = np.zeros((ng, len(graph_bins)))
    for ((i, j), v) in local_values.items():
        phi_x[i, j] = v

    if Graphs_y is None:
        phi_y = phi_x.T
    else:
        phi_y = phi_x[nx:ng, :].T
        phi_x = phi_x[0:nx, :]

    return np.dot(phi_x, phi_y).T


def sample_graphlets_probabilistic(A, k, n_samples):
    """Propabilistical sampling of n_samples of 3..k sized graphs.

    Parameters
    ----------
    A : np.array
        A binary array defining a certain graph.

    k : int
        The maximum dimension of the sampled graphlets.

    n_samples : int
        Sets the value of randomly drawn random samples,
        from sizes between 3..k

    Returns
    -------
    graphlets : generator
        Returns a generator of sampled graphlets (as pynauty graphs),
        from sizes between 3..k.

    """
    s = list(range(A.shape[0]))
    to_edge_dict_binary = lambda x, k: matrix_to_dict(x, '==', 1, False)
    for i in range(n_samples):
        index_rand = random.sample(s, random.randint(3, k))
        Q = A[index_rand, :][:, index_rand]
        yield pynauty.Graph(Q.shape[0], True, to_edge_dict_binary(Q, k))


def sample_graphlets_all_connected(A, k):
    """All the connected graphlets of size k of a given graph.

    Parameters
    ----------
    A : np.array
        A binary array defining a certain graph.

    k : int
        The maximum dimension of the sampled graphlets.

    Returns
    -------
    graphlets : generator
        Returns a generator of sampled graphlets (as pynauty graphs),
        of size k.

    """
    to_edge_dict_binary = lambda x: matrix_to_dict(x, '==', 1, False)
    for i in itertools.permutations(range(A.shape[0]), k):
        Q = A[i, :][:, i]
        if 0 not in np.sum(Q, axis=1):
            yield pynauty.Graph(Q.shape[0], True, to_edge_dict_binary(Q))
