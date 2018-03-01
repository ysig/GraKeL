"""The Jensen-Shannon Representation Alignment kernel :cite:`Bai2015AGK`."""
import collections
import math
import warnings

import numpy as np

from grakel.graph import Graph
from grakel.kernels import kernel

# Python 2/3 cross-compatibility import
from six import itervalues


class jsm(kernel):
    """Jensen-Shannon Representation Alignment kernel :cite:`Bai2015AGK`.

    Parameters
    ----------
    X,Y : *valid-graph-format*
        The pair of graphs on which the kernel is applied.

    M : int, default=3
        The maximum number of JS-representations considered.

    h : int, default=3
        The number of layers of DB representations.

    Attributes
    ----------
    _M : int
        The maximum size of the M-spheres.

    _h : int
        The maximum depth of.

    _max_mh : int
        A nonce variable holding the max of m, h.

    """

    _graph_format = "adjacency"

    def __init__(self, **kargs):
        """Initialise a lovasz_theta kernel."""
        # setup valid parameters and initialise from parent
        self._valid_parameters |= {"M", "h"}
        super(jsm, self).__init__(**kargs)

        self._M = kargs.get("M", 3)
        if type(self._M) is not int or self._M <= 0:
            raise ValueError('M must be an integer bigger than zero')

        self._h = kargs.get("h", 3)
        if type(self._h) is not int or self._h <= 0:
            raise ValueError('M must be an integer bigger than zero')

        self._max_mh = max(self._M, self._h)

    def parse_input(self, X):
        """Parse and create features for graphlet_sampling kernel.

        Parameters
        ----------
        X : iterable
            For the input to pass the test, we must have:
            Each element must be an iterable with at most three features and at
            least one. The first that is obligatory is a valid graph structure
            (adjacency matrix or edge_dictionary) while the second is
            node_labels and the third edge_labels (that correspond to the given
            graph format). A valid input also consists of graph type objects.

        Returns
        -------
        out : list
            The lovasz metrics for the given input.

        """
        if not isinstance(X, collections.Iterable):
            raise ValueError('input must be an iterable\n')
        else:
            i = 0
            out = list()
            for (idx, x) in enumerate(iter(X)):
                is_iter = isinstance(x, collections.Iterable)
                if is_iter:
                    x = list(x)
                if is_iter and len(x) in [0, 1, 2, 3]:
                    if len(x) == 0:
                        warnings.warn('Ignoring empty element ' +
                                      'on index: '+str(idx))
                        continue
                    else:
                        x = Graph(x[0], {}, {}, self._graph_format)
                elif type(x) is Graph:
                    x = Graph(x.get_adjacency_matrix(),
                              {}, {}, self._graph_format)
                else:
                    raise ValueError('each element of X must be either a ' +
                                     'graph or an iterable with at least 1 ' +
                                     'and at most 3 elements\n')
                N, D, _ = x.produce_neighborhoods(
                    self._h, "adjacency", True, self._max_mh)
                n = x.nv()
                E = get_level_edges(D, N, self._h, n)

                m_sphere_len = dict()
                JG = dict()
                for m in range(self._M):
                    m_sphere = calculate_m_sphere(D[m], n)
                    JG[m] = JS_representation(N, E, m_sphere, self._h, n)
                    m_sphere_len[m] = np.empty(shape=(n,))
                    for i in range(n):
                        m_sphere_len[m][i] = len(m_sphere[i])

                out.append((n, JG, m_sphere_len))
                i += 1

            if i == 0:
                raise ValueError('parsed input is empty')

            return out

    def pairwise_operation(self, x, y):
        """Lovasz theta kernel as proposed in :cite:`Johansson2015LearningWS`.

        Parameters
        ----------
        x, y : dict
            Subgraph samples metric dictionaries for all levels.

        Returns
        -------
        kernel : number
            The kernel value.

        """
        nx, JGx, m_sphere_lenx = x
        ny, JGy, m_sphere_leny = y

        kernel = 0
        for m in range(self._M):
            R = np.empty(shape=(nx, ny))
            for i in range(nx):
                for j in range(ny):
                    R[i, j] = np.linalg.norm(JGx[m][i, :] - JGy[m][j, :],
                                             ord=2)

            min_col, min_row = np.min(R, axis=0), np.min(R, axis=1)

            kernel += sum(np.sum(np.logical_and(
                np.logical_and(R[i, :] == min_row[i],
                               R[i, :] == min_col),
                m_sphere_leny[m] == m_sphere_lenx[m][i]))
                for i in range(nx))

        return kernel


def get_level_edges(D, N, h, nv):
    """Derive the graph edges for all the given K_expansions.

    Parameters
    ----------
    D : dict
        Distance level dictionary as produced by produce_neighboorhoods
        method of graph type objects.

    N : dict
        The neighborhood level dictionary as produced by produce_neighborhoods
        method of graph type objects.

    h : int
        The number of layers of DB representations.
        Also a valid number of levels for both dictionaries D, N.

    nv : int
        The number of vertices.

    Returns
    -------
    E : dict
        An edge dictionary for edges, with the same structure as
        the N dictionary, except having sets of touples as values.

    """
    E = {i: {j: set() for j in range(nv)} for i in range(h)}
    for level in range(h):
        for node in range(nv):
            if level > 0:
                E[level][node] = E[level-1][node].copy()
            for (i, j) in D[level]:
                if N[level][node]:
                    E[level][node].add((i, j))
    return E


def calculate_m_sphere(D, nv):
    """Calculate the m sphere.

    Paramaters
    ----------
    D : set
        Distance level dictionary at level m as derived by the produced
        from produce_neighboorhoods method of graph type objects.

    nv : int
        The number of vertices.

    Returns
    -------
    X : dict
        Dictionary of sets of nodes for each nodes corresponding
        to the all the nodes discovered by each node at level m.

    """
    X = {j: set() for j in range(nv)}
    for (v, u) in D:
        X[v].add(u)

    return X


def JS_representation(K_exp_nodes, K_exp_edges, m_sphere, h, nv):
    """m-layer Jensen Shannon representation p. 2, 3 of :cite:`Bai2015AGK`.

    Parameters
    ----------
    K_exp_nodes : dict
        The set of nodes on a given K expansion.

    K_exp_edges : dict
        The set of edges on a given K expansion.

    m_sphere ; dict
        The nodes of the current m-sphere.

    h : int
        The number of layers of DB representations.
        Also a valid number of levels for both dictionaries D, N.

    h : int
        The number of layers of DB representations.
        Also a valid number of levels for both dictionaries D, N.

    nv : int
        The number of vertices.

    Returns
    -------
    J : np.array, shape=(nv,h)
        The Jensen Shanon matrix as produced for all the vertices
        and at height h. See eq. (7) p. 3 of :cite:`Bai2015AGK`.

    """
    J = np.empty(shape=(nv, h))

    for j in range(nv):
        for k in range(h):
            G = [(K_exp_edges[k][j], len(K_exp_nodes[k][j]))]
            for q in m_sphere[j]:
                G.append((K_exp_edges[k][q], len(K_exp_nodes[k][j])))
            J[j, k] = entropy_calculation(G)

    return J


def entropy_calculation(m_neighbor_graph_set):
    r"""Entropy calculation for a single vertex.

    See eq. (1,2,3) for eq. (7) and eq. (5,6) of :cite:`Bai2015AGK`.

    Parameters
    ----------
    m_neighbor_graph_set : list
        A list of graphs touples corresponding to k_expansions of the
        m_neighbor edges and the number of vertices on this expansion,
        as :math:`\mathbf{G}^{K}_{\hat{N}^{m}_{v}}` explained in the firt
        line of p. 3 for eq. (6) of :cite:`Bai2015AGK`.
        The root node should correspond to the first element of the list.

    Returns
    -------
    JS : float
        The m-layer JS representation for a certain vertex at a certain height.

    """
    hdu = 0
    shs = 0
    snv = 0
    for (i, (edges, nv)) in enumerate(m_neighbor_graph_set):
        P = calculate_P(edges)
        hs = sum(-p*math.log(p) for p in itervalues(P))
        if i == 0:
            root_entropy = hs
        shs += hs
        hdu += nv*hs
        snv += nv
    return root_entropy + ((hdu/snv) - (shs/len(m_neighbor_graph_set)))


def calculate_P(edges):
    """Probability calculation of steady state random walks for each vertex.

    As defined after eq. (2) on p. 2 of :cite:`Bai2015AGK`.

    Parameters
    ----------
    edges : list
        A list of tuples corresponding to edges od a given graph.

    Returns
    -------
    P : dict
        For dictionary between all nodes (that have outgoing edges)
        and their probabilities.

    """
    P = dict()
    for (i, j) in edges:
        if i not in P:
            P[i] = 0
        P[i] += 1
    NP = sum(itervalues(P.values))
    for k in P.keys():
        P[k] /= NP
    return P
