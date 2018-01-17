"""The Jensen-Shannon Representation Alignment kernel :cite:`Bai2015AGK`."""
import math

import numpy as np

from grakel.graph import graph


def jsm(X, Y, M=3, h=3):
    """Jensen-Shannon Representation Alignment kernel :cite:`Bai2015AGK`.

    Parameters
    ----------
    X,Y : *valid-graph-format*
        The pair of graphs on which the kernel is applied.

    M : int, default=3
        The maximum number of JS-representations considered.

    h : int, default=3
        The number of layers of DB representations.

    Returns
    -------
    kernel : number
        The kernel value.

    """
    Gx = graph(X, {})
    Gy = graph(Y, {})
    return jsm_pair(Gx, Gy, M=M, h=h)


def jsm_pair(Gx, Gy, M=3, h=3):
    """Jensen-Shannon Representation Alignment, matrices :cite:`Bai2015AGK`.

    Parameters
    ----------
    G{x,y} : graph
        The pair of graphs on which the kernel is applied.

    M : int, default=3
        The maximum number of JS-representations considered.

    h : int, default=3
        The number of layers of DB representations.

    Returns
    -------
    kernel : number
        The kernel value.

    """
    if h <= 0:
        raise ValueError("h must be a positive integer")

    if M <= 0:
        raise ValueError("M must be a positive integer")

    Nx, Dx, _ = Gx.produce_neighborhoods(h, "adjacency", True, max(h, M))
    Ny, Dy, _ = Gy.produce_neighborhoods(h, "adjacency", True, max(h, M))
    nx = Gx.nv()
    ny = Gy.nv()

    Ex = get_level_edges(Dx, Nx, h, nx)
    Ey = get_level_edges(Dy, Ny, h, ny)

    kernel = 0
    for m in range(M):
        m_spherex = calculate_m_sphere(Dx[m], nx)
        m_spherey = calculate_m_sphere(Dy[m], ny)
        JG_mh_x = JS_representation(Nx, Ex, m_spherex, h, nx)
        JG_mh_y = JS_representation(Ny, Ey, m_spherey, h, ny)
        R = np.empty(shape=(nx, ny))
        for i in range(nx):
            for j in range(ny):
                R[i, j] = np.linalg.norm(JG_mh_x[i, :] - JG_mh_y[j, :], ord=2)

        min_row = np.min(R, axis=1)
        min_col = np.min(R, axis=0)
        k = 0
        for i in range(nx):
            for j in range(ny):
                if (R[i, j] == min_row[i] and R[i, j] == min_col[j]
                        and len(m_spherex[i]) == len(m_spherey[j])):
                    k += 1
        kernel += k

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
        hs = sum(-p*math.log(p) for p in P.values())
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
    NP = sum(P.values())
    for k in P.keys():
        P[k] /= NP
    return P
