"""The neighborhood hashing kernel as defined in :cite:`Hido2009ALG`."""

import itertools
import os

import numpy as np

from grakel.graph import graph
from grakel.tools import rotl
from grakel.tools import rotr


def neighborhood_hash(X, Y, Lx, Ly, nh_type='simple', R=3, bytes=2):
    """Neighborhood hashing kernel as proposed in :cite:`Hido2009ALG`.

    Parameters
    ----------
    X,Y : *valid-graph-format*
        The pair of graphs on which the kernel is applied.

    L{x,y} : dict
        Corresponding graph labels for vertices.

    R : int, default=3
        The maximum number of neighborhood hash.

    nh_type : str, valid_types={"simple", "count_sensitive"}, default="simple"
        The existing neighborhood hash type as defined in :cite:`Hido2009ALG`.

    bytes : int, default=2
        Byte size of hashes.

    Returns
    -------
    kernel : number
        The kernel value.

    """
    Gx = graph(X, Lx)
    Gy = graph(Y, Ly)
    return float(
        neighborhood_hash_matrix(
            {0: Gx}, {0: Gy}, nh_type=nh_type, R=R, bytes=bytes)[0, 0])


def neighborhood_hash_matrix(
        Graphs_x, Graphs_y=None, nh_type='simple', R=3, bytes=2):
    """`calculate_similarity_matrix` function as defined :cite:`Hido2009ALG`.

    Parameters
    ----------
    Graphs_{x,y} : dict, default_y=None
        Enumerative dictionary of graph type objects with keys from 0 to the
        number of values. If value of Graphs_y is None the kernel matrix is
        computed between all pairs of Graphs_x, where in another case the
        kernel_matrix rows correspond to elements of Graphs_y, and columns
        to the elements of Graphs_x.

    R : int, default=3
        The maximum number of neighborhood hash.

    nh_type : str, valid_values={'simple','count-sensitive'},
        The existing neighborhood hash type as defined in :cite:`Hido2009ALG`.

    bytes : int, default=2
        Byte size of hashes.

    Returns
    -------
    kernel_matrix : np.array
        The kernel matrix. If the Graphs_y is not None the rows correspond to
        Graphs_y and the cols to the Graphs_x (based on the given order).

    """
    if R <= 0:
        raise ValueError('R must be bigger than zero')

    if nh_type == 'simple':
        noc_f = False
        NH = lambda G: neighborhood_hash_simple(G)
    elif nh_type == 'count-sensitive':
        noc_f = True
        NH = lambda G: neighborhood_hash_count_sensitive(G)
    else:
        raise ValueError('unrecognised neighborhood hashing type')

    bytes = int(bytes)
    if bytes <= 0:
        raise ValueError('illegal number of bytes for hashing')

    h_x = len(Graphs_x.keys())
    Gs = dict()

    if Graphs_y is None:
        h_y = h_x
        g_iter = Graphs_x.values()
        ng = h_x
        pairs = [(i, j) for i in range(0, h_x) for j in range(i, h_x)]
        offset = 0
    else:
        h_y = len(Graphs_y.keys())
        g_iter = itertools.chain(Graphs_x.values(), Graphs_y.values())
        ng = h_x+h_y
        pairs = list(itertools.product(range(h_x, ng), range(0, h_x)))
        offset = h_x
    labels_hash_dict, labels_hash_set = dict(), set()
    for (i, g) in enumerate(g_iter):
        g.desired_format('adjacency')
        vertices = list(g.get_vertices())
        labels = hash_labels(
            g.get_labels(), labels_hash_dict, labels_hash_set, bytes)
        Gs[i] = tuple(radix_sort(vertices, labels)) + \
            ({v: g.neighbors(v, purpose="adjacency") for v in vertices},)
        if noc_f:
            noc = dict()
            for k in labels.keys():
                if labels[k] not in noc:
                    noc[labels[k]] = 1
                else:
                    noc[labels[k]] += 1
            Gs[i] += (noc,)

    S = np.zeros(shape=(h_x, h_y))
    for r in range(0, R):
        K = np.eye(h_x, h_y)
        for i in range(0, ng):
            Gs[i] = NH(Gs[i])
        for (i, j) in pairs:
            # targets - y - graph take the row
            K[i-offset, j] = nh_compare_labels(Gs[i], Gs[j])
        S = np.add(K, S)

    kernel_mat = np.divide(S, R)

    if Graphs_y is None:
        kernel_mat = np.triu(kernel_mat) + np.triu(kernel_mat, 1).T

    return kernel_mat


def hash_labels(labels, labels_hash_dict, labels_hash_set, bytes=2):
    """Hashes existing labels to 16-bit integer.

    Hashing without collisions and with consistency in same labels.

    Parameters
    ----------
    labels : dict
        Labels for vertices.

    labels_hash_dict : dict
        A hash table for labels.

    nh_type : str, valid_values={'simple','count-sensitive'},
        The existing neighborhood hash type as defined in :cite:`Hido2009ALG`.

    bytes : int, default=2
        Byte size of hashes.

    Returns
    -------
    new_labels : dict
        The new hashed labels for vertices.

    """
    bytes = int(bytes)
    if bytes <= 0:
        raise ValueError('illegal number of bytes for hashing')

    new_labels = dict()
    for k in labels.keys():
        if labels[k] not in labels_hash_dict:
            f = True
            while f:
                r = int.from_bytes(os.urandom(bytes), 'little')
                f = r in labels_hash_set
            labels_hash_set.add(r)
            labels_hash_dict[labels[k]] = r
            new_labels[k] = r
        else:
            new_labels[k] = labels_hash_dict[labels[k]]
    return new_labels


def radix_sort(vertices, labels):
    """Sorts vertices based on labels.

    Parameters
    ----------
    vertices : listable
        A listable of vertices.

    labels : dict
        Dictionary of labels for vertices.

    Returns
    -------
    vertices_labels : tuple, len=2
        The sorted vertices based on labels and labels for vertices.

    """
    return (sorted(list(vertices), key=lambda x: labels[x]), labels)


def neighborhood_hash_simple(G):
    """(simple) neighborhood hashing as defined in :cite:`Hido2009ALG`.

    Parameters
    ----------
    G : tuple
        A tuple of three elements consisting of vertices sorted by labels,
        vertex label dictionary and edge dictionary.

    Returns
    -------
    vertices_labels_edges : tuple
        A tuple of vertices, new_labels-dictionary and edges.

    """
    vertices, labels, edges = G
    new_labels = dict()
    for u in vertices:
        label = ROT(labels[u], 1)
        for n in edges[u]:
            label ^= labels[n]
        new_labels[u] = label
    return (vertices, new_labels, edges)


def neighborhood_hash_count_sensitive(G):
    """Count sensitive neighborhood hash as defined in :cite:`Hido2009ALG`.

    Parameters
    ----------
    G : tuple, len=3
       A tuple three elements consisting of vertices sorted by labels,
       vertex label dict, edge dict and number of occurencies dict for labels

    Returns
    -------
    vertices_labels_edges_noc : tuple
        A tuple of 4 elements consisting of vertices sorted by labels,
        vertex label dict, edge dict and number of occurencies dict.

    """
    vertices, labels, edges, noc = G
    new_labels = dict()
    new_noc = dict()
    for u in vertices:
        label = ROT(labels[u], 1)
        for n in edges[u]:
            o = noc[labels[n]]
            label ^= ROT(labels[n] ^ o, o)
        if label not in new_noc:
            new_noc[label] = 1
        else:
            new_noc[label] += 1
        new_labels[u] = label
    return (vertices, new_labels, edges, new_noc)


def ROT(l, n):
    """`rot` operation for binary numbers.

    Parameters
    ----------
    n : int
        An integer.

    l : int
        A valid number.

    Returns
    -------
    rot : int
        The result of a rot operation.

    """
    if n < 0:
        return rotr(l, -n)
    elif n > 0:
        return rotl(l, n)
    else:
        return l


def nh_compare_labels(Gx, Gy):
    """Compare labels function as defined in :cite:`Hido2009ALG`.

    Parameters
    ----------
    G_{x,y} : tuple, len=2
        Graph tuples of two elements, consisting of vertices sorted by
        (labels, vertices) and edge-labels dict.

    Returns
    -------
    kernel : Number
        The kernel value.

    """
    # get vertices
    vx, vy = iter(Gx[0]), iter(Gy[0])

    # get size of vertices
    nv_x, nv_y = len(Gx[0]), len(Gy[0])

    # get labels for vertices
    Lx, Ly = Gx[1], Gy[1]

    c = 0
    ui, uj = next(vx, None), next(vy, None)
    while (ui is not None) and (uj is not None):
        if Lx[ui] == Ly[uj]:
            c += 1
            ui = next(vx, None)
            uj = next(vy, None)
        elif Lx[ui] < Ly[uj]:
            ui = next(vx, None)
        else:
            uj = next(vy, None)

    return c/float(nv_x+nv_y-c)
