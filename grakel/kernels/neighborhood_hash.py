"""The neighborhood hashing kernel as defined in :cite:`Hido2009ALG`."""
import collections
import os
import warnings

from grakel.graph import graph
from grakel.tools import rotl
from grakel.tools import rotr
from grakel.kernels import kernel


class neighborhood_hash(kernel):
    """Neighborhood hashing kernel as proposed in :cite:`Hido2009ALG`.

    Parameters
    ----------
    R : int, default=3
        The maximum number of neighborhood hash.

    nh_type : str, valid_types={"simple", "count_sensitive"}, default="simple"
        The existing neighborhood hash type as defined in :cite:`Hido2009ALG`.

    bytes : int, default=2
        Byte size of hashes.

    Attributes
    -------
    _R : number
        The maximum number of neighborhood hash.

    _NH : function
        The neighborhood hashing function.

    _noc_f : bool
        A flag concerning the number of occurencies metric.

    _bytes : int
        Defines the byte size of hashes.

    """

    def __init__(self, **kargs):
        """Initialise a `pyramid_match` kernel."""
        self._valid_parameters |= {"R", "nh_type", "byte"}
        super(neighborhood_hash, self).__init__(**kargs)

        self._R = kargs.get("R", 3)

        if self._R <= 0:
            raise ValueError('R must be bigger than zero')

        nh_type = kargs.get("nh_type", "simple")
        if nh_type == 'simple':
            self._noc_f = False
            self._NH = lambda G: neighborhood_hash_simple(G)
        elif nh_type == 'count-sensitive':
            self._noc_f = True
            self._NH = lambda G: neighborhood_hash_count_sensitive(G)
        else:
            raise ValueError('unrecognised neighborhood hashing type')

        self._bytes = int(kargs.get("bytes", 2))
        if self._bytes <= 0:
            raise ValueError('illegal number of bytes for hashing')
        self._default_nonce_hash = int((self._bytes*8+1) * '1', 2)

    def parse_input(self, X):
        """Parse and create features for shortest_path kernel.

        Parameters
        ----------
        X : object
            For the input to pass the test, we must have:
            Each element must be an iterable with at most three features and at
            least one. The first that is obligatory is a valid graph structure
            (adjacency matrix or edge_dictionary) while the second is
            node_labels and the third edge_labels (that fitting the given graph
            format). If None the kernel matrix is calculated upon fit data.
            The test samples.

        Returns
        -------
        H : list
            A list of lists of Histograms for all levels for each graph.

        """
        if not isinstance(X, collections.Iterable):
            raise ValueError('input must be an iterable\n')
        else:
            i = 0
            gs = list()
            if self._method_calling in [1, 2]:
                self._labels_hash_dict = dict()
                self._labels_hash_set = set()
            for x in iter(X):
                if len(x) == 0:
                    warnings.warn('Ignoring empty element on index: '+str(i))
                if len(x) == 1:
                    if type(x) is graph:
                        v = list(x.get_vertices(purpose="any"))
                        L = x.get_labels(purpose="any")
                    else:
                        warnings.warn(
                            'Ignoring empty element on index: '
                            + str(i) + '\nLabels must be provided.')
                    i += 1
                elif len(x) in [2, 3]:
                    x = graph(x[0], x[1], {}, self._graph_format)
                    v = list(x.get_vertices(purpose="any"))
                    L = x.get_labels(purpose="any")
                    i += 1
                else:
                    raise ValueError('each element of X must have at least' +
                                     ' one and at most 3 elements\n')

                if self._method_calling in [1, 2]:
                    labels = self.hash_labels(L)

                g = tuple(radix_sort(v, L)) + \
                    ({n: x.neighbors(n, purpose="any") for n in v},)

                if self._noc_f:
                    noc = dict()
                    for k in labels.keys():
                        if labels[k] not in noc:
                            noc[labels[k]] = 1
                        else:
                            noc[labels[k]] += 1
                    g += (noc,)

                gr = dict()
                gr[0] = self._NH(g)
                for r in range(1, self._R):
                    gr[r] = self._NH(gr[r-1])
                gs.append(gr)

            if i == 0:
                raise ValueError('parsed input is empty')

            return gs

    def pairwise_operation(self, x, y):
        """Calculate a pairwise kernel between two elements.

        Parameters
        ----------
        x, y : dict
            Dict of len=2, tuples, consisting of vertices sorted by
            (labels, vertices) and edge-labels dict, for each radius
            from 0 .. R-1.

        Returns
        -------
        kernel : number
            The kernel value.

        """
        return sum(nh_compare_labels(x[r], y[r]) for r in range(self._R))

    def hash_labels(self, labels):
        """Hashes existing labels to 16-bit integer.

        Hashing without collisions and with consistency in same labels.

        Parameters
        ----------
        labels : dict
            Labels for vertices.

        labels_hash_dict : dict
            A hash table for labels.

        nh_type : str, valid_values={'simple','count-sensitive'},
            The existing neighborhood hash type as defined in
            :cite:`Hido2009ALG`.

        bytes : int, default=2
            Byte size of hashes.

        Returns
        -------
        new_labels : dict
            The new hashed labels for vertices.

        """
        new_labels = dict()
        for k in labels.keys():
            if labels[k] not in self._labels_hash_dict:
                if self._method_calling not in [1, 2]:
                    f = True
                    while f:
                        r = int.from_bytes(os.urandom(self._bytes), 'little')
                        f = r in self._labels_hash_set
                    self._labels_hash_set.add(r)
                    self._labels_hash_dict[labels[k]] = r
                    new_labels[k] = r
                else:
                    new_labels[k] = self._default_nonce_hash
            else:
                new_labels[k] = self._labels_hash_dict[labels[k]]
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
