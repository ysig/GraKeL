"""Neighborhood subgraph pairwise distance kernel :cite:`Costa2010FastNS`."""
import itertools
import collections
import warnings

from grakel.kernels import kernel
from grakel.graph import Graph

from grakel.kernels._c_functions import APHash


class neighborhood_subgraph_pairwise_distance(kernel):
    """The Neighborhood subgraph pairwise distance kernel.

    See :cite:`Costa2010FastNS`.

    Parameters
    ----------
    r : int, default=3
        The maximum considered radius between vertices.

    d : int, default=4
        Neighborhood depth.

    Attributes
    ----------
    _r : int
        The maximum considered radius between vertices.

    _d : int
        Neighborhood depth.


    """

    _graph_format = "dictionary"

    def __init__(self, **kargs):
        """Initialize an NSPD kernel."""
        # setup valid parameters and initialise from parent
        self._valid_parameters |= {"r", "d"}
        super(neighborhood_subgraph_pairwise_distance, self).__init__(**kargs)

        self._r = kargs.get("r", 3)
        self._d = kargs.get("d", 4)

        if self._r < 0:
            raise ValueError('r must be a positive integer')

        if self._d < 0:
            raise ValueError('d must be a positive integer')

    def parse_input(self, X):
        """Parse and create features for graphlet_sampling kernel.

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
        out : list
            A list of tuples consisting of a 2-level dictionary of the hashed
            neighborhoods from radious, vertex to the hashed values and a
            dictionary, where for each level appears a set of tuples of
            nodes connected in that level.


        """
        if not isinstance(X, collections.Iterable):
            raise ValueError('input must be an iterable\n')
        else:
            i = 0
            out = list()
            for (idx, x) in enumerate(iter(X)):
                if type(x) is Graph:
                    g = Graph(x.get_adjacency_matrix(),
                              x.get_labels(purpose="adjacency",
                                           label_type="vertex"),
                              x.get_labels(purpose="adjacency",
                                           label_type="edge"))
                elif isinstance(x, collections.Iterable) and \
                        len(x) in [0, 3]:
                    if len(x) == 0:
                        warnings.warn('Ignoring empty element' +
                                      ' on index: '+str(idx))
                        continue
                    else:
                        g = Graph(x[0], x[1], x[2])
                        g.change_format("adjacency")
                else:
                    raise ValueError('each element of X must have either ' +
                                     'a graph with labels for node and edge ' +
                                     'or 3 elements consisting of a graph ' +
                                     'type object, labels for vertices and ' +
                                     'labels for edges.')
                g.change_format(self._graph_format)
                vertices = set(g.get_vertices(purpose=self._graph_format))

                ed = g.get_edge_dictionary()

                edges = {(i, j) for i in ed.keys() for j in ed[i].keys()}
                Lv = g.get_labels(purpose=self._graph_format)
                Le = g.get_labels(purpose=self._graph_format,
                                  label_type="edge")
                N, D, D_pair = g.produce_neighborhoods(
                    self._r, purpose="dictionary",
                    with_distances=True, d=self._d)

                H = self._hash_neighborhoods(vertices, edges, Lv,
                                             Le, N, D_pair)
                out.append((H, D))
                i += 1
            if i == 0:
                raise ValueError('parsed input is empty')
            return out

    def pairwise_operation(self, x, y):
        """Neighborhood subgraph pairwise distance kernel.

        See :cite:`Costa2010FastNS`.

        Parameters
        ----------
        x, y : tuple
            A tuple same as each element of the output list of parse input.
            It consists of a 2-level dictionary of the hashed neighborhoods
            from radious, vertex to the hashed values and a dictionary,
            where for each level appears a set of tuples of nodes
            connected in that level.

        Returns
        -------
        kernel : number
            The kernel value.

        """
        Hx, Dx = x
        Hy, Dy = y

        kernel = 0

        for distance in range(self._d+1):
            if distance in Dx and distance in Dy:
                pairs = itertools.product(Dx[distance], Dy[distance])
                krd = 0
                npairs = 0
                for ((A, B), (Ap, Bp)) in pairs:
                    npairs += 1
                    for radius in range(self._r+1):
                        if Hx[radius][A] == Hy[radius][Ap]:
                            krd += int(Hx[radius][B] == Hy[radius][Bp])
                if npairs > 0:
                    # normalization by the number of pairs
                    kernel += float(krd)/npairs

        return kernel

    def _hash_neighborhoods(self, vertices, edges, Lv, Le, N, D_pair):
        """Hash all neighborhoods and all root nodes.

        Parameters
        ----------
        vertices : set
            The graph vertices.

        edges : set
            The set of edges

        N : dict
            Neighborhoods that map levels (int) to dictionaries of root node
            symbols (keys) to list of vertex symbols, which correspond to the
            neighbors, that belong to this neighborhood.

        D_pairs : dict
            A dictionary that maps edges (tuple pairs of vertex symbols) to
            element distances (int - as produced from a BFS traversal).

        Returns
        -------
        H : dict
            The hashed neighborhoods as a 2-level dict from radious,
            vertex to the hashed values.

        """
        H = {ra: dict() for ra in range(self._r+1)}
        sel = sorted(list(edges))
        for v in vertices:
            re, lv, le = sel, Lv, Le
            for radius in range(self._r, -1, -1):
                vertices = sorted(N[radius][v])
                re = {(i, j) for (i, j) in re
                      if i in vertices and j in vertices}
                lv = {v: lv[v] for v in vertices}
                le = {e: le[e] for e in edges}
                H[radius][v] = hash_graph(D_pair, vertices, re, lv, le)
        return H


def hash_graph(D, vertices, edges, glv, gle):
    """Make labels for hashing according to the proposed method.

    Produces the graph hash needed for fast comparison.

    Parameters
    ----------
    D_pairs : dict
        A dictionary that maps edges (tuple pairs of vertex symbols) to
        element distances (int - as produced from a BFS traversal).

    vertices : set
        A set of vertices.

    edges : set
        A set of edges.

    glv : dict
        Labels for vertices of the graph.

    gle : dict
        Labels for edges of the graph.

    Returns
    -------
    hash : int.
        The hash value for the given graph.

    """
    encoding = ""

    # Make labels for vertices
    Lv = dict()
    for i in vertices:
        label = sorted([(str(D[(i, j)]) + ',' + str(glv[j])) for j in vertices
                        if j != i and (i, j) in D])
        encoding += str(label)+"."
        Lv[i] = label
    encoding = encoding[:-1]+":"

    # Expand to labels for edges
    Le = dict()
    for (i, j) in edges:
        Le[(i, j)] = str(Lv[i]) + ',' + str(Lv[j]) + ',' + str(gle[(i, j)])
        encoding += str(Le[(i, j)])+"_"

    # Arash Partov hashing, as in the original
    # implementation of NSPK.
    return APHash(encoding)
