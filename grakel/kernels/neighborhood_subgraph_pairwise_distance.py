"""Neighborhood subgraph pairwise distance kernel :cite:`Costa2010FastNS`."""
import collections
import warnings

import numpy as np

from sklearn.exceptions import NotFittedError
from sklearn.utils.validation import check_is_fitted

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

    _fit_keys : dict
        A dictionary with keys from 0 to _d+1, constructed upon fit
        holding an enumeration of all the found (in the fit dataset)
        tuples of two hashes and a radius in this certain level.

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
        M : dict
            A dictionary with keys all the distances from 0 to self._d
            and values the the np.arrays with rows corresponding to the
            non-null input graphs and columns to the enumerations of tuples
            consisting of pairs of hash values and radius, from all the given
            graphs of the input (plus the fitted one's on transform).

        Z : np.array, shape=(self._d+1, n_inputs)
            A numpy array consisting of the number of distance pairs for each
            level. If no pairs appear, a value of float("Inf") is provided.

        """
        if not isinstance(X, collections.Iterable):
            raise ValueError('input must be an iterable\n')
        else:
            i = 0
            data = {j: list() for j in range(self._d+1)}
            all_keys = {j: dict() for j in range(self._d+1)}
            for (idx, x) in enumerate(iter(X)):
                is_iter = False
                if isinstance(x, collections.Iterable):
                    is_iter, x = True, list(x)
                if is_iter and len(x) in [0, 3]:
                    if len(x) == 0:
                        warnings.warn('Ignoring empty element' +
                                      ' on index: '+str(idx))
                        continue
                    else:
                        g = Graph(x[0], x[1], x[2])
                        g.change_format("adjacency")
                elif type(x) is Graph:
                    g = Graph(x.get_adjacency_matrix(),
                              x.get_labels(purpose="adjacency",
                                           label_type="vertex"),
                              x.get_labels(purpose="adjacency",
                                           label_type="edge"))

                else:
                    raise ValueError('each element of X must have either ' +
                                     'a graph with labels for node and edge ' +
                                     'or 3 elements consisting of a graph ' +
                                     'type object, labels for vertices and ' +
                                     'labels for edges.')
                g.change_format(self._graph_format)
                vertices = set(g.get_vertices(purpose=self._graph_format))

                ed = g.get_edge_dictionary()

                edges = {(j, k) for j in ed.keys() for k in ed[j].keys()}
                Lv = g.get_labels(purpose=self._graph_format)
                Le = g.get_labels(purpose=self._graph_format,
                                  label_type="edge")
                N, D, D_pair = g.produce_neighborhoods(
                    self._r, purpose="dictionary",
                    with_distances=True, d=self._d)

                H = self._hash_neighborhoods(vertices, edges, Lv,
                                             Le, N, D_pair)

                if self._method_calling == 1:
                    for d in range(self._d+1):
                        keys = all_keys[d]
                        Q = dict()
                        if d in D:
                            for (A, B) in D[d]:
                                for r in range(self._r+1):
                                    key = (H[r][A], H[r][B], r)
                                    idx = keys.get(key, None)
                                    if idx is None:
                                        idx = len(keys)
                                        keys[key] = idx
                                    Q[idx] = Q.get(idx, 0) + 1
                        data[d].append(Q)

                elif self._method_calling == 3:
                    for d in range(self._d+1):
                        keys = all_keys[d]
                        fit_keys = self._fit_keys[d]
                        Q = dict()
                        if d in D:
                            for (A, B) in D[d]:
                                for r in range(self._r+1):
                                    key = (H[r][A], H[r][B], r)
                                    idx = fit_keys.get(key, None)
                                    if idx is None:
                                        idx = keys.get(key, None)
                                        if idx is None:
                                            idx = len(keys) + len(fit_keys)
                                            keys[key] = idx
                                    Q[idx] = Q.get(idx, 0) + 1
                        data[d].append(Q)
                i += 1
            if i == 0:
                raise ValueError('parsed input is empty')

            if self._method_calling == 1:
                M = {d: np.zeros(shape=(i, len(all_keys[d])))
                     for d in range(self._d+1)}
                Z = np.zeros(shape=(self._d+1, i))
                for d in range(self._d+1):
                    for (no, dt) in enumerate(data[d]):
                        for (idx, val) in dt.items():
                            M[d][no, idx] = val
                            Z[d, no] += val
                self._fit_keys = all_keys
            elif self._method_calling == 3:
                M = {d: np.zeros(shape=(i, len(all_keys[d]) +
                                 len(self._fit_keys[d])))
                     for d in range(self._d+1)}
                Z = np.zeros(shape=(self._d+1, i))
                for d in range(self._d+1):
                    for (no, dt) in enumerate(data[d]):
                        for (idx, val) in dt.items():
                            M[d][no, idx] = val
                            Z[d, no] += val
            Z[(Z == 0)] = float("Inf")
            return (M, Z)

    def _calculate_kernel_matrix(self, Y=None):
        """Calculate the kernel matrix given a target_graph and a kernel.

        Each a matrix is calculated between all elements of Y on the rows and
        all elements of X on the columns.

        Parameters
        ----------
        Y : list, default=None
            A list of graph type objects. If None kernel is calculated between
            X and itself.

        Returns
        -------
        K : numpy array, shape = [n_targets, n_inputs]
            The kernel matrix: a calculation between all pairs of graphs
            between targets and inputs. If Y is None targets and inputs
            are the taken from self.X. Otherwise Y corresponds to targets
            and self.X to inputs.

        """
        M, Z = self.X
        if Y is None:
            S = np.zeros(shape=(Z.shape[1], Z.shape[1]))
            for d in range(self._d+1):
                S += (np.dot(M[d], M[d].T) / np.outer(Z[d, :], Z[d, :]))

        else:
            Mp, Zp = Y
            S = np.zeros(shape=(Zp.shape[1], Z.shape[1]))
            for d in range(self._d+1):
                S += (np.dot(Mp[d][:, :M[d].shape[1]], M[d].T) /
                      np.outer(Zp[d, :], Z[d, :]))

        return S

    def diagonal(self):
        """Calculate the kernel matrix diagonal of the fitted data.

        Parameters
        ----------
        None.

        Returns
        -------
        X_diag : np.array
            The diagonal of the kernel matrix, of the fitted. This consists
            of each element calculated with itself.


        """
        # Check is fit had been called
        check_is_fitted(self, ['X', '_Y'])
        try:
            check_is_fitted(self, ['_X_diag'])
        except NotFittedError:
            # Calculate diagonal of X
            Z, M = self.X
            self._X_diag = np.empty(shape=(Z.shape[1], 1))
            for d in range(self._d+1):
                self._X_diag += (M[d] * M[d].T).sum(-1) / np.square(Z[d, :])

        Z, M = self._Y
        Y_diag = np.empty(shape=(Z.shape[1], 1))
        for d in range(self._d+1):
            Y_diag += (M[d] * M[d].T).sum(-1) / np.square(Z[d, :])

        return self._X_diag, Y_diag

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
                sub_vertices = sorted(N[radius][v])
                re = {(i, j) for (i, j) in re
                      if i in sub_vertices and j in sub_vertices}
                lv = {v: lv[v] for v in sub_vertices}
                le = {e: le[e] for e in edges}
                H[radius][v] = hash_graph(D_pair, sub_vertices, re, lv, le)
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
