"""Neighborhood subgraph pairwise distance kernel :cite:`costa2010fast`."""
# Author: Ioannis Siglidis <y.siglidis@gmail.com>
# License: BSD 3 clause
import collections
import warnings

import numpy as np

from scipy.sparse import csr_matrix

from sklearn.exceptions import NotFittedError
from sklearn.utils.validation import check_is_fitted

from grakel.kernels import Kernel
from grakel.graph import Graph

from grakel.kernels._c_functions import APHash

# Python 2/3 cross-compatibility import
from six import iteritems
from six.moves import filterfalse
from builtins import range


class NeighborhoodSubgraphPairwiseDistance(Kernel):
    """The Neighborhood subgraph pairwise distance kernel.

    See :cite:`costa2010fast`.

    Parameters
    ----------
    r : int, default=3
        The maximum considered radius between vertices.

    d : int, default=4
        Neighborhood depth.

    Attributes
    ----------
    _ngx : int
        The number of graphs upon fit.

    _ngy : int
        The number of graphs upon transform.

    _fit_keys : dict
        A dictionary with keys from `0` to `_d+1`, constructed upon fit
        holding an enumeration of all the found (in the fit dataset)
        tuples of two hashes and a radius in this certain level.

    _X_level_norm_factor : dict
        A dictionary with keys from `0` to `_d+1`, that holds the self
        calculated kernel `[krg(X_i, X_i) for i=1:ngraphs_X]` for all levels.

    """

    _graph_format = "dictionary"

    def __init__(self, n_jobs=None, normalize=False, verbose=False, r=3, d=4):
        """Initialize an NSPD kernel."""
        # setup valid parameters and initialise from parent
        super(NeighborhoodSubgraphPairwiseDistance, self).__init__(
            n_jobs=n_jobs,
            normalize=normalize,
            verbose=verbose)

        self.r = r
        self.d = d
        self._initialized.update({"r": False, "d": False})

    def initialize(self):
        """Initialize all transformer arguments, needing initialization."""
        if not self._initialized["n_jobs"]:
            if self.n_jobs is not None:
                warnings.warn('no implemented parallelization for NeighborhoodSubgraphPairwiseDistance')
            self._initialized["n_jobs"] = True

        if not self._initialized["r"]:
            if type(self.r) is not int or self.r < 0:
                raise ValueError('r must be a positive integer')
            self._initialized["r"] = True

        if not self._initialized["d"]:
            if type(self.d) is not int or self.d < 0:
                raise ValueError('d must be a positive integer')
            self._initialized["d"] = True

    def parse_input(self, X):
        """Parse and create features for the NSPD kernel.

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
        M : dict
            A dictionary with keys all the distances from 0 to self.d
            and values the the np.arrays with rows corresponding to the
            non-null input graphs and columns to the enumerations of tuples
            consisting of pairs of hash values and radius, from all the given
            graphs of the input (plus the fitted one's on transform).

        """
        if not isinstance(X, collections.Iterable):
            raise TypeError('input must be an iterable\n')
        else:
            # Hold the number of graphs
            ng = 0

            # Holds all the data for combinations of r, d
            data = collections.defaultdict(dict)

            # Index all keys for combinations of r, d
            all_keys = collections.defaultdict(dict)
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
                              x.get_labels(purpose="adjacency", label_type="vertex"),
                              x.get_labels(purpose="adjacency", label_type="edge"))
                else:
                    raise TypeError('each element of X must have either ' +
                                    'a graph with labels for node and edge ' +
                                    'or 3 elements consisting of a graph ' +
                                    'type object, labels for vertices and ' +
                                    'labels for edges.')

                # Bring to the desired format
                g.change_format(self._graph_format)

                # Take the vertices
                vertices = set(g.get_vertices(purpose=self._graph_format))

                # Extract the dicitionary
                ed = g.get_edge_dictionary()

                # Convert edges to tuples
                edges = {(j, k) for j in ed.keys() for k in ed[j].keys()}

                # Extract labels for nodes
                Lv = g.get_labels(purpose=self._graph_format)
                # and for edges
                Le = g.get_labels(purpose=self._graph_format, label_type="edge")

                # Produce all the neighborhoods and the distance pairs
                # up to the desired radius and maximum distance
                N, D, D_pair = g.produce_neighborhoods(self.r, purpose="dictionary",
                                                       with_distances=True, d=self.d)

                # Hash all the neighborhoods
                H = self._hash_neighborhoods(vertices, edges, Lv, Le, N, D_pair)

                if self._method_calling == 1:
                    for d in filterfalse(lambda x: x not in D,
                                         range(self.d+1)):
                        for (A, B) in D[d]:
                            for r in range(self.r+1):
                                key = (H[r, A], H[r, B])
                                keys = all_keys[r, d]
                                idx = keys.get(key, None)
                                if idx is None:
                                    idx = len(keys)
                                    keys[key] = idx
                                data[r, d][ng, idx] = data[r, d].get((ng, idx), 0) + 1

                elif self._method_calling == 3:
                    for d in filterfalse(lambda x: x not in D,
                                         range(self.d+1)):
                        for (A, B) in D[d]:
                            # Based on the edges of the bidirected graph
                            for r in range(self.r+1):
                                keys = all_keys[r, d]
                                fit_keys = self._fit_keys[r, d]
                                key = (H[r, A], H[r, B])
                                idx = fit_keys.get(key, None)
                                if idx is None:
                                    idx = keys.get(key, None)
                                    if idx is None:
                                        idx = len(keys) + len(fit_keys)
                                        keys[key] = idx
                                data[r, d][ng, idx] = data[r, d].get((ng, idx), 0) + 1
                ng += 1
            if ng == 0:
                raise ValueError('parsed input is empty')

            if self._method_calling == 1:
                # A feature matrix for all levels
                M = dict()

                for (key, d) in filterfalse(lambda a: len(a[1]) == 0,
                                            iteritems(data)):
                    indexes, data = zip(*iteritems(d))
                    rows, cols = zip(*indexes)
                    M[key] = csr_matrix((data, (rows, cols)), shape=(ng, len(all_keys[key])),
                                        dtype=np.int64)
                self._fit_keys = all_keys
                self._ngx = ng

            elif self._method_calling == 3:
                # A feature matrix for all levels
                M = dict()

                for (key, d) in filterfalse(lambda a: len(a[1]) == 0,
                                            iteritems(data)):
                    indexes, data = zip(*iteritems(d))
                    rows, cols = zip(*indexes)
                    M[key] = csr_matrix((data, (rows, cols)),
                                        shape=(ng, len(all_keys[key]) + len(self._fit_keys[key])),
                                        dtype=np.int64)

                self._ngy = ng

            return M

    def transform(self, X, y=None):
        """Calculate the kernel matrix, between given and fitted dataset.

        Parameters
        ----------
        X : iterable
            Each element must be an iterable with at most three features and at
            least one. The first that is obligatory is a valid graph structure
            (adjacency matrix or edge_dictionary) while the second is
            node_labels and the third edge_labels (that fitting the given graph
            format).

        y : Object, default=None
            Ignored argument, added for the pipeline.

        Returns
        -------
        K : numpy array, shape = [n_targets, n_input_graphs]
            corresponding to the kernel matrix, a calculation between
            all pairs of graphs between target an features

        """
        self._method_calling = 3
        # Check is fit had been called
        check_is_fitted(self, ['X'])

        # Input validation and parsing
        if X is None:
            raise ValueError('transform input cannot be None')
        else:
            Y = self.parse_input(X)

        try:
            check_is_fitted(self, ['_X_level_norm_factor'])
        except NotFittedError:
            self._X_level_norm_factor = \
                {key: np.array(M.power(2).sum(-1))
                 for (key, M) in iteritems(self.X)}

        N = self._X_level_norm_factor
        S = np.zeros(shape=(self._ngy, self._ngx))
        for (key, Mp) in filterfalse(lambda x: x[0] not in self.X,
                                     iteritems(Y)):
            M = self.X[key]
            K = M.dot(Mp.T[:M.shape[1]]).toarray().T
            S += np.nan_to_num(K / np.sqrt(np.outer(np.array(Mp.power(2).sum(-1)), N[key])))

        self._Y = Y
        self._is_transformed = True
        if self.normalize:
            S /= np.sqrt(np.outer(*self.diagonal()))
        return S

    def fit_transform(self, X):
        """Fit and transform, on the same dataset.

        Parameters
        ----------
        X : iterable
            Each element must be an iterable with at most three features and at
            least one. The first that is obligatory is a valid graph structure
            (adjacency matrix or edge_dictionary) while the second is
            node_labels and the third edge_labels (that fitting the given graph
            format). If None the kernel matrix is calculated upon fit data.
            The test samples.

        Returns
        -------
        K : numpy array, shape = [n_input_graphs, n_input_graphs]
            corresponding to the kernel matrix, a calculation between
            all pairs of graphs between target an features

        """
        self._method_calling = 2
        self.fit(X)

        S, N = np.zeros(shape=(self._ngx, self._ngx)), dict()
        for (key, M) in iteritems(self.X):
            K = M.dot(M.T).toarray()
            K_diag = K.diagonal()
            N[key] = K_diag
            S += np.nan_to_num(K / np.sqrt(np.outer(K_diag, K_diag)))

        self._X_level_norm_factor = N

        if self.normalize:
            return S / len(self.X)
        else:
            return S

    def diagonal(self):
        """Calculate the kernel matrix diagonal of the fitted data.

        Static. Added for completeness.

        Parameters
        ----------
        None.

        Returns
        -------
        X_diag : int
            Always equal with r*d.

        Y_diag : int
            Always equal with r*d.

        """
        # constant based on normalization of krd
        check_is_fitted(self, ['X'])
        try:
            check_is_fitted(self, ['_X_diag'])
        except NotFittedError:
            # Calculate diagonal of X
            self._X_diag = len(self.X)

        try:
            check_is_fitted(self, ['_Y'])
            return self._X_diag, len(self._Y)
        except NotFittedError:
            return self._X_diag

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
        H, sel = dict(), sorted(list(edges))
        for v in vertices:
            re, lv, le = sel, Lv, Le
            for radius in range(self.r, -1, -1):
                sub_vertices = sorted(N[radius][v])
                re = {(i, j) for (i, j) in re
                      if i in sub_vertices and j in sub_vertices}
                lv = {v: lv[v] for v in sub_vertices}
                le = {e: le[e] for e in edges}
                H[radius, v] = hash_graph(D_pair, sub_vertices, re, lv, le)
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
        label = "|".join(sorted([str(D[(i, j)]) + ',' + str(glv[j])
                                 for j in vertices if (i, j) in D]))
        encoding += label + "."
        Lv[i] = label

    encoding = encoding[:-1]+":"

    # Expand to labels for edges
    for (i, j) in edges:
        encoding += Lv[i] + ',' + Lv[j] + ',' + str(gle[(i, j)]) + "_"

    # Arash Partov hashing, as in the original
    # implementation of NSPK.
    return APHash(encoding)
