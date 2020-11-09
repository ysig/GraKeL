"""Shortest path kernel as defined in :cite:`borgwardt2005shortest`."""
# Author: Ioannis Siglidis <y.siglidis@gmail.com>
# License: BSD 3 clause
import collections
import warnings

import numpy as np

from sklearn.exceptions import NotFittedError
from sklearn.utils.validation import check_is_fitted

from grakel.graph import Graph
from grakel.kernels import Kernel


class ShortestPathAttr(Kernel):
    r"""The shortest path kernel for attributes.

    The Graph labels are considered as attributes.
    The computational efficiency is decreased to :math:`O(|V|^4)`
    See :cite:`borgwardt2005shortest`.

    Parameters
    ----------
    algorithm_type : str, default={"dijkstra", "floyd_warshall", "auto"}
        Apply the dijkstra or floyd_warshall algorithm for calculating
        shortest path, or chose automatically ("auto") based on the
        current graph format ("auto").

    metric : function, default=:math:`f(x,y)=\sum_{i}x_{i}*y_{i}`,
        The metric applied between attributes of the graph labels.
        The user must provide a metric based on the format of the provided
        labels (considered as attributes).

    Attributes
    ----------
    X : list
        A list of tuples, consisting of shortest path matrices
        and their feature vectors.

    """

    def __init__(self, n_jobs=None,
                 normalize=False,
                 verbose=False,
                 algorithm_type="auto",
                 metric=np.dot):
        """Initialise a `shortest_path_attr` kernel."""
        super(ShortestPathAttr, self).__init__(
            n_jobs=n_jobs, normalize=normalize, verbose=verbose)

        self.algorithm_type = algorithm_type
        self.metric = metric
        self._initialized.update({"algorithm_type": False, "metric": False})

    def initialize(self):
        """Initialize all transformer arguments, needing initialization."""
        super(ShortestPathAttr, self).initialize()
        if not self._initialized["algorithm_type"]:
            if self.algorithm_type == "auto":
                self._graph_format = "auto"
            elif self.algorithm_type == "floyd_warshall":
                self._graph_format = "adjacency"
            elif self.algorithm_type == "dijkstra":
                self._graph_format = "dictionary"
            else:
                raise ValueError('Unsupported value ' +
                                 str(self.algorithm_type) +
                                 ' for "algorithm_type"')
            self._initialized["algorithm_type"] = True

        if not self._initialized["metric"]:
            if not callable(self.metric):
                raise TypeError('"metric" must be callable')
            self._initialized["metric"] = True

    def parse_input(self, X):
        """Parse and create features for the `shortest_path` kernel.

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
        sp_attr_tup : list
            A list of tuples of shortest path matrices and tehir attributes.

        """
        if not isinstance(X, collections.Iterable):
            raise TypeError('input must be an iterable\n')
        else:
            sp_attr_tup = list()
            ni = 0
            for (i, x) in enumerate(iter(X)):
                is_iter = isinstance(x, collections.Iterable)
                if is_iter:
                    x = list(x)
                if is_iter and len(x) in [0, 2, 3]:
                    if len(x) == 0:
                        warnings.warn('Ignoring empty element' +
                                      ' on index: '+str(i))
                        continue
                    else:
                        S, L = Graph(
                            x[0], x[1], {},
                            self._graph_format).build_shortest_path_matrix(
                                self.algorithm_type)
                elif type(x) is Graph:
                    S, L = x.build_shortest_path_matrix(self.algorithm_type)
                else:
                    raise TypeError('each element of X must be either a ' +
                                    'graph or an iterable with at least 2 ' +
                                    'and at most 3 elements\n')

                sp_attr_tup.append((S, L))
                ni += 1

            if ni == 0:
                raise ValueError('parsed input is empty')

            return sp_attr_tup

    def pairwise_operation(self, x, y):
        """Calculate shortests paths on attributes.

        Parameters
        ----------
        x, y : tuple
            Tuples of shortest path matrices and their attribute
            dictionaries.

        Returns
        -------
        kernel : number
            The kernel value.

        """
        # Initialise
        Sx, phi_x = x
        Sy, phi_y = y
        kernel = 0
        dimx = Sx.shape[0]
        dimy = Sy.shape[0]
        for i in range(dimx):
            for j in range(dimx):
                if i == j:
                    continue
                for k in range(dimy):
                    for m in range(dimy):
                        if k == m:
                            continue
                        if (Sx[i, j] == Sy[k, m] and
                                Sx[i, j] != float('Inf')):
                            kernel += self.metric(phi_x[i], phi_y[k]) *\
                                self.metric(phi_x[j], phi_y[m])

        return kernel


class ShortestPath(Kernel):
    r"""The shortest path kernel class.

    See :cite:`borgwardt2005shortest`.

    Parameters
    ----------
    algorithm_type : str, default={"dijkstra", "floyd_warshall", "auto"}
        Apply the dijkstra or floyd_warshall algorithm for calculating
        shortest path, or chose automatically ("auto") based on the
        current graph format ("auto").

    with_labels : bool, default=True, case_of_existence=(as_attributes==True)
        Calculate shortest path using graph labels.

    Attributes
    ----------
    X : dict
        A dictionary of pairs between each input graph and a bins where the
        sampled graphlets have fallen.

    _with_labels : bool
        Defines if the shortest path kernel considers also labels.

    _enum : dict
        A dictionary of graph bins holding pynauty objects

    _lt : str
        A label type needed for build shortest path function.

    _lhash : str
        A function for hashing labels, shortest paths.

    _nx : int
        Holds the number of sampled X graphs.

    _ny : int
        Holds the number of sampled Y graphs.

    _X_diag : np.array, shape=(_nx, 1)
        Holds the diagonal of X kernel matrix in a numpy array, if calculated
        (`fit_transform`).

    _phi_X : np.array, shape=(_nx, len(_graph_bins))
        Holds the features of X in a numpy array, if calculated.
        (`fit_transform`).

    Complexity
    ----------
    :math:`O(n*N*|\cup_{i}L_{i}|^{2})`, where :math:`n` the number of graph,
    :math:`N` the number of vertices of the **biggest** Graph and
    :math:`|\cup_{i}L_{i}|` the number of all distinct labels.

    """

    _graph_bins = dict()

    def __init__(self, n_jobs=None,
                 normalize=False,
                 verbose=False,
                 with_labels=True,
                 algorithm_type="auto"):
        """Initialize a `shortest_path` kernel."""
        super(ShortestPath, self).__init__(
            n_jobs=n_jobs, normalize=normalize, verbose=verbose)

        self.with_labels = with_labels
        self.algorithm_type = algorithm_type
        self._initialized.update({"with_labels": False, "algorithm_type": False})

    def initialize(self):
        """Initialize all transformer arguments, needing initialization."""
        if not self._initialized["n_jobs"]:
            if self.n_jobs is not None:
                warnings.warn('no implemented parallelization for ShortestPath')
            self._initialized["n_jobs"] = True

        if not self._initialized["algorithm_type"]:
            if self.algorithm_type == "auto":
                self._graph_format = "auto"
            elif self.algorithm_type == "floyd_warshall":
                self._graph_format = "adjacency"
            elif self.algorithm_type == "dijkstra":
                self._graph_format = "dictionary"
            else:
                raise ValueError('Unsupported "algorithm_type"')

        if not self._initialized["with_labels"]:
            if self.with_labels:
                self._lt = "vertex"
                self._lhash = lhash_labels
                self._decompose_input = decompose_input_labels
            else:
                self._lt = "none"
                self._lhash = lhash
                self._decompose_input = decompose_input

    def transform(self, X):
        """Calculate the kernel matrix, between given and fitted dataset.

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
        K : numpy array, shape = [n_targets, n_input_graphs]
            corresponding to the kernel matrix, a calculation between
            all pairs of graphs between target an features

        """
        self._method_calling = 3
        # Check is fit had been called
        check_is_fitted(self, ['X', '_nx', '_enum'])

        # Input validation and parsing
        if X is None:
            raise ValueError('transform input cannot be None')
        else:
            Y = self.parse_input(X)

        # Transform - calculate kernel matrix
        try:
            check_is_fitted(self, ['_phi_X'])
            phi_x = self._phi_X
        except NotFittedError:
            phi_x = np.zeros(shape=(self._nx, len(self._enum)))
            for i in self.X.keys():
                for j in self.X[i].keys():
                    phi_x[i, j] = self.X[i][j]
            self._phi_X = phi_x

        phi_y = np.zeros(shape=(self._ny, len(self._enum) + len(self._Y_enum)))
        for i in Y.keys():
            for j in Y[i].keys():
                phi_y[i, j] = Y[i][j]

        # store _phi_Y for independent (of normalization arg diagonal-calls)
        self._phi_Y = phi_y
        km = np.dot(phi_y[:, :len(self._enum)], phi_x.T)
        self._is_transformed = True
        if self.normalize:
            X_diag, Y_diag = self.diagonal()
            return km / np.sqrt(np.outer(Y_diag, X_diag))
        else:
            return km

    def diagonal(self):
        """Calculate the kernel matrix diagonal for fitted data.

        A funtion called on transform on a seperate dataset to apply
        normalization on the exterior.

        Parameters
        ----------
        None.

        Returns
        -------
        X_diag : np.array
            The diagonal of the kernel matrix, of the fitted data.
            This consists of kernel calculation for each element with itself.

        Y_diag : np.array
            The diagonal of the kernel matrix, of the transformed data.
            This consists of kernel calculation for each element with itself.

        """
        # Check is fit and transform had been called
        try:
            check_is_fitted(self, ['_phi_X'])
        except NotFittedError:
            check_is_fitted(self, ['X'])
            # calculate feature matrices.
            phi_x = np.zeros(shape=(self._nx, len(self._enum)))

            for i in self.X.keys():
                for j in self.X[i].keys():
                    phi_x[i, j] = self.X[i][j]
                    # Transform - calculate kernel matrix
            self._phi_X = phi_x

        try:
            check_is_fitted(self, ['X_diag'])
        except NotFittedError:
            # Calculate diagonal of X
            self._X_diag = np.sum(np.square(self._phi_X), axis=1)
            self._X_diag = np.reshape(self._X_diag, (self._X_diag.shape[0], 1))

        try:
            check_is_fitted(self, ['_phi_Y'])
            # Calculate diagonal of Y
            Y_diag = np.sum(np.square(self._phi_Y), axis=1)
            return self._X_diag, Y_diag
        except NotFittedError:
            return self._X_diag

    def fit_transform(self, X, y=None):
        """Fit and transform, on the same dataset.

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
        self._method_calling = 2
        self.fit(X)

        # calculate feature matrices.
        phi_x = np.zeros(shape=(self._nx, len(self._enum)))

        for i in self.X.keys():
            for j in self.X[i].keys():
                phi_x[i, j] = self.X[i][j]

        # Transform - calculate kernel matrix
        self._phi_X = phi_x
        km = np.dot(phi_x, phi_x.T)

        self._X_diag = np.diagonal(km)
        if self.normalize:
            return np.divide(km, np.sqrt(np.outer(self._X_diag, self._X_diag)))
        else:
            return km

    def parse_input(self, X):
        """Parse and create features for "shortest path" kernel.

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
        sp_counts : dict
            A dictionary that for each vertex holds the counts of shortest path
            tuples.

        """
        if not isinstance(X, collections.Iterable):
            raise TypeError('input must be an iterable\n')
            # Not a dictionary
        else:
            i = -1
            sp_counts = dict()
            if self._method_calling == 1:
                self._enum = dict()
            elif self._method_calling == 3:
                self._Y_enum = dict()
            for (idx, x) in enumerate(iter(X)):
                is_iter = isinstance(x, collections.Iterable)
                if is_iter:
                    x = list(x)
                if is_iter and (len(x) == 0 or
                                (len(x) == 1 and not self.with_labels) or
                                len(x) in [2, 3]):
                    if len(x) == 0:
                        warnings.warn('Ignoring empty element on index: '
                                      + str(idx))
                        continue
                    elif len(x) == 1:
                        spm_data = Graph(x[0], {}, {}, self._graph_format
                                         ).build_shortest_path_matrix(self.algorithm_type,
                                                                      labels=self._lt)
                    else:
                        spm_data = Graph(x[0], x[1], {}, self._graph_format
                                         ).build_shortest_path_matrix(self.algorithm_type,
                                                                      labels=self._lt)
                elif type(x) is Graph:
                    spm_data = x.build_shortest_path_matrix(self.algorithm_type, labels=self._lt)
                else:
                    raise TypeError('each element of X must have at least' +
                                    ' one and at most 3 elements\n')
                i += 1

                S, L = self._decompose_input(spm_data)
                sp_counts[i] = dict()
                for u in range(S.shape[0]):
                    for v in range(S.shape[1]):
                        if u == v or S[u, v] == float("Inf"):
                            continue
                        label = self._lhash(S, u, v, *L)
                        if label not in self._enum:
                            if self._method_calling == 1:
                                idx = len(self._enum)
                                self._enum[label] = idx
                            elif self._method_calling == 3:
                                if label not in self._Y_enum:
                                    idx = len(self._enum) + len(self._Y_enum)
                                    self._Y_enum[label] = idx
                                else:
                                    idx = self._Y_enum[label]
                        else:
                            idx = self._enum[label]
                        if idx in sp_counts[i]:
                            sp_counts[i][idx] += 1
                        else:
                            sp_counts[i][idx] = 1

            if i == -1:
                raise ValueError('parsed input is empty')

            if self._method_calling == 1:
                self._nx = i+1
            elif self._method_calling == 3:
                self._ny = i+1
            return sp_counts


def lhash(S, u, v, *args):
    return S[u, v]


def decompose_input(a):
    return (a, [])


def lhash_labels(S, u, v, *args):
    return (args[0][u], args[0][v], S[u, v])


def decompose_input_labels(args):
    return (args[0], args[1:])
