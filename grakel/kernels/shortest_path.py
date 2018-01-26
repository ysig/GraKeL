"""Shortest path kernel as defined in :cite:`Borgwardt2005ShortestpathKO`."""
import collections
import warnings

import numpy as np

from sklearn.exceptions import NotFittedError
from sklearn.utils.validation import check_is_fitted

from grakel.graph import graph
from grakel.kernels import kernel


class shortest_path_attr(kernel):
    r"""The shortest path kernel for attributes.

    See :cite:`Borgwardt2005ShortestpathKO`.

    Parameters
    ----------
    algorithm_type : str, default={"dijkstra", "floyd_warshall", "auto"}
        Apply the dijkstra or floyd_warshall algorithm for calculating
        shortest path, or chose automatically ("auto") based on the
        current graph format ("auto").

    as_attributes : bool, default=False
        The labels are considered as attributes. The computational
        efficiency is decreased to :math:`O(|V|^4)`

    attribute_kernel : function, default=:math:`f(x,y)=\sum_{i}x_{i}*y_{i}`,
        The kernel applied between attributes of the graph labels.
        The user must provide a kernel based on the format of the provided
        labels (considered as attributes).

    Attributes
    ----------
    X : list
        A list of tuples, consisting of shortest path matrices
        and their feature vectors.

    """

    def __init__(self, **kargs):
        """Initialise a subtree_wl kernel."""
        self._valid_parameters |= {"algorithm_type",
                                   "as_attributes",
                                   "attribute_kernel"}
        super(shortest_path_attr, self).__init__(**kargs)

        self._algorithm_type = kargs.get("algorithm_type", "auto")

        if self._algorithm_type == "auto":
            self._graph_format = "auto"
        elif self._algorithm_type == "floyd_warshall":
            self._algorithm_type = "adjacency"
        elif self._algorithm_type == "dijkstra":
            self._algorithm_type = "dictionary"
        else:
            raise ValueError('Unsupported "algorithm_type"')

        self._attribute_kernel = kargs.get(
            "attribute_kernel", lambda x, y: np.dot(x, y))

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
        sp_attr_tup : list
            A list of tuples of shortest path matrices and tehir attributes.

        """
        if not isinstance(X, collections.Iterable):
            raise ValueError('input must be an iterable\n')
        else:
            i = 0
            sp_attr_tup = list()
            for x in iter(X):
                if len(x) == 0:
                    warnings.warn('Ignoring empty element on index: '+str(i))
                if len(x) == 1:
                    if type(x) is graph:
                        S, L = x.build_shortest_path_matrix(
                                    self._algorithm_type)
                    else:
                        warnings.warn(
                            'Ignoring empty element on index: '
                            + str(i) + '\nLabels must be provided.')
                    i += 1
                elif len(x) in [2, 3]:
                    S, L = graph(
                        x[0], x[1], {},
                        self._graph_format).build_shortest_path_matrix(
                            self._algorithm_type)
                    i += 1
                else:
                    raise ValueError('each element of X must have at least' +
                                     ' one and at most 3 elements\n')

                sp_attr_tup.append((S, L))
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
                                kernel += \
                                    self._attribute_kernel(
                                        phi_x[i], phi_y[k]) *\
                                    self. attribute_kernel(
                                        phi_x[j], phi_y[m])

            return kernel


class shortest_path(kernel):
    r"""The shortest path kernel class.

    See :cite:`Borgwardt2005ShortestpathKO`.

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

    def __init__(self, **kargs):
        """Initialise a subtree_wl kernel."""
        self._valid_parameters |= {"with_labels",
                                   "algorithm_type"}
        super(shortest_path, self).__init__(**kargs)

        if kargs.get("with_labels", True):
            self._lt = "vertex"
            self._lhash = lambda S, u, v, *args: \
                (args[0][u], args[0][v], S[u, v])
        else:
            self._lt = "none"
            self._lhash = lambda S, u, v, *args: S[u, v]

        self._algorithm_type = kargs.get("algorithm_type", "auto")

        if self._algorithm_type == "auto":
            self._graph_format = "auto"
        elif self._algorithm_type == "floyd_warshall":
            self._algorithm_type = "adjacency"
        elif self._algorithm_type == "dijkstra":
            self._algorithm_type = "dictionary"
        else:
            raise ValueError('Unsupported "algorithm_type"')

    def transform(self, X):
        """Calculate the kernel matrix, between given and fitted dataset.

        Paramaters
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
            check_is_fitted(self, ['phi_X'])
            phi_x = self._phi_X
        except NotFittedError:
            phi_x = np.zeros(shape=(self._nx, len(self._enum)))
            for i in self.X.keys():
                for j in self.X[i].keys():
                    phi_x[i, j] = self.X[i][j]

        phi_y = np.zeros(shape=(self._ny, len(self._enum) + len(self._Y_enum)))
        for i in Y.keys():
            for j in Y[i].keys():
                phi_y[i, j] = Y[i][j]

        # store _phi_Y for independent (of normalization arg diagonal-calls)
        self._phi_Y = phi_y
        km = np.dot(phi_y[:, :len(self._enum)], phi_x.T)
        if self._normalize:
            X_diag, Y_diag = self.diagonal()
            km /= np.sqrt(np.dot(Y_diag, X_diag.T))
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
        check_is_fitted(self, ['_phi_X', '_phi_Y'])
        try:
            check_is_fitted(self, ['X_diag'])
        except NotFittedError:
            # Calculate diagonal of X
            self._X_diag = np.sum(np.square(self._phi_X), axis=1)
            self._X_diag = np.reshape(self._X_diag,
                                      (self._X_diag.shape[0], 1))
        # Calculate diagonal of Y
        Y_diag = np.sum(np.square(self._phi_Y), axis=1)
        return self._X_diag, np.reshape(Y_diag, (Y_diag.shape[0], 1))

    def fit_transform(self, X):
        """Fit and transform, on the same dataset.

        Paramaters
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

        self._X_diag = np.diagonal(km).reshape(km.shape[0], 1)
        if self._normalize:
            return np.divide(km,
                             np.sqrt(np.multiply(self._X_diag.T,
                                                 self._X_diag)))
        else:
            return km
        return km

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
        sp_counts : dict
            A dictionary that for each vertex holds the counts of shortest_path
            tuples.

        """
        if not isinstance(X, collections.Iterable):
            raise ValueError('input must be an iterable\n')
            # Not a dictionary
        else:
            i = -1
            sp_counts = dict()
            if self._method_calling == 1:
                self._enum = dict()
            elif self._method_calling == 3:
                self._Y_enum = dict()
            for x in iter(X):
                if len(x) == 0:
                    warnings.warn('Ignoring empty element on index: '+str(i))
                if len(x) == 1:
                    if type(x) is graph:
                        S, *L = x.build_shortest_path_matrix(
                                    self._algorithm_type,
                                    labels=self._lt)
                    else:
                        S, *L = graph(
                            x[0], {}, {},
                            self._graph_format).build_shortest_path_matrix(
                                self._algorithm_type,
                                labels=self._lt)
                    i += 1
                elif len(x) in [2, 3]:
                    S, *L = graph(
                        x[0], x[1], {},
                        self._graph_format).build_shortest_path_matrix(
                            self._algorithm_type,
                            labels=self._lt)
                    i += 1
                else:
                    raise ValueError('each element of X must have at least' +
                                     ' one and at most 3 elements\n')

                sp_counts[i] = dict()
                for u in range(S.shape[0]):
                    for v in range(S.shape[0]):
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
