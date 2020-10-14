"""The weisfeiler lehman optimal assignment kernel :cite:`kriege2016valid`."""
# Author: Giannis Nikolentzos <nikolentzos@lix.polytechnique.fr>
# License: BSD 3 clause
import collections
import warnings

import numpy as np

from sklearn.exceptions import NotFittedError
from sklearn.utils.validation import check_is_fitted

from grakel.graph import Graph
from grakel.kernels import Kernel

from scipy.sparse import lil_matrix

# Python 2/3 cross-compatibility import
from six import iteritems
from six import itervalues


class WeisfeilerLehmanOptimalAssignment(Kernel):
    """Compute the Weisfeiler Lehman Optimal Assignment Kernel.

    See :cite:`kriege2016valid`.

    Parameters
    ----------
    n_iter : int, default=5
        The number of iterations.

    Attributes
    ----------
    X : dict
     Holds a list of fitted subkernel modules.

    sparse : bool
        Defines if the data will be stored in a sparse format.
        Sparse format is slower, but less memory consuming and in some cases the only solution.

    _nx : number
        Holds the number of inputs.

    _n_iter : int
        Holds the number, of iterations.

    _hierarchy : dict
        A hierarchy produced by the WL relabeling procedure.

    _inv_labels : dict
        An inverse dictionary, used for relabeling on each iteration.

    """

    _graph_format = "dictionary"

    def __init__(self, n_jobs=None, verbose=False,
                 normalize=False, n_iter=5, sparse=False):
        """Initialise a `weisfeiler_lehman` kernel."""
        super(WeisfeilerLehmanOptimalAssignment, self).__init__(
            n_jobs=n_jobs, verbose=verbose, normalize=normalize)

        self.n_iter = n_iter
        self.sparse = sparse
        self._initialized.update({"n_iter": False, 'sparse': True})

    def initialize(self):
        """Initialize all transformer arguments, needing initialization."""
        super(WeisfeilerLehmanOptimalAssignment, self).initialize()

        if not self._initialized["n_iter"]:
            if type(self.n_iter) is not int or self.n_iter <= 0:
                raise TypeError("'n_iter' must be a positive integer")
            self._n_iter = self.n_iter + 1
            self._initialized["n_iter"] = True
        if not self._initialized["sparse"]:
            if self.sparse not in [False, True]:
                TypeError('sparse could be False, True')
            self._initialized["sparse"] = False

    def parse_input(self, X):
        """Parse input for weisfeiler lehman optimal assignment.

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
        Hs : numpy array, shape = [n_input_graphs, hierarchy_size]
            An array where the rows contain the histograms of the graphs.

        """
        if self._method_calling not in [1, 2]:
            raise ValueError('method call must be called either from fit ' +
                             'or fit-transform')
        elif hasattr(self, '_X_diag'):
            # Clean _X_diag value
            delattr(self, '_X_diag')

        # Input validation and parsing
        if not isinstance(X, collections.Iterable):
            raise TypeError('input must be an iterable\n')
        else:
            nx = 0
            Gs_ed, L, distinct_values = dict(), dict(), set()
            for (idx, x) in enumerate(iter(X)):
                is_iter = isinstance(x, collections.Iterable)
                if is_iter:
                    x = list(x)
                if is_iter and (len(x) == 0 or len(x) >= 2):
                    if len(x) == 0:
                        warnings.warn('Ignoring empty element on index: '
                                      + str(idx))
                        continue
                    else:
                        if len(x) > 2:
                            extra = tuple()
                            if len(x) > 3:
                                extra = tuple(x[3:])
                            x = Graph(x[0], x[1], x[2], graph_format=self._graph_format)
                            extra = (x.get_labels(purpose=self._graph_format,
                                                  label_type="edge", return_none=True), ) + extra
                        else:
                            x = Graph(x[0], x[1], {}, graph_format=self._graph_format)
                            extra = tuple()

                elif type(x) is Graph:
                    x.desired_format(self._graph_format)
                else:
                    raise TypeError('each element of X must be either a ' +
                                    'graph object or a list with at least ' +
                                    'a graph like object and node labels ' +
                                    'dict \n')
                Gs_ed[nx] = x.get_edge_dictionary()
                L[nx] = x.get_labels(purpose="dictionary")
                distinct_values |= set(itervalues(L[nx]))
                nx += 1
            if nx == 0:
                raise ValueError('parsed input is empty')

        # Save the number of "fitted" graphs.
        self._nx = nx

        # Initialize hierarchy
        self._hierarchy = dict()
        self._hierarchy['root'] = dict()
        self._hierarchy['root']['parent'] = None
        self._hierarchy['root']['children'] = list()
        self._hierarchy['root']['w'] = 0
        self._hierarchy['root']['omega'] = 0

        # get all the distinct values of current labels
        WL_labels_inverse = dict()

        # assign a number to each label
        label_count = 0
        for dv in sorted(list(distinct_values)):
            WL_labels_inverse[dv] = label_count
            self._insert_into_hierarchy(label_count, 'root')
            label_count += 1

        # Initalize an inverse dictionary of labels for all iterations
        self._inv_labels = dict()
        self._inv_labels[0] = WL_labels_inverse

        for j in range(nx):
            new_labels = dict()
            for k in L[j].keys():
                new_labels[k] = WL_labels_inverse[L[j][k]]
            L[j] = new_labels

        for i in range(1, self._n_iter):
            new_previous_label_set, WL_labels_inverse, L_temp = set(), dict(), dict()
            for j in range(nx):
                # Find unique labels and sort
                # them for both graphs
                # Keep for each node the temporary
                L_temp[j] = dict()
                for v in Gs_ed[j].keys():
                    credential = str(L[j][v]) + "," + \
                        str(sorted([L[j][n] for n in Gs_ed[j][v].keys()]))
                    L_temp[j][v] = credential
                    new_previous_label_set.add((credential, L[j][v]))

            label_list = sorted(list(new_previous_label_set), key=lambda tup: tup[0])
            for dv, previous_label in label_list:
                WL_labels_inverse[dv] = label_count
                self._insert_into_hierarchy(label_count, previous_label)
                label_count += 1

            # Recalculate labels
            for j in range(nx):
                new_labels = dict()
                for k in L_temp[j].keys():
                    new_labels[k] = WL_labels_inverse[L_temp[j][k]]
                L[j] = new_labels
            self._inv_labels[i] = WL_labels_inverse

        # Compute the vector representation of each graph
        if self.sparse:
            Hs = lil_matrix((nx, len(self._hierarchy)))
        else:
            Hs = np.zeros((nx, len(self._hierarchy)))
        for j in range(nx):
            for k in L[j].keys():
                current_label = L[j][k]
                while self._hierarchy[current_label]['parent'] is not None:
                    Hs[j, current_label] += self._hierarchy[current_label]['omega']
                    current_label = self._hierarchy[current_label]['parent']

        return Hs

    def _insert_into_hierarchy(self, label, previous_label):
        """Inserts a label into the hierarchy.

        Parameters
        ----------
        label : int
            The label to insert into the hierarchy.

        previous_label : int
            The previous label of the node.

        """
        self._hierarchy[label] = dict()
        self._hierarchy[label]['parent'] = previous_label
        self._hierarchy[label]['children'] = list()
        self._hierarchy[label]['w'] = self._hierarchy[previous_label]['w'] + 1
        self._hierarchy[label]['omega'] = 1
        self._hierarchy[previous_label]['children'].append(label)

    def fit_transform(self, X, y=None):
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

        y : Object, default=None
            Ignored argument, added for the pipeline.

        Returns
        -------
        K : numpy array, shape = [n_targets, n_input_graphs]
            corresponding to the kernel matrix, a calculation between
            all pairs of graphs between target an features

        """
        self._method_calling = 2
        self._is_transformed = False
        self.initialize()
        if X is None:
            raise ValueError('transform input cannot be None')
        else:
            self.X = self.parse_input(X)

        # Compute the histogram intersection kernel
        K = np.zeros((self._nx, self._nx))
        if self.sparse:
            for i in range(self._nx):
                for j in range(i, self._nx):
                    K[i, j] = np.sum(self.X[i, :].minimum(self.X[j, :]))
                    K[j, i] = K[i, j]
        else:
            for i in range(self._nx):
                for j in range(i, self._nx):
                    K[i, j] = np.sum(np.min(self.X[np.ix_([i, j]), :], axis=1))
                    K[j, i] = K[i, j]

        self._X_diag = np.diagonal(K)
        if self.normalize:
            old_settings = np.seterr(divide='ignore')
            K = np.nan_to_num(np.divide(K, np.sqrt(np.outer(self._X_diag, self._X_diag))))
            np.seterr(**old_settings)
        return K

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
        check_is_fitted(self, ['X', '_nx', '_hierarchy', '_inv_labels'])

        # Input validation and parsing
        if X is None:
            raise ValueError('transform input cannot be None')
        else:
            if not isinstance(X, collections.Iterable):
                raise ValueError('input must be an iterable\n')
            else:
                nx = 0
                distinct_values = set()
                Gs_ed, L = dict(), dict()
                for (i, x) in enumerate(iter(X)):
                    is_iter = isinstance(x, collections.Iterable)
                    if is_iter:
                        x = list(x)
                    if is_iter and len(x) in [0, 2, 3]:
                        if len(x) == 0:
                            warnings.warn('Ignoring empty element on index: '
                                          + str(i))
                            continue

                        elif len(x) in [2, 3]:
                            x = Graph(x[0], x[1], {}, self._graph_format)
                    elif type(x) is Graph:
                        x.desired_format("dictionary")
                    else:
                        raise ValueError('each element of X must have at ' +
                                         'least one and at most 3 elements\n')
                    Gs_ed[nx] = x.get_edge_dictionary()
                    L[nx] = x.get_labels(purpose="dictionary")

                    # Hold all the distinct values
                    distinct_values |= set(
                        v for v in itervalues(L[nx])
                        if v not in self._inv_labels[0])
                    nx += 1
                if nx == 0:
                    raise ValueError('parsed input is empty')

        # get all the distinct values of new labels
        WL_labels_inverse = dict()

        # assign a number to each label
        label_count = sum([len(self._inv_labels[i]) for i in range(len(self._inv_labels))])
        for dv in sorted(list(distinct_values)):
            WL_labels_inverse[dv] = label_count
            self._insert_into_hierarchy(label_count, 'root')
            label_count += 1

        for j in range(nx):
            new_labels = dict()
            for (k, v) in iteritems(L[j]):
                if v in self._inv_labels[0]:
                    new_labels[k] = self._inv_labels[0][v]
                else:
                    new_labels[k] = WL_labels_inverse[v]
            L[j] = new_labels

        for i in range(1, self._n_iter):
            L_temp, new_previous_label_set = dict(), set()
            for j in range(nx):
                # Find unique labels and sort them for both graphs
                # Keep for each node the temporary
                L_temp[j] = dict()
                for v in Gs_ed[j].keys():
                    credential = str(L[j][v]) + "," + \
                        str(sorted([L[j][n] for n in Gs_ed[j][v].keys()]))
                    L_temp[j][v] = credential
                    if credential not in self._inv_labels[i]:
                        new_previous_label_set.add((credential, L[j][v]))

            # Calculate the new label_set
            WL_labels_inverse = dict()
            if len(new_previous_label_set) > 0:
                for dv, previous_label in sorted(list(new_previous_label_set), key=lambda tup: tup[0]):
                    WL_labels_inverse[dv] = label_count
                    self._insert_into_hierarchy(label_count, previous_label)
                    label_count += 1

            # Recalculate labels
            for j in range(nx):
                new_labels = dict()
                for (k, v) in iteritems(L_temp[j]):
                    if v in self._inv_labels[i]:
                        new_labels[k] = self._inv_labels[i][v]
                    else:
                        new_labels[k] = WL_labels_inverse[v]
                L[j] = new_labels

        # Compute the vector representation of each graph
        if self.sparse:
            Hs = lil_matrix((nx, len(self._hierarchy)))
        else:
            Hs = np.zeros((nx, len(self._hierarchy)))
        for j in range(nx):
            for k in L[j].keys():
                current_label = L[j][k]
                while self._hierarchy[current_label]['parent'] is not None:
                    Hs[j, current_label] += self._hierarchy[current_label]['omega']
                    current_label = self._hierarchy[current_label]['parent']

        self.Y = Hs

        # Compute the histogram intersection kernel
        K = np.zeros((nx, self._nx))
        if self.sparse:
            for i in range(self._nx):
                for j in range(i, self._nx):
                    K[i, j] = np.sum(Hs[i, :self.X.shape[1]].minimum(self.X[j, :]))
        else:
            for i in range(nx):
                for j in range(self._nx):
                    K[i, j] = np.sum(np.min([Hs[i, :self.X.shape[1]], self.X[j, :]], axis=0))

        self._is_transformed = True
        if self.normalize:
            X_diag, Y_diag = self.diagonal()
            old_settings = np.seterr(divide='ignore')
            K = np.nan_to_num(np.divide(K, np.sqrt(np.outer(Y_diag, X_diag))))
            np.seterr(**old_settings)

        return K

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
        # Check if fit had been called
        check_is_fitted(self, ['X'])
        try:
            check_is_fitted(self, ['_X_diag'])
            if self._is_transformed:
                Y_diag = np.zeros(self.Y.shape[0])
                for i in range(self.Y.shape[0]):
                    Y_diag[i] = np.sum(np.min(self.Y[np.ix_([i, i]), :], axis=1))
        except NotFittedError:
            # Calculate diagonal of X
            if self._is_transformed:
                self._X_diag = np.zeros(self.X.shape[0])
                for i in range(self.X.shape[0]):
                    self._X_diag[i] = np.sum(np.min(self.X[np.ix_([i, i]), :], axis=1))

                Y_diag = np.zeros(self.Y.shape[0])
                for i in range(self.Y.shape[0]):
                    Y_diag[i] = np.sum(np.min(self.Y[np.ix_([i, i]), :], axis=1))
            else:
                # case sub kernel is only fitted
                self._X_diag = np.zeros(self.X.shape[0])
                for i in range(self.X.shape[0]):
                    self._X_diag[i] = np.sum(np.min(self.X[np.ix_([i, i]), :], axis=1))

        if self._is_transformed:
            return self._X_diag, Y_diag
        else:
            return self._X_diag


def efit(object, data):
    """Fit an object on data."""
    object.fit(data)


def efit_transform(object, data):
    """Fit-Transform an object on data."""
    return object.fit_transform(data)


def etransform(object, data):
    """Transform an object on data."""
    return object.transform(data)
