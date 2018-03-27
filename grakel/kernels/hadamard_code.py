"""The hadamard code kernel as defined in :cite:`icpram16`."""
import collections
import warnings

import numpy as np

from math import ceil
from numpy import log2

from scipy.linalg import hadamard

from sklearn.exceptions import NotFittedError
from sklearn.utils.validation import check_is_fitted

from grakel.graph import Graph
from grakel.kernels import kernel

# Python 2/3 cross-compatibility import
from six import iteritems
from six import itervalues

default_executor = lambda fn, *eargs, **ekargs: fn(*eargs, **ekargs)


class hadamard_code(kernel):
    """The simple Hadamard code kernel, as proposed in :cite:`icpram16`.

    Parameters
    ----------
    base_kernel : `grakel.kernels.kernel` or tuple
        If tuple it must consist of a valid kernel object and a
        dictionary of parameters. General parameters concerning
        normalization, concurrency, .. will be ignored, and the
        ones of given on `__init__` will be passed in case it is needed.

    rho : int, condition_of_appearance: hc_type=="shortened", default=-1
        The size of each single bit arrays. If -1 is chosen r is calculated as
        the biggest possible that satisfies an equal division.

    L : int, condition_of_appearance: hc_type=="shortened", default=4
        The number of bytes to store the bitarray of each label.

    niter : int, default=5
        The number of iterations.

    Attributes
    ----------
    _niter : int
        The number of iterations.

    _shortened : bool
        A flag signifying if the shortened version of the algorithm
        will be used.

    _add : function
        A function setted relevant to the version of the algorithm
        for adding hashed labels.

    _get : function
        A function setted relevant to the version of the algorithm
        for getting a label element.

    _base_kernel : function
        A void function that initializes a base kernel object.

    """

    _graph_format = "auto"

    def __init__(self, executor=default_executor, verbose=False,
                 normalize=False, niter=5, base_kernel=None):
        """Initialise a `hadamard_code` kernel."""
        super(hadamard_code, self).__init__(
            executor=executor, verbose=verbose, normalize=normalize)

        self.niter = niter
        self.base_kernel = base_kernel
        self.initialized_ = {"niter": False, "base_kernel": False}

    def initialize_(self):
        """Initialize all transformer arguments, needing initialization."""
        if not self.initialized_["base_kernel"]:
            base_kernel = self.base_kernel
            if base_kernel is not None:
                if type(base_kernel) is type and issubclass(base_kernel, kernel):
                    params = dict()
                else:
                    try:
                        base_kernel, params = base_kernel
                    except Exception:
                        raise TypeError('Base kernel was not formulated in '
                                        'the correct way. '
                                        'Check documentation.')

                    if not (type(base_kernel) is type and
                            issubclass(base_kernel, kernel)):
                        raise TypeError('The first argument must be a valid '
                                        'grakel.kernel.kernel Object')
                    if type(params) is not dict:
                        raise ValueError('If the second argument of base '
                                         'kernel exists, it must be a diction'
                                         'ary between parameters names and '
                                         'values')
                    params.pop("normalize", None)

                params["normalize"] = False
                params["verbose"] = self.verbose
                params["executor"] = self.executor
                self._base_kernel = lambda *args: base_kernel(**params)
            else:
                raise ValueError('Upon initialization base_kernel cannot be '
                                 'None')
            self.initialized_["base_kernel"] = True

        if not self.initialized_["niter"]:
            if type(self.niter) is not int or self.niter <= 0:
                raise TypeError("'niter' must be a positive integer")
            self.initialized_["niter"] = True

    def parse_input(self, X):
        """Parse input for weisfeiler lehman.

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
        base_kernel : object
            Returns base_kernel. Only if called from `fit` or `fit_transform`.

        K : np.array
            Returns the kernel matrix. Only if called from `transform` or
            `fit_transform`.

        """
        if self._base_kernel is None:
            raise ValueError('User must provide a base_kernel')
        # Input validation and parsing
        if not isinstance(X, collections.Iterable):
            raise TypeError('input must be an iterable\n')
        else:
            nx = 0
            if self._method_calling in [1, 2]:
                nl = 0
                labels_enum = dict()
                base_kernel = dict()
                for kidx in range(self.niter):
                    base_kernel[kidx] = self._base_kernel()
            elif self._method_calling == 3:
                nl = len(self._labels_enum)
                labels_enum = dict(self._labels_enum)
                base_kernel = self.X
            inp = list()
            neighbors = list()
            for (idx, x) in enumerate(iter(X)):
                is_iter = False
                if isinstance(x, collections.Iterable):
                    x, is_iter = list(x), True
                if is_iter and len(x) in [0, 2, 3]:
                    if len(x) == 0:
                        warnings.warn('Ignoring empty element on index: '
                                      + str(idx))
                        continue
                    else:
                        x = Graph(x[0], x[1], {},
                                  graph_format=self._graph_format)
                elif type(x) is not Graph:
                    raise TypeError('each element of X must be either a ' +
                                    'graph object or a list with at least ' +
                                    'a graph like object and node labels ' +
                                    'dict \n')

                label = x.get_labels(purpose='any')
                inp.append((x.get_graph_object(), label))
                neighbors.append(x.get_edge_dictionary())
                for v in set(itervalues(label)):
                    if v not in labels_enum:
                        labels_enum[v] = nl
                        nl += 1
                nx += 1
            if nx == 0:
                raise ValueError('parsed input is empty')

        # Calculate the hadamard matrix
        H = hadamard(int(2**(ceil(log2(nl)))))

        # Intial labeling of vertices based on their corresponding Hadamard
        # code (i-th row of the Hadamard matrix) where i is the i-th label on
        # enumeration
        new_graphs = list()
        for (obj, label) in inp:
            new_labels = dict()
            for (k, v) in iteritems(label):
                new_labels[k] = H[labels_enum[v], :]
            new_graphs.append((obj, label))

        # Add the zero iteration element
        if self._method_calling == 1:
            base_kernel[0].fit(new_graphs)
        elif self._method_calling == 2:
            K = base_kernel[0].fit_transform(new_graphs)
        elif self._method_calling == 3:
            K = base_kernel[0].transform(new_graphs)

        # Main
        for i in range(1, self.niter):
            graphs = new_graphs
            new_graphs = list()
            for ((obj, old_labels), neighbor) in zip(graphs, neighbors):
                # Find unique labels and sort them for both graphs
                # Keep for each node the temporary
                new_labels = dict()
                for (k, ns) in iteritems(neighbor):
                    new_labels[k] = old_labels[k]
                    for q in ns:
                        new_labels[k] = np.add(new_labels[k], old_labels[q])
                new_graphs.append((obj, new_labels))

            # calculate kernel
            if self._method_calling == 1:
                base_kernel[i].fit(new_graphs)
            elif self._method_calling == 2:
                K += base_kernel[i].fit_transform(new_graphs)
            elif self._method_calling == 3:
                K += base_kernel[i].transform(new_graphs)

        if self._method_calling == 1:
            self._labels_enum = labels_enum
            return base_kernel
        elif self._method_calling == 2:
            self._labels_enum = labels_enum
            return K, base_kernel
        elif self._method_calling == 3:
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
        check_is_fitted(self, ['X'])

        # Input validation and parsing
        if X is None:
            raise ValueError('transform input cannot be None')
        else:
            km = self.parse_input(X)

        if self.normalize:
            X_diag, Y_diag = self.diagonal()
            km /= np.sqrt(np.outer(Y_diag, X_diag))

        return km

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

        y : None
            There is no need of a target in a transformer, yet the pipeline API
            requires this parameter.

        Returns
        -------
        K : numpy array, shape = [n_targets, n_input_graphs]
            corresponding to the kernel matrix, a calculation between
            all pairs of graphs between target an features

        """
        self._method_calling = 2
        self.initialize_()
        if X is None:
            raise ValueError('transform input cannot be None')
        else:
            km, self.X = self.parse_input(X)

        self._X_diag = np.reshape(np.diagonal(km), (km.shape[0], 1))
        if self.normalize:
            return np.divide(km,
                             np.sqrt(np.outer(self._X_diag, self._X_diag)))
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
        # Check if fit had been called
        check_is_fitted(self, ['X'])
        try:
            check_is_fitted(self, ['_X_diag'])
            Y_diag = self.X[0].diagonal()[1]
            for i in range(1, self.niter):
                Y_diag += self.X[i].diagonal()[1]
        except NotFittedError:
            # Calculate diagonal of X
            X_diag, Y_diag = self.X[0].diagonal()
            X_diag.flags.writeable = True
            for i in range(1, self.niter):
                x, y = self.X[i].diagonal()
                X_diag += x
                Y_diag += y
            self._X_diag = X_diag

        return self._X_diag, Y_diag
