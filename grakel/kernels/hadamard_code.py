"""The hadamard code kernel as defined in :cite:`icpram16`."""
import collections
import warnings

import numpy as np

from math import ceil
from numpy import log2

from scipy.linalg import hadamard

from sklearn.exceptions import NotFittedError
from sklearn.utils.validation import check_is_fitted
from sklearn.externals import joblib

from grakel.graph import Graph
from grakel.kernels import Kernel

# Python 2/3 cross-compatibility import
from six import iteritems
from six import itervalues


class HadamardCode(Kernel):
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

    def __init__(self, n_jobs=None, verbose=False,
                 normalize=False, niter=5, base_kernel=None):
        """Initialise a `hadamard_code` kernel."""
        super(HadamardCode, self).__init__(
            n_jobs=n_jobs, verbose=verbose, normalize=normalize)

        self.niter = niter
        self.base_kernel = base_kernel
        self.initialized_.update({"niter": False, "base_kernel": False})

    def initialize_(self):
        """Initialize all transformer arguments, needing initialization."""
        super(HadamardCode, self).initialize_()

        if not self.initialized_["base_kernel"]:
            base_kernel = self.base_kernel
            if base_kernel is not None:
                if type(base_kernel) is type and issubclass(base_kernel, Kernel):
                    params = dict()
                else:
                    try:
                        base_kernel, params = base_kernel
                    except Exception:
                        raise TypeError('Base kernel was not formulated in '
                                        'the correct way. '
                                        'Check documentation.')

                    if not (type(base_kernel) is type and
                            issubclass(base_kernel, Kernel)):
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
                params["n_jobs"] = None
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
        """Parse input and create features, while initializing and/or calculating sub-kernels.

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
            nx, labels = 0, list()
            if self._method_calling in [1, 2]:
                nl, labels_enum, base_kernel = 0, dict(), dict()
                for kidx in range(self.niter):
                    base_kernel[kidx] = self._base_kernel()
            elif self._method_calling == 3:
                nl, labels_enum, base_kernel = len(self._labels_enum), dict(self._labels_enum), self.X
            inp = list()
            neighbors = list()
            for (idx, x) in enumerate(iter(X)):
                is_iter = False
                if isinstance(x, collections.Iterable):
                    x, is_iter = list(x), True
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
                            extra = (x.get_labels(purpose='any',
                                                  label_type="edge", return_none=True), ) + extra
                        else:
                            x = Graph(x[0], x[1], {}, graph_format=self._graph_format)
                            extra = tuple()
                elif type(x) is Graph:
                    el = x.get_labels(purpose=self._graph_format, label_type="edge", return_none=True)
                    if el is None:
                        extra = tuple()
                    else:
                        extra = (el, )
                else:
                    raise TypeError('each element of X must be either a ' +
                                    'graph object or a list with at least ' +
                                    'a graph like object and node labels ' +
                                    'dict \n')

                label = x.get_labels(purpose='any')
                inp.append((x.get_graph_object(), extra))
                neighbors.append(x.get_edge_dictionary())
                labels.append(label)
                for v in set(itervalues(label)):
                    if v not in labels_enum:
                        labels_enum[v] = nl
                        nl += 1
                nx += 1
            if nx == 0:
                raise ValueError('parsed input is empty')

        # Calculate the hadamard matrix
        H = hadamard(int(2**(ceil(log2(nl)))))

        def generate_graphs(labels):
            # Intial labeling of vertices based on their corresponding Hadamard code (i-th row of the
            # Hadamard matrix) where i is the i-th label on enumeration
            new_graphs, new_labels = list(), list()
            for ((obj, extra), label) in zip(inp, labels):
                new_label = dict()
                for (k, v) in iteritems(label):
                    new_label[k] = H[labels_enum[v], :]
                new_graphs.append((obj, {i: tuple(j) for (i, j) in iteritems(new_label)}) + extra)
                new_labels.append(new_label)

            yield new_graphs
            # Main
            for i in range(1, self.niter):
                new_graphs, labels, new_labels = list(), new_labels, list()
                for ((obj, extra), neighbor, old_label) in zip(inp, neighbors, labels):
                    # Find unique labels and sort them for both graphs and keep for each node
                    # the temporary
                    new_label = dict()
                    for (k, ns) in iteritems(neighbor):
                        new_label[k] = old_label[k]
                        for q in ns:
                            new_label[k] = np.add(new_label[k], old_label[q])
                    new_labels.append(new_label)
                    new_graphs.append((obj, {i: tuple(j) for (i, j) in iteritems(new_label)}) +
                                      extra)
                yield new_graphs

        if self._method_calling in [1, 2]:
            base_kernel = {i: self._base_kernel() for i in range(self.niter)}

        if self._parallel is None:
            # Add the zero iteration element
            if self._method_calling == 1:
                for (i, g) in enumerate(generate_graphs(labels)):
                    base_kernel[i].fit(g)
            elif self._method_calling == 2:
                K = np.sum((base_kernel[i].fit_transform(g) for (i, g)
                           in enumerate(generate_graphs(labels))), axis=0)
            elif self._method_calling == 3:
                # Calculate the kernel matrix without parallelization
                K = np.sum((self.X[i].transform(g) for (i, g)
                           in enumerate(generate_graphs(labels))), axis=0)

        else:
            if self._method_calling == 1:
                self._parallel(joblib.delayed(efit)(base_kernel[i], g)
                               for (i, g) in enumerate(generate_graphs(labels)))
            elif self._method_calling == 2:
                # Calculate the kernel marix with parallelization
                K = np.sum(self._parallel(joblib.delayed(efit_transform)(base_kernel[i], g) for (i, g)
                           in enumerate(generate_graphs(labels))), axis=0)
            elif self._method_calling == 3:
                # Calculate the kernel marix with parallelization
                K = np.sum(self._parallel(joblib.delayed(etransform)(self.X[i], g) for (i, g)
                           in enumerate(generate_graphs(labels))), axis=0)

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

        self._is_transformed = True
        if self.normalize:
            X_diag, Y_diag = self.diagonal()
            old_settings = np.seterr(divide='ignore')
            km /= np.sqrt(np.outer(Y_diag, X_diag))
            km = np.nan_to_num(km)
            np.seterr(**old_settings)

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
        self._is_transformed = False
        self.initialize_()
        if X is None:
            raise ValueError('transform input cannot be None')
        else:
            km, self.X = self.parse_input(X)

        self._X_diag = np.diagonal(km)
        if self.normalize:
            old_settings = np.seterr(divide='ignore')
            km = np.nan_to_num(np.divide(km, np.sqrt(np.outer(self._X_diag, self._X_diag))))
            np.seterr(**old_settings)
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
            if self._is_transformed:
                Y_diag = self.X[0].diagonal()[1]
                for i in range(1, self.niter):
                    Y_diag += self.X[i].diagonal()[1]
        except NotFittedError:
            # Calculate diagonal of X
            if self._is_transformed:
                X_diag, Y_diag = self.X[0].diagonal()
                # X_diag is considered a mutable and should not affect the kernel matrix itself.
                X_diag.flags.writeable = True
                for i in range(1, self.niter):
                    x, y = self.X[i].diagonal()
                    X_diag += x
                    Y_diag += y
                self._X_diag = X_diag
            else:
                # case sub kernel is only fitted
                X_diag = self.X[0].diagonal()
                # X_diag is considered a mutable and should not affect the kernel matrix itself.
                X_diag.flags.writeable = True
                for i in range(1, self.niter):
                    x = self.X[i].diagonal()
                    X_diag += x
                self._X_diag = X_diag

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
