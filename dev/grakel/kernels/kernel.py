"""The main class file representing a kernel."""
import collections
import warnings
import copy

import numpy as np

from sklearn.base import BaseEstimator
from sklearn.base import TransformerMixin
from sklearn.exceptions import NotFittedError
from sklearn.utils.validation import check_is_fitted

from grakel.graph import Graph

# Python 2/3 cross-compatibility import
from six import iteritems

default_executor = lambda fn, *eargs, **ekargs: fn(*eargs, **ekargs)


class Kernel(BaseEstimator, TransformerMixin):
    """A general class for graph kernels.

    At default a kernel is considered as pairwise. Doing so the coder that
    adds a new kernel, possibly only needs to overwrite the attributes:
    `parse_input` and `pairwise_operation` on the new kernel object.

    Parameters
    ----------
    executor : Executor, optional
        Defines a general executor that can be applied for concurrent
        computations.

    normalize : bool, optional
        Normalize the output of the graph kernel.

    verbose : bool, optional
        Define if messages will be printed on stdout.

    Attributes
    ----------
    X : list
        Stores the input that occurs from parse input, on fit input data.
        Default format of the list objects is `grakel.graph.graph`.

    _graph_format : str
        Stores in which type the graphs will need to be stored.

    _verbose : bool
        Defines if two print arguments on stdout.

    _executor : Executor
        A general executor that can be applied for concurrent computations.

    _normalize : bool
        Defines if normalization will be applied on the kernel matrix.

    _valid_parameters : set
        Holds the default valid parameters names for initialization.

    _method_calling : int
        An inside enumeration defines which method calls another method.
            - 1 stands for fit
            - 2 stands for fit_transform
            - 3 stands for transform

    """

    X = None
    _graph_format = "dictionary"
    _method_calling = 0

    def __init__(self,
                 executor=default_executor,
                 normalize=False,
                 verbose=False):
        """`__init__` for `kernel` object."""
        self.verbose = verbose
        self.executor = default_executor
        self.normalize = normalize
        self.initialized_ = dict()

    def fit(self, X, y=None):
        """Fit a dataset, for a transformer.

        Parameters
        ----------
        X : iterable
            Each element must be an iterable with at most three features and at
            least one. The first that is obligatory is a valid graph structure
            (adjacency matrix or edge_dictionary) while the second is
            node_labels and the third edge_labels (that fitting the given graph
            format). The train samples.

        y : None
            There is no need of a target in a transformer, yet the pipeline API
            requires this parameter.

        Returns
        -------
        self : object
        Returns self.

        """
        self._is_transformed = False
        self._method_calling = 1

        # Parameter initialization
        self.initialize_()

        # Input validation and parsing
        if X is None:
            raise ValueError('`fit` input cannot be None')
        else:
            self.X = self.parse_input(X)

        # Return the transformer
        return self

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
            raise ValueError('`transform` input cannot be None')
        else:
            Y = self.parse_input(X)

        # Transform - calculate kernel matrix
        km = self._calculate_kernel_matrix(Y)
        self._Y = Y

        # Self transform must appear before the diagonal call on normilization
        self._is_transformed = True
        if self.normalize:
            X_diag, Y_diag = self.diagonal()
            km /= np.sqrt(np.outer(Y_diag, X_diag))
        return km

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
        self.fit(X)

        # Transform - calculate kernel matrix
        km = self._calculate_kernel_matrix()

        self._X_diag = np.diagonal(km)
        if self.normalize:
            return km / np.sqrt(np.outer(self._X_diag, self._X_diag))
        else:
            return km

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
        if Y is None:
            K = np.zeros(shape=(len(self.X), len(self.X)))
            cache = list()
            for (i, x) in enumerate(self.X):
                K[i, i] = self.executor(self.pairwise_operation, x, x)
                for (j, y) in enumerate(cache):
                    K[j, i] = self.executor(self.pairwise_operation, y, x)
                cache.append(x)
            K = np.triu(K) + np.triu(K, 1).T

        else:
            K = np.zeros(shape=(len(Y), len(self.X)))
            for (j, y) in enumerate(Y):
                for (i, x) in enumerate(self.X):
                    K[j, i] = self.executor(self.pairwise_operation, y, x)
        return K

    def diagonal(self):
        """Calculate the kernel matrix diagonal of the fit/transformed data.

        Parameters
        ----------
        None.

        Returns
        -------
        X_diag : np.array
            The diagonal of the kernel matrix between the fitted data.
            This consists of each element calculated with itself.

        Y_diag : np.array
            The diagonal of the kernel matrix, of the transform.
            This consists of each element calculated with itself.

        """
        # Check is fit had been called
        check_is_fitted(self, ['X'])
        try:
            check_is_fitted(self, ['_X_diag'])
        except NotFittedError:
            # Calculate diagonal of X
            self._X_diag = np.empty(shape=(len(self.X),))
            for (i, x) in enumerate(self.X):
                self._X_diag[i] = self.executor(self.pairwise_operation, x, x)

        try:
            # If transform has happened return both diagonals
            check_is_fitted(self, ['_Y'])
            Y_diag = np.empty(shape=(len(self._Y),))
            for (i, y) in enumerate(self._Y):
                Y_diag[i] = self.executor(self.pairwise_operation, y, y)

            return self._X_diag, Y_diag
        except NotFittedError:
            # Else just return both X_diag
            return self._X_diag

    def parse_input(self, X):
        """Parse the given input and raise errors if it is invalid.

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
        Xp : list
            List of graph type objects.

        """
        if not isinstance(X, collections.Iterable):
            raise TypeError('input must be an iterable\n')
        else:
            Xp = list()
            for (i, x) in enumerate(iter(X)):
                is_iter = isinstance(x, collections.Iterable)
                if is_iter:
                    x = list(x)
                if is_iter and len(x) in [0, 1, 2, 3]:
                    if len(x) == 0:
                        warnings.warn('Ignoring empty element' +
                                      'on index: '+str(i)+'..')
                        continue
                    elif len(x) == 1:
                        Xp.append(Graph(x[0], {}, {},
                                        self._graph_format))
                    elif len(x) == 2:
                        Xp.append(Graph(x[0], x[1], {}, self._graph_format))
                    else:
                        Xp.append(Graph(x[0], x[1], x[2], self._graph_format))
                elif type(x) is Graph:
                    Xp.append(x)
                else:
                    raise TypeError('Each element of X must have at least ' +
                                    'one and at most 3 elements.\n')
            if len(Xp) == 0:
                raise ValueError('Parsed input is empty.')
            return Xp

    def initialize_(self):
        """Initialize all transformer arguments, needing initialisation."""
        pass

    def pairwise_operation(self, x, y):
        """Calculate a pairwise kernel between two elements.

        Parameters
        ----------
        x, y : Object
            Objects as occur from parse_input.

        Returns
        -------
        kernel : number
            The kernel value.

        """
        raise NotImplementedError('Pairwise operation is not implemented!')

    def set_params(self, **params):
        """Call the parent method."""
        if len(self.initialized_):
            # Copy the parameters
            params = copy.deepcopy(params)

            # Iterate over the parameters
            for key, value in iteritems(params):
                key, delim, sub_key = key.partition('__')
                if delim:
                    if sub_key in self.initialized_:
                        self.initialized_[sub_key] = False
                elif key in self.initialized_:
                    self.initialized_[key] = False

        # Set parameters
        super(Kernel, self).set_params(**params)
