"""The main class file representing a kernel."""
import collections
import warnings

import numpy as np

from sklearn.base import BaseEstimator
from sklearn.base import TransformerMixin
from sklearn.exceptions import NotFittedError
from sklearn.utils.validation import check_is_fitted

from grakel.graph import graph

default_verbose_value = False
default_normalize_value = False


class kernel(BaseEstimator, TransformerMixin):
    """A general class for graph kernels.

    At default a kernel is considered as pairwise. Doing so the coder that
    adds a new kernel, possibly only needs to overwrite the attributes:
    `parse_input` and `pairwise_operation` on the new kernel object.

    Parameters
    ----------
    executor : Executor, optional
        Defines

    normalize : bool, optional
        Normalize the output of the graph kernel.

    verbose : bool, optional
        Define if messages will be printed on stdout.

    Attributes
    ----------
    X : list
        Stores the input that occurs from parse input, on fit input data.
        Default format of the list objects is `grakel.graph.graph`.

    pairwise_operation : function
        A kernel between two objects as occuring with the same
        type as X (as occuring from `parse_input`).

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

    _valid_parameters = {"executor",
                         "normalize",
                         "verbose"}
    _method_calling = 0

    def __init__(self, **kargs):
        """`__init__` for `kernel` object."""
        self._verbose = kargs.get("verbose", default_verbose_value)

        self._executor = kargs.get("executor",
                                   lambda f, *args, **kargs: f(*args, **kargs))

        self._normalize = kargs.get("normalize", default_normalize_value)

        unrecognised_args = set(kargs.keys()) - self._valid_parameters

        if len(unrecognised_args) > 0:
            warnings.warn(
                'Ignoring unrecognised arguments: ' +
                ', '.join('"' + str(arg) + '"' for arg in unrecognised_args))

    def fit(self, X, y=None):
        """Fit a dataset, for a transformer.

        Parameters
        ----------
        X : iterable
            Each element must be an iterable with at most three features and at
            least one. The first that is obligatory is a valid graph structure
            (adjacency matrix or edge_dictionary) while the second is
            node_labels and the third edge_labels (that fitting the given grap
            format). The train samples.

        y : None
            There is no need of a target in a transformer, yet the pipeline API
            requires this parameter.

        Returns
        -------
        self : object
        Returns self.

        """
        self._method_calling = 1
        # Input validation and parsing
        if X is None:
            raise ValueError('fit input cannot be None')
        else:
            self.X = self.parse_input(X)

        # Return the transformer
        return self

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
        check_is_fitted(self, ['X'])

        # Input validation and parsing
        if X is None:
            raise ValueError('transform input cannot be None')
        else:
            Y = self.parse_input(X)

        # Transform - calculate kernel matrix
        km = self._calculate_kernel_matrix(Y)
        self._Y = Y
        if self._normalize:
            X_diag, Y_diag = self.diagonal()
            km /= np.sqrt(np.dot(Y_diag, X_diag.T))
        return km

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

        # Transform - calculate kernel matrix
        km = self._calculate_kernel_matrix()

        self._X_diag = np.diagonal(km).reshape(km.shape[0], 1)
        if self._normalize:
            return np.divide(km,
                             np.sqrt(np.multiply(self._X_diag.T,
                                                 self._X_diag)))
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
                K[i, i] = self._executor(self.pairwise_operation, x, x)
                for (j, y) in enumerate(cache):
                    K[j, i] = self._executor(self.pairwise_operation, y, x)
                cache.append(x)

            K = np.triu(K) + np.triu(K, 1).T

        else:
            K = np.zeros(shape=(len(Y), len(self.X)))
            for (j, y) in enumerate(Y):
                for (i, x) in enumerate(self.X):
                    K[j, i] = self._executor(self.pairwise_operation, x, y)

        return K

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
            self._X_diag = np.empty(shape=(len(self.X), 1))
            for (i, x) in enumerate(self.X):
                self._X_diag[i] = self._executor(self.pairwise_operation, x, x)

        Y_diag = np.empty(shape=(len(self._Y), 1))
        for (i, y) in enumerate(self._Y):
            Y_diag[i] = self._executor(self.pairwise_operation, y, y)

        return self._X_diag, Y_diag

    def parse_input(self, X):
        """Parse the given input and raise errors if it is invalid.

        Parameters
        ----------
        X : object
            For the input to pass the test, we must have:
            Each element must be an iterable with at most three features and at
            least one. The first that is obligatory is a valid graph structure
            (adjacency matrix or edge_dictionary) while the second is
            node_labels and the third edge_labels (that correspond to the given
            graph format). A valid input also consists of graph type objects.
            If None the kernel matrix is calculated upon fit data.
            The test samples.

        Returns
        -------
        Xp : list
            List of graph type objects.

        """
        if not isinstance(X, collections.Iterable):
            raise ValueError('input must be an iterable\n')
        else:
            Xp = list()
            for (i, x) in enumerate(iter(X)):
                if len(x) == 0:
                    warnings.warn('Ignoring empty element on index: '+str(i))
                if len(x) == 1:
                    if type(x) is graph:
                        Xp.append(x)
                    else:
                        Xp.append(graph(x[0], {}, {},
                                        self._graph_format))
                elif len(x) == 2:
                    Xp.append(graph(x[0], x[1], {}, self._graph_format))
                elif len(x) == 3:
                    Xp.append(graph(x[0], x[1], x[2], self._graph_format))
                else:
                    raise ValueError('each element of X must have at least' +
                                     ' one and at most 3 elements\n')
            if len(Xp) == 0:
                raise ValueError('parsed input is empty')
            return Xp

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
        return 0
