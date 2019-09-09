"""The svm theta kernel as defined in :cite:`johansson2014global`."""
# Author: Ioannis Siglidis <y.siglidis@gmail.com>
# License: BSD 3 clause
import collections
import warnings

import numpy as np

from scipy.linalg import eigvalsh
from sklearn.svm import OneClassSVM
from sklearn.utils import check_random_state

from grakel.kernels import Kernel
from grakel.graph import Graph
from grakel.tools import distribute_samples

positive_eigenvalue_limit = float("+1e-6")
min_weight = float("1e-10")


def _inner(x, y):
    return np.inner(x, y)[0, 0]


class SvmTheta(Kernel):
    """Calculate the SVM theta kernel.

    See :cite:`johansson2014global`.

    Parameters
    ----------
    X,Y : *valid-graph-format*
        The pair of graphs on which the kernel is applied.

    n_samples : int, default=50
        Number of samples.

    subsets_size_range : tuple, len=2, default=(2,8)
        (min, max) size of the vertex set of sampled subgraphs.

    metric : function (number, number -> number), default=:math:`f(x,y)=x*y`
        The applied metric between the svm_theta numbers of the two graphs.

    random_state :  RandomState or int, default=None
        A random number generator instance or an int to initialize a RandomState as a seed.

    Attributes
    ----------
    random_state_ : RandomState
        A RandomState object handling all randomness of the class.

    """

    _graph_format = "adjacency"

    def __init__(self, n_jobs=None, normalize=False,
                 verbose=False, random_state=None, n_samples=50,
                 subsets_size_range=(2, 8), metric=_inner):
        """Initialise a lovasz_theta kernel."""
        # setup valid parameters and initialise from parent
        super(SvmTheta, self).__init__(n_jobs=n_jobs,
                                       normalize=normalize,
                                       verbose=verbose)

        self.n_samples = n_samples
        self.subsets_size_range = subsets_size_range
        self.metric = metric
        self.random_state = random_state
        self._initialized.update({"n_samples": False, "subsets_size_range": False,
                                  "metric": False, "random_state": False})

    def initialize(self):
        """Initialize all transformer arguments, needing initialization."""
        super(SvmTheta, self).initialize()
        if not self._initialized["n_samples"]:
            if self.n_samples <= 0 or type(self.n_samples) is not int:
                raise TypeError('n_samples must an integer be bigger '
                                'than zero')
            self._initialized["n_samples"] = True

        if not self._initialized["subsets_size_range"]:
            if (type(self.subsets_size_range) is not tuple
                    or len(self.subsets_size_range) != 2
                    or any(type(i) is not int for i in self.subsets_size_range)
                    or self.subsets_size_range[0] > self.subsets_size_range[1]
                    or self.subsets_size_range[0] <= 0):
                raise TypeError('subsets_size_range subset size range'
                                'must be a tuple of two integers in '
                                'increasing order, bigger than 1')
            self._initialized["subsets_size_range"] = True

        if not self._initialized["metric"]:
            if not callable(self.metric):
                raise TypeError('metric between arguments' +
                                'must be a function')
            self._initialized["metric"] = True

        if not self._initialized["random_state"]:
            self.random_state_ = check_random_state(self.random_state)
            self._initialized["random_state"] = True

    def parse_input(self, X):
        """Parse and create features for svm_theta kernel.

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
        out : list
            The lovasz metrics for the given input.

        """
        if not isinstance(X, collections.Iterable):
            raise TypeError('input must be an iterable\n')
        else:
            i = 0
            out = list()
            for (idx, x) in enumerate(iter(X)):
                is_iter = False
                if isinstance(x, collections.Iterable):
                    x, is_iter = list(x), True
                if is_iter and len(x) in [0, 1, 2, 3]:
                    if len(x) == 0:
                        warnings.warn('Ignoring empty element ' +
                                      'on index: '+str(idx))
                        continue
                    else:
                        x = Graph(x[0], {}, {}, self._graph_format)
                elif type(x) is not Graph:
                    raise TypeError('each element of X must be either a ' +
                                    'graph or an iterable with at least 1 ' +
                                    'and at most 3 elements\n')
                i += 1
                A = x.get_adjacency_matrix()
                dual_coeffs = _calculate_svm_theta_(A)
                out.append(self._calculate_svm_theta_levels_(A, dual_coeffs))

            if i == 0:
                raise ValueError('parsed input is empty')

            return out

    def pairwise_operation(self, x, y):
        """Lovasz theta kernel as proposed in :cite:`johansson2014global`.

        Parameters
        ----------
        x, y : dict
            Subgraph samples metric dictionaries for all levels.

        Returns
        -------
        kernel : number
            The kernel value.

        """
        return self.metric(x, y)

    def _calculate_svm_theta_levels_(self, A, dual_coefs):
        """Calculate the svm_theta by levels for amaximum number of samples.

        Parameters
        ----------
        A : np.array, ndim=2
            The adjacency matrix.

        Returns
        -------
        level_values : np.array, shape=(num_of_levels, 1)
            Returns for all levels the mean of lovasz_numbers for
            sampled subsets.

        """
        # Calculate subsets
        n = A.shape[0]
        samples_on_subsets = distribute_samples(n, self.subsets_size_range, self.n_samples)

        # Calculate level dictionary with lovasz values
        phi = np.zeros(shape=(self.subsets_size_range[1] -
                              self.subsets_size_range[0]+1, 1))
        for (i, level) in enumerate(range(self.subsets_size_range[0],
                                          self.subsets_size_range[1]+1)):
            v = samples_on_subsets.get(level, None)
            if v is not None:
                level_values = list()
                for k in range(v):
                    if level <= n:
                        indexes = self.random_state_.choice(n, level, replace=False)
                    else:
                        indexes = range(n)
                    # calculate the metrix value for that level
                    level_values.append(np.sum(dual_coefs[indexes]))
                phi[i] = np.mean(level_values)

        return phi


def _calculate_svm_theta_(A):
    """Calculate the svm theta for the given graph.

    Parameters
    ----------
    A: np.array, ndim=2
        A square numpy array corresponding to the adjacency matrix.

    Returns
    -------
    svm_theta: float
        Returns the svm theta number.

    """
    K = (A > min_weight).astype(float)
    np.fill_diagonal(K, .0)
    min_eigv = eigvalsh(K, lower=False, eigvals=(0, 0))[0]
    if min_eigv < 0 and abs(min_eigv) > positive_eigenvalue_limit:
        K /= -min_eigv
        d = K.diagonal()
        d.setflags(write=True)
        d += 1.

    svm = OneClassSVM(kernel="precomputed")
    svm.fit(K)

    alphas = np.zeros(shape=(A.shape[0],))
    np.put(alphas, svm.support_, svm._dual_coef_[0])
    return alphas
