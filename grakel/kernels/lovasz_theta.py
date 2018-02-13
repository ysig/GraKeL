"""The lovasz theta kernel as defined in :cite:`Johansson2015LearningWS`."""
import itertools
import collections
import warnings

import numpy as np

from grakel.kernels import kernel
from grakel.graph import Graph
from grakel.tools import distribute_samples

from cvxopt.base import matrix
from cvxopt.base import spmatrix
from cvxopt.solvers import sdp
from cvxopt.solvers import options

options['show_progress'] = False


class lovasz_theta(kernel):
    """Lovasz theta kernel as proposed in :cite:`Johansson2015LearningWS`.

    Parameters
    ----------
    X,Y : *valid-graph-format*
        The pair of graphs on which the kernel is applied.

    n_samples : int, default=50
        The number of samples.

    subsets_size_range : tuple, len=2, default=(2,8)
        (min, max) size of the vertex set of sampled subgraphs.

    metric : function (number, number -> number), default=:math:`f(x,y) = x*y`
        The applied metric between the lovasz_theta numbers of the two graphs.

    Attributes
    ----------
    _n_samples : int
        Number of samples drawn for the computation of lovasz theta.

    _ssr : tuple, len=2
        A tuple containing two integers designating the minimum and the maximum
        size of the vertex set of considered subgraphs.

    _metric : function (number, number -> number)
        The applied metric between the lovasz_theta numbers of the two graphs.

    """

    _graph_format = "adjacency"

    def __init__(self, **kargs):
        """Initialise a lovasz_theta kernel."""
        # setup valid parameters and initialise from parent
        self._valid_parameters |= {"n_samples", "random_seed",
                                   "subsets_size_range", "metric"}
        super(lovasz_theta, self).__init__(**kargs)

        self._n_samples = kargs.get("n_samples", 50)
        if self._n_samples <= 0 or type(self._n_samples) is not int:
            raise ValueError('n_samples must an integer be bigger than zero')

        self._ssr = kargs.get("subsets_size_range", (2, 8))
        if (type(self._ssr) is not tuple or len(self._ssr) != 2 or
                any(type(i) is not int for i in self._ssr) or
                self._ssr[0] > self._ssr[1]):
            raise ValueError('subsets_size_range subset size range must ' +
                             'be a tuple of two integers in increasing order')

        self._metric = kargs.get("metric", lambda x, y: x*y)
        if not callable(self._metric):
            raise ValueError('metric between arguments must be a function')

        np.random.seed(kargs.get("random_seed", 6578909))

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
        out : list
            The lovasz metrics for the given input.

        """
        if not isinstance(X, collections.Iterable):
            raise ValueError('input must be an iterable\n')
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
                    raise ValueError('each element of X must be either a ' +
                                     'graph or an iterable with at least 1 ' +
                                     'and at most 3 elements\n')
                i += 1
                out.append(
                    self._calculate_subgraph_samples_metric_dictionary_(
                        x.get_adjacency_matrix()))

            if i == 0:
                raise ValueError('parsed input is empty')

            return out

    def pairwise_operation(self, x, y):
        """Lovasz theta kernel as proposed in :cite:`Johansson2015LearningWS`.

        Parameters
        ----------
        x, y : dict
            Subgraph samples metric dictionaries for all levels.

        Returns
        -------
        kernel : number
            The kernel value.

        """
        kernel = 0

        for level in x.keys():
            if level in y and bool(x[level]) and bool(y[level]):
                    Z = len(x[level])*len(y[level])
                    kernel += sum(self._metric(a, b) for (a, b) in
                                  itertools.product(x[level], y[level]))/Z

        return kernel

    def _calculate_subgraph_samples_metric_dictionary_(self, A):
        """Calculate a graph metric.

        Calculates a graph metric as the lovasz theta kernel or the svm-theta
        kernel, on a set of randomly sampled subgraphs, producing a dictionary
        of subgraph levels and sets.

        Parameters
        ----------
        A : np.array, ndim=2
            The adjacency matrix.

        Returns
        -------
        level_values : dict
            Returns a dictionary with levels (subsite size) and the lovasz
            value of all the sampled subgraphs.

        """
        # Calculate subsets
        n = A.shape[0]
        samples_on_subsets = distribute_samples(n, self._ssr, self._n_samples)

        # Calculate level dictionary with lovasz values
        level_values = collections.defaultdict(list)
        for (level, v) in samples_on_subsets.items():
            subsets = set()
            for k in range(v):
                while True:
                    indexes = np.random.choice(n, level, replace=False)
                    tup_indexes = tuple(indexes)
                    if tup_indexes not in subsets:
                        subsets.add(tup_indexes)
                        break
                # calculate the metrix value for that level
                level_values[level].append(
                    self._internal_metric(A[:, indexes][indexes, :]))

        return level_values

    def _internal_metric(self, A):
        """Calculate the internal metric.

        Parameters
        ----------
        A : np.array, ndim=2
            The adjacency matrix.

        Returns
        -------
        metric: Number
            Returns the metric number.

        """
        return _calculate_lovasz_theta_(A)


def _calculate_lovasz_theta_(A):
    """Calculate the lovasz theta for the given graph.

    Parameters
    ----------
    A : np.array, ndim=2
        The adjacency matrix.

    Returns
    -------
    lovasz_theta: float
        Returns the lovasz theta number.

    """
    if A.shape[0] == 1:
        return 1.0
    else:
        tf_adjm = A > 0
        i_list, j_list = np.nonzero(np.triu(tf_adjm.astype(int), k=1))

        nv = A.shape[0]
        ne = len(i_list)

        x_list = list()
        e_list = list()
        for (e, (i, j)) in enumerate(zip(i_list, j_list)):
            e_list.append(int(e)), x_list.append(int(i*nv+j))
            if tf_adjm[i, j]:
                e_list.append(int(e)), x_list.append(int(i*nv+j))

    # Add on the last row, diagonal elements
    e_list = e_list+(nv*[ne])
    x_list = x_list+[int(i*nv + i) for i in range(nv)]

    # initialise g sparse (to values -1, based on two list that
    # define index and one that defines shape
    g_sparse = spmatrix(-1, x_list, e_list, (nv*nv, ne+1))

    # Initialise optimization parameters
    h = matrix(-1.0, (nv, nv))
    c = matrix([0.0]*ne + [1.0])

    # Solve the convex optimization problem
    sol = sdp(c, Gs=[g_sparse], hs=[h])

    return sol['x'][ne]
