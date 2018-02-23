"""The lovasz theta kernel as defined in :cite:`Johansson2015LearningWS`."""
import collections
import warnings

import random
import numpy as np

from grakel.kernels import kernel
from grakel.graph import Graph
from grakel.tools import distribute_samples

from math import sqrt
from numpy import pad
from numpy.linalg import LinAlgError
from numpy.linalg import norm
from scipy.linalg import cholesky
from scipy.linalg import eigvalsh
from scipy.linalg import solve
from cvxopt.base import matrix
from cvxopt.base import spmatrix
from cvxopt.solvers import sdp
from cvxopt.solvers import options

options['show_progress'] = False
min_weight = float("1e-10")
angle_precision = float("1e-6")
tolerance = float("1e-1")


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

    base_kernel : function (np.1darray, np.1darray -> number),
                  default=:math:`f(x,y) = x*y`
        The applied metric between the lovasz_theta numbers of subgraphs.

    Attributes
    ----------
    _n_samples : int
        Number of samples drawn for the computation of lovasz theta.

    _ssr : tuple, len=2
        A tuple containing two integers designating the minimum and the maximum
        size of the vertex set of considered subgraphs.

    _d : int,
        The maximum matrix dimension of fit plus 1. Signifies the number
        of features assigned for lovasz labelling.

    _base_kernel : function (number, number -> number)
        The applied base_kernel between features of the mean lovasz_theta
        numbers samples for all levels of the two graphs.

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

        self._base_kernel = kargs.get("base_kernel", lambda x, y: x.T.dot(y))
        if not callable(self._base_kernel):
            raise ValueError('base_kernel between arguments ' +
                             'must be a function')

        rs = kargs.get("random_seed", 6578909)
        random.seed(rs)
        np.random.seed(rs)

    def parse_input(self, X):
        """Parse and create features for lovasz_theta kernel.

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
            raise ValueError('input must be an iterable\n')
        else:
            i = 0
            adjm = list()
            max_dim = 0
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
                A = x.get_adjacency_matrix()
                adjm.append(A)
                max_dim = max(max_dim, A.shape[0])

            if self._method_calling == 1:
                self._d = max_dim + 1

            out = list()
            for A in adjm:
                X, t = _calculate_lovasz_embeddings_(A)
                U = _calculate_lovasz_labelling_(X, t, self._d)
                out.append(self._calculate_MEC_(U))
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
        return self._base_kernel(x, y)

    def _calculate_MEC_(self, U):
        """Calculate the minimum enclosing cone for given U.

        Calculates the minimum a graph metric as the lovasz theta kernel or the
        svm-theta kernel, on a set of randomly sampled subgraphs, producing a
        dictionary of subgraph levels and sets.

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
        n = U.shape[1]
        samples_on_subsets = distribute_samples(n, self._ssr, self._n_samples)

        # Calculate level dictionary with lovasz values
        phi = np.zeros(shape=(self._ssr[1]-self._ssr[0]+1, 1))
        for (i, level) in enumerate(range(self._ssr[0], self._ssr[1]+1)):
            v = samples_on_subsets.get(level, None)
            if v is not None:
                level_values = list()
                for k in range(v):
                    indexes = np.random.choice(n, level, replace=False)
                    # calculate the metrix value for that level
                    level_values.append(_minimum_cone_(U[:, indexes]))
                phi[i] = np.mean(level_values)

        return phi


def _calculate_lovasz_embeddings_(A):
    """Calculate the lovasz embeddings for the given graph.

    Parameters
    ----------
    A : np.array, ndim=2
        The adjacency matrix.

    Returns
    -------
    lovasz_theta : float
        Returns the lovasz theta number.

    optimization_solution : np.array
        Returns the slack variable solution
        of the primal optimization program.

    """
    if A.shape[0] == 1:
        return 1.0
    else:
        tf_adjm = (np.abs(A) <= min_weight)
        np.fill_diagonal(tf_adjm, False)
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

    return np.array(sol['ss'][0]), sol['x'][ne]


def _calculate_lovasz_labelling_(X, t, d):
    """Calculate the lovasz labelling for the extracted embeddings.

    Parameters
    ----------
    X : np.array, ndim=2
        The optimization solution of the sdp program.

    t : float
        Returns the lovasz theta number.

    d : int
        The number of features size.

    Returns
    -------
    lovasz_theta : float
        Returns the lovasz theta number.

    optimization_solution : np.array
        Returns the slack variable solution
        of the primal optimization program.

    """
    n = X.shape[0]

    try:
        V = cholesky(X)
    except LinAlgError:
        x = X.diagonal()
        x.setflags(write=True)
        x += 2*abs(eigvalsh(X, lower=False, eigvals=(0, 0))[0])
        V = cholesky(X)

    V = pad(V, [(0, d-n), (0, 0)], mode='constant', constant_values=0)

    c = np.zeros(shape=(d,))
    c[-1] = 1

    C = np.outer(c, np.ones(shape=(n,)))

    U = 1/sqrt(t)*(C+V)
    return U


def _minimum_cone_(U):
    """Calculate the minimum cone.

    Computes the angle and center of the minimum cone, with its point
    in the origin, enclosing the vectors in A.

    Parameters
    ----------
    U : np.array, ndim=2
        The vectors that will be enclosed.

    Returns
    -------
    angle_cosine : float
        Returns the cosine of the minimum angle.

    """
    n = U.shape[1]

    P = np.random.permutation(n) - 1
    R = np.array([], dtype=int)
    c, _ = _b_minidisk_(U, P, R)

    with np.errstate(divide='ignore'):
        c /= norm(c, 2)

    t = min(np.dot(U.T, c))
    if t > 1. and t < angle_precision:
        t = 1.

    elif t < -1. and t > -angle_precision:
        t = -1.

    return t


def _b_minidisk_(A, P, R):
    """Calculate the minidisk.

    Implements Welzl's algorithm (*move-to-front*)

    Parameters
    ----------
    A : np.array, ndim=2
        The vectors that will be enclosed.

    P : np.array, ndim=1
        A random permutation of indeces.

    R : np.array, ndim=1
        A subset of vectors that are enclosed.

    Returns
    -------
    c : np.array, ndim=1
        The center vector C of the cone as it is being minimized.

    r : int
        The cone radius as it is beeing minimized.

    """
    d, nP, nR = A.shape[0], P.shape[0], R.shape[0]

    # original algorithm
    if nP == 0 or nR == d+1:
        if nR == 0:
            c, r = np.zeros(shape=(d,)), 0
        else:
            c, r = _fitball_(A[:, R])
    else:
        p = P[random.randint(0, nP - 1)]
        P_prime = np.delete(P, np.where(P == p))
        c, r = _b_minidisk_(A, P_prime, R)
        if norm(A[:, p] - c, 2) - r > tolerance:
            # if not inside ball
            if p not in R:
                R_prime = pad(R, [(0, 1)], mode='constant', constant_values=p)
                c, r = _b_minidisk_(A, P_prime, R_prime)
    return c, r


def _fitball_(A):
    """Fit the minimum ball.

    Parameters
    ----------
    A : np.array, ndim=2
        The vectors that will be enclosed inside the ball.

    Returns
    -------
    c : np.array, ndim=1
        The center vector C of the ball.

    r : int
        The ball radius.

    """
    d = A.shape[0]
    n = A.shape[1]

    if n == 1:
        c, r = A[:, 0], 0
    else:
        Q = A-np.outer(A[:, 0], np.ones(shape=(n, 1)))
        B = 2*np.dot(Q.T, Q)
        b = B.diagonal()/2

        L = solve(B[1:, :][:, 1:], b[1:])
        L = pad(L, [(1, 0)], mode='constant', constant_values=0)

        C = np.zeros(shape=(d,))

        for i in range(1, n):
            C = C + L[i]*Q[:, i]

        r = np.sqrt(np.dot(C, C))
        c = C + A[:, 1]

    return c, r
