"""Multiscale Laplacian Graph Kernel as defined in :cite:`kondor2016multiscale`."""
# Author: Ioannis Siglidis <y.siglidis@gmail.com>
# License: BSD 3 clause
# Python 2/3 cross-compatibility import
from __future__ import print_function

import collections
import warnings
import numpy as np

from numbers import Real
from math import exp

from sklearn.utils import check_random_state

from numpy.linalg import eig
from numpy.linalg import inv
from numpy.linalg import multi_dot
from numpy.linalg import eigvals

from grakel.graph import Graph
from scipy.sparse.csgraph import laplacian

from grakel.kernels import Kernel

positive_eigenvalue_limit = float("+1e-6")


class MultiscaleLaplacian(Kernel):
    """Laplacian Graph Kernel as proposed in :cite:`kondor2016multiscale`.

    Parameters
    ----------
    random_state :  RandomState or int, default=None
        A random number generator instance or an int to initialize a RandomState as a seed.

    L : int, default=3
        The number of neighborhoods.

    gamma : Real, default=0.01
        A smoothing parameter of float value.

    heta : float, default=0.01
        A smoothing parameter of float value.

    P : int, default=10
        Restrict the maximum number of eigenvalues, taken on eigenvalue decomposition.

    n_samples : int, default=50
        The number of vertex samples.

    Attributes
    ----------
    random_state_ : RandomState
        A RandomState object handling all randomness of the class.

    _data_level : dict
        A dictionary containing the feature basis information needed
        for each level calculation on transform.

    """

    _graph_format = "adjacency"

    def __init__(self,
                 n_jobs=None,
                 normalize=False, verbose=False,
                 random_state=None,
                 L=3,
                 P=10,
                 gamma=0.01,
                 heta=0.01,
                 n_samples=50):
        """Initialise a `multiscale_laplacian` kernel."""
        super(MultiscaleLaplacian, self).__init__(
            n_jobs=n_jobs,
            normalize=normalize,
            verbose=verbose)

        self.random_state = random_state
        self.gamma = gamma
        self.heta = heta
        self.L = L
        self.P = P
        self.n_samples = n_samples
        self._initialized.update({"random_state": False, "gamma": False,
                                  "heta": False, "L": False, "n_samples": False, "P": False})

    def initialize(self):
        """Initialize all transformer arguments, needing initialization."""
        super(MultiscaleLaplacian, self).initialize()

        if not self._initialized["random_state"]:
            self.random_state_ = check_random_state(self.random_state)
            self._initialized["random_state"] = True

        if not self._initialized["gamma"]:
            if not isinstance(self.gamma, Real):
                raise TypeError('gamma must be a real number')
            elif self.gamma == .0:
                warnings.warn('with zero gamma the calculation may crash')
            elif self.gamma < 0:
                raise TypeError('gamma must be positive')
            self._initialized["gamma"] = True

        if not self._initialized["heta"]:
            if not isinstance(self.heta, Real):
                raise TypeError('heta must be a real number')
            elif self.heta == .0:
                warnings.warn('with zero heta the calculation may crash')
            elif self.heta < 0:
                raise TypeError('heta must be positive')
            self._initialized["heta"] = True

        if not self._initialized["L"]:
            if type(self.L) is not int:
                raise TypeError('L must be an integer')
            elif self.L < 0:
                raise TypeError('L must be positive')
            self._initialized["L"] = True

        if not self._initialized["n_samples"]:
            if type(self.n_samples) is not int or self.n_samples <= 0:
                raise TypeError('n_samples must be a positive integer')
            self._initialized["n_samples"] = True

        if not self._initialized["P"]:
            if type(self.P) is not int or self.P <= 0:
                raise TypeError('P must be a positive integer')
            self._initialized["P"] = True

    def parse_input(self, X):
        """Fast ML Graph Kernel.

        See supplementary material :cite:`kondor2016multiscale`, algorithm 1.

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
            A list of tuples with S matrices inverses
            and their 4th-root determinants.

        """
        if not isinstance(X, collections.Iterable):
            raise TypeError('input must be an iterable\n')
        else:
            ng = 0
            out = list()
            data = dict()
            neighborhoods = dict()
            for (idx, x) in enumerate(iter(X)):
                is_iter = False
                if isinstance(x, collections.Iterable):
                    is_iter, x = True, list(x)
                if is_iter and len(x) in [0, 2, 3]:
                    if len(x) == 0:
                        warnings.warn('Ignoring empty element ' +
                                      'on index: '+str(idx))
                        continue
                    else:
                        x = Graph(x[0], x[1], {}, self._graph_format)
                elif type(x) is Graph:
                    x.desired_format(self._graph_format)
                else:
                    raise TypeError('each element of X must be either a '
                                    'graph or an iterable with at least 1 '
                                    'and at most 3 elements\n')
                phi_d = x.get_labels()
                A = x.get_adjacency_matrix()
                try:
                    phi = np.array([list(phi_d[i]) for i in range(A.shape[0])])
                except TypeError:
                    raise TypeError('Features must be iterable and castable '
                                    'in total to a numpy array.')

                Lap = laplacian(A).astype(float)
                _increment_diagonal_(Lap, self.heta)
                data[ng] = {0: A, 1: phi, 2: inv(Lap)}
                neighborhoods[ng] = x
                ng += 1

            if ng == 0:
                raise ValueError('parsed input is empty')

            # Define a function for calculating the S's of subgraphs of each iteration
            def calculate_C(k, j, l):
                if type(neighborhoods[k]) is Graph:
                    neighborhoods[k] = neighborhoods[k].produce_neighborhoods(
                        r=self.L, sort_neighbors=False)

                indexes = neighborhoods[k][l][j]
                L = laplacian(data[k][0][indexes, :][:, indexes]).astype(float)
                _increment_diagonal_(L, self.heta)
                U = data[k][1][indexes, :]
                S = multi_dot((U.T, inv(L), U))
                _increment_diagonal_(S, self.gamma)

                return (inv(S), np.sum(np.log(np.real(eigvals(S)))))

            if self._method_calling == 1:
                V = [(k, j) for k in range(ng)
                     for j in range(data[k][0].shape[0])]

                ns = min(len(V), self.n_samples)

                self.random_state_.shuffle(V)
                vs = V[:ns]
                phi_k = np.array([data[k][1][j, :] for (k, j) in vs])

                # w the eigen vectors, v the eigenvalues
                K = phi_k.dot(phi_k.T)

                # Calculate eigenvalues
                v, w = eig(K)
                v, w = np.real(v), np.real(w.T)

                # keep only the positive
                vpos = np.argpartition(v, -self.P)[-self.P:]
                vpos = vpos[np.where(v[vpos] > positive_eigenvalue_limit)]

                # ksi.shape = (k, Ns) * (Ns, P)
                ksi = w[vpos].dot(phi_k).T / np.sqrt(v[vpos])
                for j in range(ng):
                    # (n_samples, k) * (k, P)
                    data[j][1] = data[j][1].dot(ksi)
                self._data_level = {0: ksi}
                for l in range(1, self.L+1):
                    # Take random samples from all the vertices of all graphs
                    self.random_state_.shuffle(V)
                    vs = V[:ns]

                    # Compute the reference subsampled Gram matrix
                    K_proj = {k: np.zeros(shape=(data[k][0].shape[0], ns)) for k in range(ng)}
                    K, C = np.zeros(shape=(len(vs), len(vs))), dict()
                    for (m, (k, j)) in enumerate(vs):
                        C[m] = calculate_C(k, j, l)
                        K_proj[k][j, m] = K[m, m] = self.pairwise_operation(C[m], C[m])
                        for (s, (k2, j2)) in enumerate(vs):
                            if s < m:
                                K[s, m] = K[m, s] \
                                        = K_proj[k2][j2, m] \
                                        = K_proj[k][j, s] \
                                        = self.pairwise_operation(C[s], C[m])
                            else:
                                break

                    # Compute the kernels of the relations of the reference to everything else
                    for (k, j) in V[ns:]:
                        for (m, _) in enumerate(vs):
                            K_proj[k][j, m] = self.pairwise_operation(C[m], calculate_C(k, j, l))

                    # w the eigen vectors, v the eigenvalues
                    v, w = eig(K)
                    v, w = np.real(v), np.real(w.T)

                    # keep only the positive
                    vpos = np.argpartition(v, -self.P)[-self.P:]
                    vpos = vpos[np.where(v[vpos] > positive_eigenvalue_limit)]

                    # Q shape=(k, P)
                    Q = w[vpos].T / np.sqrt(v[vpos])
                    for j in range(ng):
                        # (n, ns) * (ns, P)
                        data[j][1] = K_proj[j].dot(Q)
                    self._data_level[l] = (C, Q)

            elif self._method_calling == 3:
                ksi = self._data_level[0]
                for j in range(ng):
                    # (n, k) * (k, P)
                    data[j][1] = data[j][1].dot(ksi)

                for l in range(1, self.L+1):
                    C, Q = self._data_level[l]
                    for j in range(ng):
                        K_proj = np.zeros(shape=(data[j][0].shape[0], len(C)))
                        for n in range(data[j][0].shape[0]):
                            for m in range(len(C)):
                                K_proj[n, m] = self.pairwise_operation(C[m], calculate_C(j, n, l))
                        data[j][1] = K_proj.dot(Q)

            # Apply the final calculation of S.
            for k in range(ng):
                S = multi_dot((data[k][1].T, data[k][2], data[k][1]))
                _increment_diagonal_(S, self.gamma)
                out.append((inv(S), np.sum(np.log(np.real(eigvals(S))))))

            return out

    def pairwise_operation(self, x, y):
        """FLG calculation for the fast multiscale laplacian.

        Parameters
        ----------
        x, y : tuple
            An np.array of inverse and the log determinant of S
            (for the calculation of S matrices see the algorithm 1
             of the supplement material in cite:`kondor2016multiscale`).

        Returns
        -------
        kernel : number
            The FLG core kernel value.

        """
        S_inv_x, log_det_x = x
        S_inv_y, log_det_y = y

        # Calculate the result in term of logs
        log_detS = -np.sum(np.log(np.real(eigvals(S_inv_x + S_inv_y))))
        logr = (log_detS - 0.5*(log_det_x + log_det_y))/2.0

        if logr < -30:
            return .0
        else:
            return exp(logr)


def _increment_diagonal_(A, value):
    """Increment the diagonal of an array by a value.

    Parameters
    ----------
    A : np.array
        The array whose diagonal will be extracted.

    value : number
        The value that will be incremented on the diagonal.


    Returns
    -------
    None.

    """
    d = A.diagonal()
    d.setflags(write=True)
    d += value
