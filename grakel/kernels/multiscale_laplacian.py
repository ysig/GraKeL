"""Multiscale Laplacian Graph Kernel as defined in :cite:`Kondor2016TheML`."""
# Python 2/3 cross-compatibility import
from __future__ import print_function

import collections
import warnings
import numpy as np
import time

from numbers import Real
from math import sqrt

from numpy.linalg import det
from numpy.linalg import eig
from scipy.linalg import inv
from numpy.linalg import multi_dot

from grakel.graph import Graph
from scipy.sparse.csgraph import laplacian

from grakel.kernels import Kernel

# Python 2/3 cross-compatibility import
from six import iteritems

positive_eigenvalue_limit = float("+1e-6")


class MultiscaleLaplacianFast(Kernel):
    """Laplacian Graph Kernel as proposed in :cite:`Kondor2016TheML_fast`.

    Parameters
    ----------
    L : int, default=3
        The number of neighborhoods.

    gamma : Real, default=0.01
        A smoothing parameter of float value.

    heta : float, default=0.01
        A smoothing parameter of float value.

    N : int, default=50
        The number of vertex samples.

    Attributes
    ----------
    L : int
        The number of neihborhoods.

    gamma : Real
        A smoothing parameter for calculation of S matrices.

    heta : float
        A smoothing parameter for calculation of S matrices.

    N : int
        The number of vertex samples.

    ksi : np.array, len(shape)=2
        The total ksi transformation of phi's as produced on fit.

    """

    _graph_format = "adjacency"

    def __init__(self,
                 n_jobs=None,
                 normalize=False, verbose=False,
                 random_seed=42,
                 L=3,
                 gamma=0.01,
                 heta=0.01,
                 N=50):
        """Initialise a `multiscale_laplacian` kernel."""
        super(MultiscaleLaplacianFast, self).__init__(
            n_jobs=n_jobs,
            normalize=normalize,
            verbose=verbose)

        self.random_seed = random_seed
        self.gamma = gamma
        self.heta = heta
        self.L = L
        self.N = N
        self.initialized_.update({"random_seed": False, "gamma": False,
                                  "heta": False, "L": False, "N": False})

    def initialize_(self):
        """Initialize all transformer arguments, needing initialization."""
        super(MultiscaleLaplacianFast, self).initialize_()

        if not self.initialized_["random_seed"]:
            np.random.seed(self.random_seed)
            self.initialized_["random_seed"] = True

        if not self.initialized_["gamma"]:
            if not isinstance(self.gamma, Real):
                raise TypeError('gamma must be a real number')
            elif self.gamma == .0:
                warnings.warn('with zero gamma the calculation may crash')
            elif self.gamma < 0:
                raise TypeError('gamma must be positive')
            self.initialized_["gamma"] = True

        if not self.initialized_["heta"]:
            if not isinstance(self.heta, Real):
                raise TypeError('heta must be a real number')
            elif self.heta == .0:
                warnings.warn('with zero heta the calculation may crash')
            elif self.heta < 0:
                raise TypeError('heta must be positive')
            self.initialized_["heta"] = True

        if not self.initialized_["L"]:
            if type(self.L) is not int:
                raise TypeError('L must be an integer')
            elif self.L < 0:
                raise TypeError('L must be positive')
            self.initialized_["L"] = True

        if not self.initialized_["N"]:
            if type(self.N) is not int or self.N <= 0:
                raise TypeError('N must be a positive integer')
            self.initialized_["N"] = True

    def parse_input(self, X):
        """Fast ML Graph Kernel.

        See supplementary material :cite:`Kondor2016TheML`, algorithm 1.

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
                elif type(x) is not Graph:
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

            if self._method_calling == 1:
                V = [(k, j) for k in range(ng)
                     for j in range(data[k][0].shape[0])]

                ns = min(len(V), self.N)

                np.random.shuffle(V)
                vs = V[:ns]
                phi_k = np.array([data[k][1][j, :] for (k, j) in vs])

                # w the eigen vectors, v the eigenvalues
                v, w = eig(phi_k.dot(phi_k.T))
                v, w = np.real(v), np.real(w.T)

                # keep only the positive
                vpos = np.where(v > positive_eigenvalue_limit)

                # ksi shape=(k, P)
                ksi = np.dot(w[vpos], phi_k).T / np.sqrt(v[vpos])

                self._ksi = ksi
                for j in range(ng):
                    # (n, k) * (k, P)
                    data[j][1] = data[j][1].dot(ksi)

                for l in range(1, self.L+1):
                    np.random.shuffle(V)
                    vs = V[:ns]
                    K = np.zeros(shape=(len(vs), len(vs)))
                    C = dict()
                    for (m, (k, j)) in enumerate(vs):
                        if type(neighborhoods[k]) is Graph:
                            neighborhoods[k] = neighborhoods[k].produce_neighborhoods(
                                r=self.L, sort_neighbors=False)

                        indexes = neighborhoods[k][l][j]
                        L = laplacian(data[k][0][indexes, :][:, indexes]).astype(float)
                        _increment_diagonal_(L, self.heta)
                        U = data[k][1][indexes, :].T
                        S = multi_dot((U, inv(L), U.T))
                        _increment_diagonal_(S, self.gamma)
                        C[m] = (inv(S)/2.0, sqrt(sqrt(det(S))))

                    for m in range(len(C)):
                        K[m, m] = self.pairwise_operation(C[m], C[m])
                        for k in range(m + 1, len(C)):
                            K[m, k] = K[k, m] = self.pairwise_operation(C[m], C[k])

                    phi_k = np.array([data[k][1][j, :] for (k, j) in vs])

                    # w the eigen vectors, v the eigenvalues
                    v, w = eig(phi_k.dot(phi_k.T))
                    v, w = np.real(v), np.real(w.T)

                    # keep only the positive
                    vpos = np.where(v > positive_eigenvalue_limit)

                    # ksi shape=(k, P)
                    ksi = w[vpos].dot(phi_k).T / np.sqrt(v[vpos])

                    self._ksi = self._ksi.dot(ksi)
                    for j in range(ng):
                        # (n, k) * (k, P)
                        data[j][1] = data[j][1].dot(ksi)

            elif self._method_calling == 3:
                for j in range(ng):
                    # (n, k) * (k, P)
                    data[j][1] = data[j][1].dot(self._ksi)

            for k in range(ng):
                S = multi_dot((data[k][1].T, data[k][2], data[k][1]))
                _increment_diagonal_(S, self.gamma)
                out.append((inv(S)/2.0, sqrt(sqrt(det(S)))))

            return out

    def pairwise_operation(self, x, y):
        """FLG calculation for the fast multiscale laplacian.

        Parameters
        ----------
        x, y : tuple
            An np.array of inverse (divided by 2) and fourth root of the
            determinant of S matrices (for the calculation of S matrices
            see the algorithm 1 of the supplement material in
            cite:`Kondor2016TheML`).

        Returns
        -------
        kernel : number
            The FLG core kernel value.

        """
        S_inv_x, det_x_q = x
        S_inv_y, det_y_q = y

        # Calculate the kernel nominator
        k_nom = np.sqrt(det(inv(S_inv_x + S_inv_y)))

        # Calculate the kernel denominator
        k_denom = det_x_q * det_y_q

        return k_nom/k_denom


class MultiscaleLaplacian(Kernel):
    """Laplacian Graph Kernel as proposed in :cite:`Kondor2016TheML`.

    Parameters
    ----------
    L : int, default=3
        The number of neighborhoods.

    gamma : Real, default=0.01
        A small softening parameter of float value.

    heta : float, default=0.01
        A smoothing parameter of float value.

    Attributes
    ----------
    L : int
        The number of neighborhoods.

    gamma : Real
        A smoothing parameter for calculation of S matrices.

    heta : float
        A smoothing parameter for calculation of S matrices.

    """

    _graph_format = "adjacency"

    def __init__(self,
                 n_jobs=None,
                 normalize=False,
                 verbose=False,
                 L=3,
                 gamma=0.01,
                 heta=0.01):
        """Initialise a `multiscale_laplacian` kernel."""
        super(MultiscaleLaplacian, self).__init__(n_jobs=n_jobs,
                                                  normalize=normalize,
                                                  verbose=verbose)

        self.gamma = gamma
        self.heta = heta
        self.L = L
        self.initialized_.update({"gamma": False, "heta": False, "L": False})

    def initialize_(self):
        """Initialize all transformer arguments, needing initialization."""
        super(MultiscaleLaplacian, self).initialize_()

        if not self.initialized_["gamma"]:
            if not isinstance(self.gamma, Real):
                raise TypeError('gamma must be a real number')
            elif self.gamma == .0:
                warnings.warn('with zero gamma the calculation may crash')
            elif self.gamma < 0:
                raise TypeError('gamma must be a positive integer')
            self.initialized_["gamma"] = True

        if not self.initialized_["heta"]:
            if not isinstance(self.heta, Real):
                raise TypeError('heta must be a real number')

            if self.heta == .0:
                warnings.warn('with zero heta the calculation may crash')
            elif self.heta < 0:
                raise TypeError('heta must be positive')
            self.initialized_["heta"] = True

        if not self.initialized_["L"]:
            if type(self.L) is not int:
                raise TypeError('L must be an integer')
            elif self.L < 0:
                raise TypeError('L must be positive')
            self.initialized_["L"] = True

    def parse_input(self, X):
        """Parse and create features for multiscale_laplacian kernel.

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
            Tuples consisting of the Adjacency matrix, phi, phi_outer
            dictionary of neihborhood indexes and inverse laplacians
            up to level self.L and the inverse Laplacian of A.

        """
        if not isinstance(X, collections.Iterable):
            raise TypeError('input must be an iterable\n')
        else:
            ng = 0
            out = list()
            start = time.time()
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
                elif type(x) is not Graph:
                    x.desired_format(self._graph_format)
                else:
                    raise TypeError('each element of X must be either a ' +
                                    'graph or an iterable with at least 1 ' +
                                    'and at most 3 elements\n')
                ng += 1
                phi_d = x.get_labels()
                A = x.get_adjacency_matrix()
                N = x.produce_neighborhoods(r=self.L, sort_neighbors=False)
                try:
                    phi = np.array([list(phi_d[i]) for i in range(A.shape[0])])
                except TypeError:
                    raise TypeError('Features must be iterable and castable ' +
                                    'in total to a numpy array.')
                phi_outer = np.dot(phi, phi.T)

                Lap = laplacian(A).astype(float)
                _increment_diagonal_(Lap, self.heta)
                L = inv(Lap)

                Q = dict()
                for level in range(1, self.L+1):
                    Q[level] = dict()
                    for (key, item) in iteritems(N[level]):
                        Q[level][key] = dict()
                        Q[level][key]["n"] = np.array(item)
                        if len(item) < A.shape[0]:
                            laplac = laplacian(A[item, :][:, item]).astype(float)
                            _increment_diagonal_(laplac, self.heta)
                            laplac = inv(laplac)
                        else:
                            laplac = L
                        Q[level][key]["l"] = laplac

                out.append((A, phi, phi_outer, Q, L))

            if self.verbose:
                print("Preprocessing took:", time.time() - start, "s.")
            if ng == 0:
                raise ValueError('parsed input is empty')

            return out

    def pairwise_operation(self, x, y):
        """ML kernel as proposed in :cite:`Kondor2016TheML`..

        Parameters
        ----------
        x, y : tuple
            Tuple consisting of A, phi, neighborhoods up to self.L and the
            laplacian of A.

        Returns
        -------
        kernel : number
            The kernel value.

        """
        # Extract components
        Ax, phi_x, a, Qx, Lx = x
        Ay, phi_y, d, Qy, Ly = y

        nx, ny = Ax.shape[0], Ay.shape[0]

        # Create the gram matrix
        b = np.dot(phi_x, phi_y.T)
        c = b.T
        gram_matrix = np.vstack([np.hstack([a, b]), np.hstack([c, d])])

        # a lambda that calculates indexes inside the gram matrix
        # and the corresponindg laplacian given a node and a level

        for level in range(1, self.L+1):
            gram_matrix_n = np.empty(shape=gram_matrix.shape)

            for i in range(nx):
                qi = Qx[level][i]

                # xx
                for j in range(i, nx):
                    qj = Qx[level][j]
                    idx_ij = np.append(qi["n"], qj["n"])
                    extracted_gm = gram_matrix[idx_ij, :][:, idx_ij]

                    gram_matrix_n[i, j] =\
                        self._generalized_FLG_core_(qi["l"], qj["l"], extracted_gm)

                # xy
                for j in range(i, ny):
                    qj = Qy[level][j]
                    idx_ij = np.append(qi["n"], qj["n"] + nx)
                    extracted_gm = gram_matrix[idx_ij, :][:, idx_ij]

                    gram_matrix_n[i, j + nx] =\
                        self._generalized_FLG_core_(qi["l"], qj["l"], extracted_gm)

            for i in range(ny):
                idx = i + nx
                qi = Qy[level][i]
                qi_n = qi["n"] + nx

                # yx
                for j in range(i, nx):
                    qj = Qx[level][j]
                    idx_ij = np.append(qi_n, qj["n"])
                    extracted_gm = gram_matrix[idx_ij, :][:, idx_ij]

                    gram_matrix_n[idx, j] =\
                        self._generalized_FLG_core_(qi["l"], qj["l"], extracted_gm)

                # yy
                for j in range(i, ny):
                    qj = Qy[level][j]
                    idx_ij = np.append(qi_n, qj["n"] + nx)
                    extracted_gm = gram_matrix[idx_ij, :][:, idx_ij]

                    gram_matrix_n[idx, j + nx] =\
                        self._generalized_FLG_core_(qi["l"], qj["l"], extracted_gm)

            gram_matrix = np.triu(gram_matrix) + np.triu(gram_matrix, 1).T

        return self._generalized_FLG_core_(Lx, Ly, gram_matrix)

    def _generalized_FLG_core_(self, Lx, Ly, gram_matrix):
        """FLG core calculation for the multiscale gaussian.

        Parameters
        ----------
        L_{x,y} (np.array)
            Inverse laplacians of graph {x,y}.

        Gram_matrix : np.array
            The corresponding gram matrix for the two graphs.

        Returns
        -------
        kernel : number
            The FLG core kernel value.

        """
        nx = Lx.shape[0]

        # w the eigen vectors, v the eigenvalues
        v, w = eig(gram_matrix)
        v, w = np.real(v), np.real(w.T)

        # keep only the positive
        vpos = np.where(v > positive_eigenvalue_limit)[0]

        k = .0
        if vpos.shape[0] > 0:
            # calculate the Q matrix
            Q = np.square(v[vpos]) * w[vpos].T
            Qx, Qy = Q[:nx], Q[nx:]

            # Calculate the S matrices
            Sx = multi_dot((Qx.T, Lx, Qx))
            Sy = multi_dot((Qy.T, Ly, Qy))

            _increment_diagonal_(Sx, self.gamma)
            _increment_diagonal_(Sy, self.gamma)

            # A small lambda to calculate ^ 1/4
            def quatre(x):
                return sqrt(sqrt(x))

            # Caclulate the kernel nominator
            k_nom = sqrt(det(inv(inv(Sx)/2.0 + inv(Sy)/2.0)))

            # Caclulate the kernel denominator
            k_denom = quatre(det(Sx)*det(Sy))

            k = k_nom/k_denom
        return k


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
