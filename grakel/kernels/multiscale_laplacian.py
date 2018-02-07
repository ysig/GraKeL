"""Multiscale Laplacian Graph Kernel as defined in :cite:`Kondor2016TheML`."""
import collections
import warnings

import numpy as np

from numpy.linalg import det
from numpy.linalg import eig
from numpy.linalg import inv
from numpy.linalg import multi_dot

from grakel.graph import Graph
from grakel.graph import laplacian

from grakel.kernels import kernel

heta = 0.01


class multiscale_laplacian(kernel):
    """Laplacian Graph Kernel as proposed in :cite:`Kondor2016TheML`.

    Parameters
    ----------
    L : int, default=3
        Number of neighborhoods.

    gamma : float, default=0.01
        A small softening parameter of float value.

    Attributes
    ----------
    kernel : number
        The kernel value.

    """

    def __init__(self, **kargs):
        """Initialise a `multiscale_laplacian` kernel."""
        self._valid_parameters |= {"L", "gamma"}
        super(multiscale_laplacian, self).__init__(**kargs)

        self._gamma = kargs.get("gamma", 0.01)
        if self._gamma > 0.05:
            warnings.warn('gamma to big')
        elif self._gamma == .0:
            warnings.warn('with zero gamma the calculation may crash')
        elif self._gamma < 0:
            raise ValueError('gamma must be positive')

        self._L = kargs.get("L", 3)
        if type(self._L) is not int:
            raise ValueError('L must be an integer')
        elif self._L < 0:
            raise ValueError('L must be positive')

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
                if isinstance(x, collections.Iterable) and \
                        len(x) in [0, 2, 3]:
                    if len(x) == 0:
                        warnings.warn('Ignoring empty element ' +
                                      'on index: '+str(idx))
                        continue
                    else:
                        x = Graph(x[0], x[1], {}, self._graph_format)
                elif type(x) is not Graph:
                    x.desired_format(self._graph_format)
                else:
                    raise ValueError('each element of X must be either a ' +
                                     'graph or an iterable with at least 1 ' +
                                     'and at most 3 elements\n')
                i += 1
                phi = x.get_labels()
                A = x.get_adjacency_matrix()
                N = x.produce_neighborhoods(r=self._L)
                Lap = x.laplacian()
                out.append((A, phi, N, Lap))

            if i == 0:
                raise ValueError('parsed input is empty')

            return out

    def pairwise_operation(self, x, y):
        """Lovasz theta kernel as proposed in :cite:`Johansson2015LearningWS`.

        Parameters
        ----------
        x, y : tuple
            Tuple consisting of A, phi, neighborhoods up to self._L and the
            laplacian of A.

        Returns
        -------
        kernel : number
            The kernel value.

        """
        # Extract components
        Ax, phi_x, Nx, Lx = x
        Ay, phi_y, Ny, Ly = y

        # create the gram matrix
        nx, ny = Ax.shape[0], Ay.shape[0]
        gram_matrix_size = nx + ny
        gram_matrix = np.empty(shape=(gram_matrix_size, gram_matrix_size))

        # initialise the gram matrix
        pick = lambda i: list(phi_x[i]) if i < nx else list(phi_y[i-nx])
        for i in range(0, gram_matrix_size):
            vec_a = pick(i)
            for j in range(0, gram_matrix_size):
                vec_b = pick(i)
                gram_matrix[i, j] = np.dot(vec_a, vec_b)

        # a lambda that calculates indexes inside the gram matrix
        # and the corresponindg laplacian given a node and a level
        pick = lambda node, level:\
            (Nx[level][node],
                laplacian(Ax[Nx[level][node], :][:, Nx[level][node]]))\
            if node < nx else\
            ([idx+nx for idx in Ny[level][node-nx]],
                laplacian(Ay[Ny[level][node-nx], :][:, Ny[level][node-nx]]))

        for l in range(1, self._L+1):
            gram_matrix = np.empty(shape=(gram_matrix_size, gram_matrix_size))
            for i in range(0, gram_matrix_size):
                # calculate the correct indexes of neighbors
                # and the corresponding laplacian
                (idx_i, La) = pick(i, l)
                for j in range(0, gram_matrix_size):
                    (idx_j, Lb) = pick(j, l)

                    # calculate the corresponding gram matrix
                    idx_ij = idx_i + idx_j
                    extracted_gm = gram_matrix[idx_ij, :][:, idx_ij]

                    # calculate the core for this step
                    gram_matrix[i, j] =\
                        self._generalized_FLG_core_(La, Lb, extracted_gm)

        return self._generalized_FLG_core_(Lx, Ly, gram_matrix)

    def _generalized_FLG_core_(self, Lx, Ly, gram_matrix):
        """FLG core calculation for the multiscale gaussian.

        Parameters
        ----------
        L_{x,y} (np.array)
            Laplacians of graph {x,y}.

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
        v, w = np.real(v), np.real(w)

        # keep only the positive
        vpos = np.where(v > 0)

        # calculate the Q matrix
        Q = np.multiply(np.square(v[vpos]), w.T[vpos].T)
        Qx = Q[0:nx, :]
        Qy = Q[nx:, :]

        # Calculate the S matrices
        Tx = multi_dot((Qx.T, inv(np.add(Lx, heta*np.eye(Lx.shape[0]))), Qx))
        Ty = multi_dot((Qy.T, inv(np.add(Ly, heta*np.eye(Ly.shape[0]))), Qy))

        Sx = Tx + self._gamma*np.eye(Tx.shape[0])
        Sy = Ty + self._gamma*np.eye(Ty.shape[0])

        # A small lambda to calculate ^ 1/4
        quatre = lambda x: np.sqrt(np.sqrt(x))

        # !!! Overflow problem!: det(inv(Sx)+inv(Sy))=0.0
        # Need to solve. Not all executions are fatal
        # Caclulate the kernel nominator
        k_nom = np.sqrt(det(abs(np.multiply(inv(inv(Sx)+inv(Sy)), 2))))

        # Caclulate the kernel denominator
        k_denom = quatre(abs(det(Sx)*det(Sy)))

        return k_nom/k_denom
