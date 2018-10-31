"""The pyramid match kernel as in :cite:`nikolentzos2017matching`."""
# Author: Ioannis Siglidis <y.siglidis@gmail.com>
# License: BSD 3 clause
import collections
import warnings

import numpy as np

from itertools import chain

from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigs

from grakel.graph import Graph
from grakel.kernels import Kernel

# Python 2/3 cross-compatibility import
from six import itervalues
from six import iteritems


class PyramidMatch(Kernel):
    """Pyramid match kernel class.

    Kernel defined in :cite:`nikolentzos2017matching`

    Parameters
    ----------
    with_labels : bool, default=True
        A flag that determines if the kernel computation will consider labels.

    L : int, default=4
        Pyramid histogram level.

    d : int, default=6
        The dimension of the hypercube.

    Attributes
    ----------
    _num_labels : int
        The number of distinct labels, on the fit data.

    _labels : dict
        A dictionary of label enumeration, made from fitted data.

    """

    _graph_format = "adjacency"

    def __init__(self, n_jobs=None,
                 normalize=False,
                 verbose=False,
                 with_labels=True,
                 L=4,
                 d=6):
        """Initialise a `pyramid_match` kernel."""
        super(PyramidMatch, self).__init__(n_jobs=n_jobs,
                                           normalize=normalize,
                                           verbose=verbose)

        self.with_labels = with_labels
        self.L = L
        self.d = d
        self._initialized.update({"d": False, "L": False, "with_labels": False})

    def initialize(self):
        """Initialize all transformer arguments, needing initialization."""
        super(PyramidMatch, self).initialize()

        if not self._initialized["with_labels"]:
            if type(self.with_labels) != bool:
                raise TypeError('with labels must be a boolean variable')
            self._initialized["with_labels"] = True

        if not self._initialized["L"]:
            if type(self.L) is not int or self.L < 0:
                raise TypeError('L: the number of levels must be an integer '
                                'bigger equal to 0')
            self._initialized["L"] = True

        if not self._initialized["d"]:
            if type(self.d) is not int or self.d < 1:
                raise TypeError('d: hypercube dimension must be an '
                                'integer bigger than 1')
            self._initialized["d"] = True

    def parse_input(self, X):
        """Parse and create features for pyramid_match kernel.

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
        H : list
            A list of lists of Histograms for all levels for each graph.

        """
        if not isinstance(X, collections.Iterable):
            raise TypeError('input must be an iterable\n')
        else:
            i = 0
            Us = []
            if self.with_labels:
                Ls = []
            for (idx, x) in enumerate(iter(X)):
                is_iter = isinstance(x, collections.Iterable)
                if is_iter:
                    x = list(x)
                if is_iter and (len(x) == 0 or (len(x) >= 1 and not self.with_labels) or
                                (len(x) >= 2 and self.with_labels)):
                    if len(x) == 0:
                        warnings.warn('Ignoring empty element on index: ' + str(idx))
                        continue
                    elif not self.with_labels:
                        x = Graph(x[0], {}, {}, self._graph_format)
                    else:
                        x = Graph(x[0], x[1], {}, self._graph_format)
                elif not type(x) is Graph:
                    raise TypeError('each element of X must be either a graph object or a list with '
                                    'at least a graph like object and node labels dict \n')
                A = x.get_adjacency_matrix()
                if self.with_labels:
                    L = x.get_labels(purpose="adjacency")
                i += 1
                if A.shape[0] == 0:
                    Us.append(np.zeros((1, self.d)))
                else:
                    # Perform eigenvalue decomposition.
                    # Rows of matrix U correspond to vertex representations
                    # Embed vertices into the d-dimensional space
                    if A.shape[0] > self.d+1:
                        # If size of graph smaller than d, pad with zeros
                        Lambda, U = eigs(csr_matrix(A, dtype=np.float),
                                         k=self.d, ncv=10*self.d)
                        idx = Lambda.argsort()[::-1]
                        U = U[:, idx]
                    else:
                        Lambda, U = np.linalg.eig(A)
                        idx = Lambda.argsort()[::-1]
                        U = U[:, idx]
                        U = U[:, :self.d]
                    # Replace all components by their absolute values
                    U = np.absolute(U)
                    Us.append((A.shape[0], U))
                if self.with_labels:
                    Ls.append(L)

        if i == 0:
            raise ValueError('parsed input is empty')

        if self.with_labels:
            # Map labels to values between 0 and |L|-1
            # where |L| is the number of distinct labels
            if self._method_calling in [1, 2]:
                self._num_labels = 0
                self._labels = set()
                for L in Ls:
                    self._labels |= set(itervalues(L))
                self._num_labels = len(self._labels)
                self._labels = {l: i for (i, l) in enumerate(self._labels)}
                return self._histogram_calculation(Us, Ls, self._labels)

            elif self._method_calling == 3:
                labels = set()
                for L in Ls:
                    labels |= set(itervalues(L))
                rest_labels = labels - set(self._labels.keys())
                nouveau_labels = dict(chain(iteritems(self._labels),
                                      ((j, i) for (i, j) in enumerate(rest_labels, len(self._labels)))))
                return self._histogram_calculation(Us, Ls, nouveau_labels)
        else:
            return self._histogram_calculation(Us)

    def _histogram_calculation(self, Us, *args):
        """Calculate histograms.

        Parameters
        ----------
        Us : list
            List of tuples with the first element corresponding to the
            number of vertices of a graph and the second to it's
            corresponding to vertex embeddings on the d-dimensional space.

        Ls : list, optional
            List of labels corresponding to each graph.
            If provided the histograms are calculated with labels.

        Labels : dict, optional
            A big dictionary with enumeration of labels.

        Returns
        -------
        Hs : list
            List of histograms for each graph.

        """
        Hs = list()
        if len(args) == 0:
            for (i, (n, u)) in enumerate(Us):
                du = list()
                if n > 0:
                    for j in range(self.L):
                        # Number of cells along each dimension at level j
                        k = 2**j
                        # Determines the cells in which each vertex lies
                        # along each dimension since nodes lie in the unit
                        # hypercube in R^d
                        D = np.zeros((self.d, k))
                        T = np.floor(u*k)
                        T[np.where(T == k)] = k-1
                        for p in range(u.shape[0]):
                            if p >= n:
                                break
                            for q in range(u.shape[1]):
                                # Identify the cell into which the i-th
                                # vertex lies and increase its value by 1
                                D[q, int(T[p, q])] += 1
                        du.append(D)
                Hs.append(du)

        elif len(args) > 0:
            Ls = args[0]
            Labels = args[1]
            num_labels = len(Labels)
            for (i, ((n, u), L)) in enumerate(zip(Us, Ls)):
                du = list()
                if n > 0:
                    for j in range(self.L):
                        # Number of cells along each dimension at level j
                        k = 2**j
                        # To store the number of vertices that are assigned
                        # a specific label and lie in each of the 2^j cells
                        # of each dimension at level j
                        D = np.zeros((self.d*num_labels, k))
                        T = np.floor(u*k)
                        T[np.where(T == k)] = k-1
                        for p in range(u.shape[0]):
                            if p >= n:
                                break
                            for q in range(u.shape[1]):
                                # Identify the cell into which the i-th
                                # vertex lies and increase its value by 1
                                D[Labels[L[p]]*self.d + q, int(T[p, q])] += 1
                        du.append(D)
                Hs.append(du)
        return Hs

    def pairwise_operation(self, x, y):
        """Calculate a pairwise kernel between two elements.

        Parameters
        ----------
        x, y : dict
            Histograms as produced by `parse_input`.

        Returns
        -------
        kernel : number
            The kernel value.

        """
        k = 0
        if len(x) != 0 and len(y) != 0:
            intersec = np.zeros(self.L)
            for (p, xp, yp) in zip(range(self.L), x, y):
                # Calculate histogram intersection
                # (eq. 6 in :cite:`nikolentzos2017matching`)
                if xp.shape[0] < yp.shape[0]:
                    xpp, ypp = xp, yp[:xp.shape[0], :]
                elif yp.shape[0] < xp.shape[0]:
                    xpp, ypp = xp[:yp.shape[0], :], yp
                else:
                    xpp, ypp = xp, yp
                intersec[p] = np.sum(np.minimum(xpp, ypp))
                k += intersec[self.L-1]
                for p in range(self.L-1):
                    # Computes the new matches that occur at level p.
                    # These matches weight less than those that occur at
                    # higher levels (e.g. p+1 level)
                    k += (1.0/(2**(self.L-p-1)))*(intersec[p]-intersec[p+1])
        return k
