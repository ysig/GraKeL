"""The pyramid match kernel as in :cite:`Nikolentzos2017MatchingNE`."""
import collections
import warnings

import numpy as np

from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigs

from grakel.graph import Graph
from grakel.kernels import kernel


class pyramid_match(kernel):
    """Pyramid match kernel class.

    Kernel defined in :cite:`Nikolentzos2017MatchingNE`

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
    _L : int
        Defines the histogram level of the pyramid.

    _d : int
        The dimension of the hypercube.

    _with_labels : bool
        Defines if to use labels in the calculation of the `pyramid_match`
        kernel.

    _num_labels : int
        The number of distinct labels, on the fit data.

    _labels : dict
        A dictionary of label enumeration, made from fitted data.

    """

    _graph_format = "adjacency"

    def __init__(self, **kargs):
        """Initialise a `pyramid_match` kernel."""
        self._valid_parameters |= {"L", "d", "with_labels"}
        super(pyramid_match, self).__init__(**kargs)

        self._with_labels = kargs.get("with_labels", True)
        if type(self._with_labels) != bool:
            raise ValueError('with labels must be a boolean variable')

        self._L = kargs.get("L", 4)
        if self._L < 0:
            raise ValueError('L: the number of levels must be bigger equal' +
                             ' to 0')

        self._d = kargs.get("d", 6)
        if self._d < 1:
            raise ValueError('d: hypercube dimension must be bigger than 1')

    def parse_input(self, X):
        """Parse and create features for shortest_path kernel.

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
        H : list
            A list of lists of Histograms for all levels for each graph.

        """
        if not isinstance(X, collections.Iterable):
            raise ValueError('input must be an iterable\n')
        else:
            i = 0
            Us = []
            if self._with_labels:
                Ls = []
            for (idx, x) in enumerate(iter(X)):
                if isinstance(x, collections.Iterable)\
                        and len(x) in [0, 2, 3]:
                    if len(x) == 0:
                        warnings.warn('Ignoring empty element on index: ' +
                                      str(idx))
                        continue
                    else:
                        x = Graph(x[0], x[1], {}, self._graph_format)
                elif not type(x) is Graph:
                    raise ValueError('each element of X must be either a ' +
                                     'graph object or a list with at least ' +
                                     'a graph like object and node labels ' +
                                     'dict \n')
                A = x.get_adjacency_matrix()
                if self._with_labels:
                    L = x.get_labels(purpose="adjacency")
                i += 1
                if A.shape[0] == 0:
                    Us.append(np.zeros((1, self._d)))
                else:
                    # Perform eigenvalue decomposition.
                    # Rows of matrix U correspond to vertex representations
                    # Embed vertices into the d-dimensional space
                    if A.shape[0] > self._d+1:
                        # If size of graph smaller than d, pad with zeros
                        Lambda, U = eigs(csr_matrix(A, dtype=np.float),
                                         k=self._d, ncv=10*self._d)
                        idx = Lambda.argsort()[::-1]
                        U = U[:, idx]
                    else:
                        Lambda, U = np.linalg.eig(A)
                        idx = Lambda.argsort()[::-1]
                        U = U[:, idx]
                        U = U[:, :self._d]
                    # Replace all components by their absolute values
                    U = np.absolute(U)
                    Us.append((A.shape[0], U))
                if self._with_labels:
                    Ls.append(L)

        if i == 0:
            raise ValueError('parsed input is empty')

        if self._method_calling in [1, 2]:
            # Map labels to values between 0 and |L|-1
            # where |L| is the number of distinct labels
            if self._with_labels:
                self._num_labels = 0
                self._labels = set()
                for L in Ls:
                    self._labels |= set(L.values())
                self._num_labels = len(self._labels)
                self._labels = {l: i for (i, l) in enumerate(self._labels)}

        # Calculate histograms
        if self._with_labels:
            return self._histogram_calculation(Us, Ls)
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

        Returns
        -------
        Hs : list
            List of histograms for each graph.

        """
        Hs = list()
        if len(args) == 0:
            for (i, (n, u)) in enumerate(Us):
                if n > 0:
                    du = list()
                    for j in range(self._L):
                        # Number of cells along each dimension at level j
                        k = 2**j
                        # Determines the cells in which each vertex lies
                        # along each dimension since nodes lie in the unit
                        # hypercube in R^d
                        D = np.zeros((self._d, k))
                        T = np.floor(Us[i]*k)
                        T[np.where(T == k)] = k-1
                        for p in range(u.shape[0]):
                            if p >= n:
                                break
                            for q in range(u.shape[1]):
                                # Identify the cell into which the i-th
                                # vertex lies and increase its value by 1
                                D[q, int(T[p, q])] = D[q, int(T[p, q])] + 1
                            du.append(D)
                    Hs.append(du)

        elif len(args) > 0:
            Ls = args[0]
            for (i, ((n, u), L)) in enumerate(zip(Us, Ls)):
                du = list()
                if n > 0:
                    for j in range(self._L):
                        # Number of cells along each dimension at level j
                        k = 2**j
                        # To store the number of vertices that are assigned
                        # a specific label and lie in each of the 2^j cells
                        # of each dimension at level j
                        D = np.zeros((self._d*self._num_labels, k))
                        T = np.floor(u*k)
                        T[np.where(T == k)] = k-1
                        for p in range(u.shape[0]):
                            if p >= n:
                                break
                            for q in range(u.shape[1]):
                                # Identify the cell into which the i-th
                                # vertex lies and increase its value by 1
                                D[self._labels[L[p]]*self._d + q,
                                  int(T[p, q])] = \
                                    D[self._labels[L[p]]*self._d + q,
                                      int(T[p, q])] + 1
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
        intersec = np.zeros(self._L)
        for (p, xp, yp) in zip(range(self._L), x, y):
            # Calculate histogram intersection
            # (eq. 6 in :cite:`Nikolentzos2017MatchingNE`)
            intersec[p] = np.sum(np.minimum(xp, yp))
            k += intersec[self._L-1]
            for p in range(self._L-1):
                # Computes the new matches that occur at level p.
                # These matches weight less than those that occur at
                # higher levels (e.g. p+1 level)
                k += (1.0/(2**(self._L-p-1)))*(intersec[p]-intersec[p+1])
        return k
