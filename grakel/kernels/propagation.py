"""The propagation kernel as defined in :cite:`Neumann2015PropagationKE`."""
import collections
import random
import warnings

import numpy as np

from grakel.graph import Graph
from grakel.kernels import kernel

# Python 2/3 cross-compatibility import
from six import itervalues

default_random_seed_value = 1235476565


class propagation(kernel):
    r"""The Propagation kernel for fully labeled graphs.

    See :cite:`Neumann2015PropagationKE`: Algorithms 1, 3, p. 216, 221.

    Parameters
    ----------
    t_max : int, default=5
        Maximum number of iterations.

    w : int, default=10
        Bin width.

    M : str, default="TV"
        The preserved distance metric (on local sensitive hashing):
            - "H": hellinger
            - "L1": l1-norm
            - "L2": l2-norm
            - "TV": total-variation

    base_kernel : function (np.array, np.array -> number),
        default=:math:`f(x,y)=\sum_{i} x_{i}*y_{i}`
        A base_kernel between two 1-dimensional numpy arrays of numbers
        that outputs a number.

    Attributes
    ----------
    _enum_labels : dict
        Holds the enumeration of the input labels.

    _parent_labels : set
        Holds a set of the input labels.

    _M : str
        The preserved distance metric (on local sensitive hashing).

    _tmax : int
        Holds the maximum number of iterations.

    _w : int
        Holds the bin width.

    _base_kernel : function (np.array, np.array -> number)
        A base_kernel between two 1-dimensional numpy arrays of numbers
        that outputs a number.

    """

    _graph_format = "adjacency"

    def __init__(self, **kargs):
        """Initialise a subtree_wl kernel."""
        self._valid_parameters |= {"random_seed", "base_kernel"
                                   "M", "t_max", "w"}
        super(propagation, self).__init__(**kargs)

        random.seed(int(kargs.get("random_seed", default_random_seed_value)))

        self._M = kargs.get("M", "TV")
        if type(self._M) is not str or self._M not in ["H", "L1", "L2", "TV"]:
            raise ValueError('Metric type must be a str, one of "H", "L1", ' +
                             '"L2", "TV".')

        self._tmax = kargs.get("t_max", 5)
        if type(self._tmax) is not int or self._tmax <= 0:
            raise ValueError('The number of iterations must be a ' +
                             'positive integer.')

        self._w = kargs.get("w", 10)
        if type(self._w) is not int or self._w <= 0:
            raise ValueError('The bin width must be a positive integer.')

        self._base_kernel = \
            kargs.get("base_kernel",
                      lambda x, y: np.sum(np.array(x)*np.array(y)))
        if not callable(self._base_kernel):
            raise ValueError('The base kernel must be callable.')

    def pairwise_operation(self, x, y):
        """Calculate the kernel value between two elements.

        Parameters
        ----------
        x, y: list
            Inverse label dictionaries.

        Returns
        -------
        kernel : number
            The kernel value.

        """
        k = .0
        for t in range(self._tmax):
            max_idx = min(len(x[t]), len(y[t]))
            k += self._base_kernel(x[t][:max_idx], y[t][:max_idx])
        return k

    def parse_input(self, X):
        """Parse and create features for graphlet_sampling kernel.

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
        local_values : dict
            A dictionary of pairs between each input graph and a bins where the
            sampled graphlets have fallen.

        """
        if not isinstance(X, collections.Iterable):
            raise ValueError('input must be an iterable\n')
        else:
            i = -1
            transition_matrix = dict()
            labels = set()
            L = list()
            for (idx, x) in enumerate(iter(X)):
                is_iter = isinstance(x, collections.Iterable)
                if is_iter:
                    x = list(x)
                if is_iter and len(x) in [0, 2, 3, 4]:
                    if len(x) == 0:
                        warnings.warn('Ignoring empty element on ' +
                                      'index: '+str(idx))
                        continue
                    if len(x) == 2 and type(x[0]) is Graph:
                        g, T = x
                    else:
                        g = Graph(x[0], x[1], {}, self._graph_format)
                        if len(x) == 4:
                            T = x[3]
                        else:
                            T = None
                elif type(x) is Graph:
                    g, T = x, None
                else:
                    raise ValueError('Each element of X must be either a ' +
                                     'Graph or an iterable with at least 2 ' +
                                     'and at most 4 elements\n')

                if T is not None:
                    if T.shape[0] != T.shape[1]:
                        raise ValueError('Transition matrix on index' +
                                         ' ' + str(idx) + 'must be ' +
                                         'a square matrix.')
                    if T.shape[0] != g.nv():
                        raise ValueError('Propagation matrix must ' +
                                         'have the same dimension ' +
                                         'as the number of vertices.')
                else:
                    T = g.get_adjacency_matrix()
                T = (T > 0).astype(int)

                i += 1
                transition_matrix[i] = T
                label = g.get_labels(purpose='adjacency')
                labels |= set(itervalues(label))
                L.append((g.nv(), label))

            if i == -1:
                raise ValueError('Parsed input is empty')

            # The number of parsed graphs
            n = i+1

            # enumerate labels
            if self._method_calling == 1:
                enum_labels = {l: i for (i, l)
                               in enumerate(sorted(list(labels)))}
                self._enum_labels = enum_labels
                self._parent_labels = labels
            elif self._method_calling == 3:
                new_elements = labels - self._parent_labels
                if len(new_elements) > 0:
                    new_enum_labels = {l: i for (i, l) in
                                       enumerate(len(self._enum_labels),
                                       sorted(list(new_elements)))}
                    enum_labels = dict(self._enum_labels, **new_enum_labels)
                else:
                    enum_labels = self._enum_labels

            # make a matrix for all graphs that contains label vectors
            P = dict()
            for (k, (nv, label)) in enumerate(L):
                P[k] = np.zeros(shape=(nv, len(enum_labels)))
                for j in range(nv):
                    P[k][j, enum_labels[label[j]]] = 1

            # feature vectors
            phi = {k: dict() for k in range(n)}
            for t in range(self._tmax):
                # for each graph hash P and produce the feature vectors
                for k in range(n):
                    phi[k][t] = np.bincount(calculate_LSH(P[k],
                                                          self._w, self._M))

                # calculate the Propagation matrix if needed
                if t < self._tmax-1:
                    for k in range(n):
                        P[k] = np.dot(transition_matrix[k], P[k])

            return [phi[k] for k in range(n)]


def calculate_LSH(X, w, M):
    """Calculate Local Sensitive Hashing needed for propagation kernels.

    See :cite:`Neumann2015PropagationKE`, p.12.
    Parameters
    ----------
    X : np.array
        A float array of shape (N, D) with N vertices and D features.

    w : int
        Bin width.

    M : str
        The preserved distance metric:
           - "H": hellinger
           - "L1": l1-norm
           - "L2": l2-norm
           - "TV": total-variation

    Returns
    -------
    lsh : np.array.
        The local sensitive hash coresponding to each vertex.

    """
    if len(X.shape) != 2:
        raise ValueError('X must be an N x D array')
    D = X.shape[1]

    if M not in ["H", "L1", "L2", "TV"]:
        raise ValueError('Invalid metric type')

    if M is "H":
        X = np.sqrt(X)

    # simple normal
    u = np.random.randn(D)
    if M in ["TV", "L1"]:
        # cauchy
        u = np.divide(u, np.random.randn(D))

    # random offset
    b = w*np.random.rand()

    # hash
    h = np.floor((np.dot(X, u)+b)/w)
    return np.unique(h, return_inverse=True)[1]
