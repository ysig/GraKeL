"""The propagation kernel as defined in :cite:`neumann2015propagation`."""
# Author: Ioannis Siglidis <y.siglidis@gmail.com>
# License: BSD 3 clause
import collections
import warnings

import numpy as np

from itertools import chain
from collections import Counter
from numbers import Real

from sklearn.utils import check_random_state

from grakel.graph import Graph
from grakel.kernels import Kernel

# Python 2/3 cross-compatibility import
from six import itervalues
from six import iteritems
from six.moves import filterfalse


def _dot(x, y):
    return sum(x[k]*y[k] for k in x)


class Propagation(Kernel):
    r"""The Propagation kernel for fully labeled graphs.

    See :cite:`neumann2015propagation`: Algorithms 1, 3, p. 216, 221.

    Parameters
    ----------
    t_max : int, default=5
        Maximum number of iterations.

    w : int, default=0.01
        Bin width.

    M : str, default="TV"
        The preserved distance metric (on local sensitive hashing):
            - "H": hellinger
            - "TV": total-variation

    metric : function (Counter, Counter -> number),
        default=:math:`f(x,y)=\sum_{i} x_{i}*y_{i}`
        A metric between two 1-dimensional numpy arrays of numbers that outputs a number.
        It must consider the case where the keys of y are not in x, when different features appear
        at transform.

    random_state :  RandomState or int, default=None
        A random number generator instance or an int to initialize a RandomState as a seed.

    Attributes
    ----------
    _enum_labels : dict
        Holds the enumeration of the input labels.

    _parent_labels : set
        Holds a set of the input labels.

    random_state_ : RandomState
        A RandomState object handling all randomness of the class.

    """

    _graph_format = "adjacency"
    attr_ = False

    def __init__(self,
                 n_jobs=None,
                 verbose=False,
                 normalize=False,
                 random_state=None,
                 metric=_dot,
                 M="TV",
                 t_max=5,
                 w=0.01):
        """Initialise a propagation kernel."""
        super(Propagation, self).__init__(n_jobs=n_jobs,
                                          verbose=verbose,
                                          normalize=normalize)

        self.random_state = random_state
        self.M = M
        self.t_max = t_max
        self.w = w
        self.metric = metric
        self._initialized.update({"M": False, "t_max": False, "w": False,
                                  "random_state": False, "metric": False})

    def initialize(self):
        """Initialize all transformer arguments, needing initialization."""
        super(Propagation, self).initialize()

        if not self._initialized["random_state"]:
            self.random_state_ = check_random_state(self.random_state)
            self._initialized["random_state"] = True

        if not self._initialized["metric"]:
            if (type(self.M) is not str or
                    (self.M not in ["H", "TV"] and not self.attr_) or
                    (self.M not in ["L1", "L2"] and self.attr_)):
                if self.attr_:
                    raise TypeError('Metric type must be a str, one of "L1", "L2"')
                else:
                    raise TypeError('Metric type must be a str, one of "H", "TV"')

            if not self.attr_:
                self.take_sqrt_ = self.M == "H"

            self.take_cauchy_ = self.M in ["TV", "L1"]
            self._initialized["metric"] = True

        if not self._initialized["t_max"]:
            if type(self.t_max) is not int or self.t_max <= 0:
                raise TypeError('The number of iterations must be a ' +
                                'positive integer.')
            self._initialized["t_max"] = True

        if not self._initialized["w"]:
            if not isinstance(self.w, Real) and self.w <= 0:
                raise TypeError('The bin width must be a positive number.')
            self._initialized["w"] = True

        if not self._initialized["metric"]:
            if not callable(self.metric):
                raise TypeError('The base kernel must be callable.')
            self._initialized["metric"] = True

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
        return sum(self.metric(x[t], y[t]) for t in range(self.t_max))

    def parse_input(self, X):
        """Parse and create features for the propation kernel.

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
                        raise TypeError('Transition matrix on index' +
                                        ' ' + str(idx) + 'must be ' +
                                        'a square matrix.')
                    if T.shape[0] != g.nv():
                        raise TypeError('Propagation matrix must ' +
                                        'have the same dimension ' +
                                        'as the number of vertices.')
                else:
                    T = g.get_adjacency_matrix()

                i += 1
                transition_matrix[i] = (T.T / np.sum(T, axis=1)).T
                label = g.get_labels(purpose='adjacency')
                try:
                    labels |= set(itervalues(label))
                except TypeError:
                    raise TypeError('For a non attributed kernel, labels should be hashable.')
                L.append((g.nv(), label))

            if i == -1:
                raise ValueError('Parsed input is empty')

            # The number of parsed graphs
            n = i+1

            # enumerate labels
            if self._method_calling == 1:
                enum_labels = {l: i for (i, l) in enumerate(list(labels))}
                self._enum_labels = enum_labels
                self._parent_labels = labels
            elif self._method_calling == 3:
                new_elements = labels - self._parent_labels
                if len(new_elements) > 0:
                    new_enum_labels = iter((l, i) for (i, l) in
                                           enumerate(list(new_elements), len(self._enum_labels)))
                    enum_labels = dict(chain(iteritems(self._enum_labels), new_enum_labels))
                else:
                    enum_labels = self._enum_labels

            # make a matrix for all graphs that contains label vectors
            P, data, indexes = dict(), list(), [0]
            for (k, (nv, label)) in enumerate(L):
                data += [(indexes[-1] + j, enum_labels[label[j]]) for j in range(nv)]
                indexes.append(indexes[-1] + nv)

            # Initialise the on hot vector
            rows, cols = zip(*data)
            P = np.zeros(shape=(indexes[-1], len(enum_labels)))
            P[rows, cols] = 1
            dim_orig = len(self._enum_labels)

            # feature vectors
            if self._method_calling == 1:
                # simple normal
                self._u, self._b, self._hd = list(), list(), list()
                for t in range(self.t_max):
                    u = self.random_state_.randn(len(enum_labels))

                    if self.take_cauchy_:
                        # cauchy
                        u = np.divide(u, self.random_state_.randn(len(enum_labels)))

                    self._u.append(u)
                    # random offset
                    self._b.append(self.w*self.random_state_.rand())

                phi = {k: dict() for k in range(n)}
                for t in range(self.t_max):
                    # for hash all graphs inside P and produce the feature vectors
                    hashes = self.calculate_LSH(P, self._u[t], self._b[t])
                    hd = dict((j, i) for i, j in enumerate(set(np.unique(hashes))))
                    self._hd.append(hd)
                    features = np.vectorize(lambda i: hd[i])(hashes)

                    # Accumulate the results.
                    for k in range(n):
                        phi[k][t] = Counter(features[indexes[k]:indexes[k+1]])

                    # calculate the Propagation matrix if needed
                    if t < self.t_max-1:
                        for k in range(n):
                            start, end = indexes[k:k+2]
                            P[start:end, :] = np.dot(transition_matrix[k], P[start:end, :])

                return [phi[k] for k in range(n)]

            elif (self._method_calling == 3 and dim_orig >= len(enum_labels)):
                phi = {k: dict() for k in range(n)}
                for t in range(self.t_max):
                    # for hash all graphs inside P and produce the feature vectors
                    hashes = self.calculate_LSH(P, self._u[t], self._b[t])
                    hd = dict(chain(
                            iteritems(self._hd[t]),
                            iter((j, i) for i, j in enumerate(
                                    filterfalse(lambda x: x in self._hd[t],
                                                np.unique(hashes)),
                                    len(self._hd[t])))))

                    features = np.vectorize(lambda i: hd[i])(hashes)

                    # Accumulate the results.
                    for k in range(n):
                        phi[k][t] = Counter(features[indexes[k]:indexes[k+1]])

                    # calculate the Propagation matrix if needed
                    if t < self.t_max-1:
                        for k in range(n):
                            start, end = indexes[k:k+2]
                            P[start:end, :] = np.dot(transition_matrix[k], P[start:end, :])

                return [phi[k] for k in range(n)]

            else:
                cols = np.array(cols)
                vertices = np.where(cols < dim_orig)[0]
                vertices_p = np.where(cols >= dim_orig)[0]
                nnv = len(enum_labels) - dim_orig
                phi = {k: dict() for k in range(n)}
                for t in range(self.t_max):
                    # hash all graphs inside P and produce the feature vectors
                    hashes = self.calculate_LSH(P[vertices, :dim_orig],
                                                self._u[t], self._b[t])

                    hd = dict(chain(
                            iteritems(self._hd[t]),
                            iter((j, i) for i, j in enumerate(
                                    filterfalse(lambda x: x in self._hd[t],
                                                np.unique(hashes)),
                                    len(self._hd[t])))))

                    features = np.vectorize(lambda i: hd[i], otypes=[int])(hashes)

                    # for each the new labels graph hash P and produce the feature vectors
                    u = self.random_state_.randn(nnv)
                    if self.take_cauchy_:
                        # cauchy
                        u = np.divide(u, self.random_state_.randn(nnv))

                    u = np.hstack((self._u[t], u))

                    # calculate hashes for the remaining
                    hashes = self.calculate_LSH(P[vertices_p, :], u, self._b[t])
                    hd = dict(chain(iteritems(hd), iter((j, i) for i, j in enumerate(hashes, len(hd)))))

                    features_p = np.vectorize(lambda i: hd[i], otypes=[int])(hashes)

                    # Accumulate the results
                    for k in range(n):
                        A = Counter(features[np.logical_and(
                            indexes[k] <= vertices, vertices <= indexes[k+1])])
                        B = Counter(features_p[np.logical_and(
                            indexes[k] <= vertices_p, vertices_p <= indexes[k+1])])
                        phi[k][t] = A + B

                    # calculate the Propagation matrix if needed
                    if t < self.t_max-1:
                        for k in range(n):
                            start, end = indexes[k:k+2]
                            P[start:end, :] = np.dot(transition_matrix[k], P[start:end, :])

                        Q = np.all(P[:, dim_orig:] > 0, axis=1)
                        vertices = np.where(~Q)[0]
                        vertices_p = np.where(Q)[0]

                return [phi[k] for k in range(n)]

    def calculate_LSH(self, X, u, b):
        """Calculate Local Sensitive Hashing needed for propagation kernels.

        See :cite:`neumann2015propagation`, p.12.

        Parameters
        ----------
        X : np.array
            A float array of shape (N, D) with N vertices and D features.

        u : np.array, shape=(D, 1)
            A projection vector.

        b : float
            An offset (times w).

        Returns
        -------
        lsh : np.array.
            The local sensitive hash coresponding to each vertex.

        """
        if self.take_sqrt_:
            X = np.sqrt(X)

        # hash
        return np.floor((np.dot(X, u)+b)/self.w)


class PropagationAttr(Propagation):
    r"""The Propagation kernel for fully attributed graphs.

    See :cite:`neumann2015propagation`: Algorithms 1, 3, p. 216, 221.

    Parameters
    ----------
    t_max : int, default=5
        Maximum number of iterations.

    w : int, default=0.01
        Bin width.

    M : str, default="TV"
        The preserved distance metric (on local sensitive hashing):
            - "L1": l1-norm
            - "L2": l2-norm

    metric : function (np.array, np.array -> number),
        default=:math:`f(x,y)=\sum_{i} x_{i}*y_{i}`
        A metric between two 1-dimensional numpy arrays of numbers
        that outputs a number.

    Attributes
    ----------
    M : str
        The preserved distance metric (on local sensitive hashing).

    tmax : int
        Holds the maximum number of iterations.

    w : int
        Holds the bin width.

    metric : function (np.array, np.array -> number)
        A metric between two 1-dimensional numpy arrays of numbers
        that outputs a number.

    """

    _graph_format = "adjacency"
    attr_ = True

    def __init__(self,
                 n_jobs=None,
                 verbose=False,
                 normalize=False,
                 random_state=None,
                 metric=_dot,
                 M="L1",
                 t_max=5,
                 w=4):
        """Initialise a propagation kernel."""
        super(PropagationAttr, self).__init__(n_jobs=n_jobs,
                                              verbose=verbose,
                                              normalize=normalize,
                                              random_state=random_state,
                                              metric=metric,
                                              M=M,
                                              t_max=t_max,
                                              w=w)

    def initialize(self):
        """Initialize all transformer arguments, needing initialization."""
        super(PropagationAttr, self).initialize()

    def parse_input(self, X):
        """Parse and create features for the attributed propation kernel.

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
            # The number of parsed graphs
            n = 0
            transition_matrix = dict()
            indexes = [0]
            Attr = list()
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
                        raise TypeError('Transition matrix on index' +
                                        ' ' + str(idx) + 'must be ' +
                                        'a square matrix.')
                    if T.shape[0] != g.nv():
                        raise TypeError('Propagation matrix must ' +
                                        'have the same dimension ' +
                                        'as the number of vertices.')
                else:
                    T = g.get_adjacency_matrix()

                nv = g.nv()
                transition_matrix[n] = (T.T / np.sum(T, axis=1)).T
                attr = g.get_labels(purpose="adjacency")
                try:
                    attributes = np.array([attr[j] for j in range(nv)])
                except TypeError:
                    raise TypeError('All attributes of a single graph should have the same dimension.')

                Attr.append(attributes)
                indexes.append(indexes[-1] + nv)
                n += 1
            try:
                P = np.vstack(Attr)
            except ValueError:
                raise ValueError('Attribute dimensions should be the same, for all graphs')

            if self._method_calling == 1:
                self._dim = P.shape[1]
            else:
                if self._dim != P.shape[1]:
                    raise ValueError('transform attribute vectors should'
                                     'have the same dimension as in fit')

            if n == 0:
                raise ValueError('Parsed input is empty')

            # feature vectors
            if self._method_calling == 1:
                # simple normal
                self._u, self._b, self._hd = list(), list(), list()
                for t in range(self.t_max):
                    u = self.random_state_.randn(self._dim)
                    if self.take_cauchy_:
                        # cauchy
                        u = np.divide(u, self.random_state_.randn(self._dim))

                    self._u.append(u)
                    # random offset
                    self._b.append(self.w*self.random_state_.randn(self._dim))

                phi = {k: dict() for k in range(n)}
                for t in range(self.t_max):
                    # for hash all graphs inside P and produce the feature vectors
                    hashes = self.calculate_LSH(P, self._u[t], self._b[t]).tolist()

                    hd = {j: i for i, j in enumerate({tuple(l) for l in hashes})}
                    self._hd.append(hd)

                    features = np.array([hd[tuple(l)] for l in hashes])

                    # Accumulate the results.
                    for k in range(n):
                        phi[k][t] = Counter(features[indexes[k]:indexes[k+1]].flat)

                    # calculate the Propagation matrix if needed
                    if t < self.t_max-1:
                        for k in range(n):
                            start, end = indexes[k:k+2]
                            P[start:end, :] = np.dot(transition_matrix[k], P[start:end, :])

                return [phi[k] for k in range(n)]

            if self._method_calling == 3:
                phi = {k: dict() for k in range(n)}
                for t in range(self.t_max):
                    # for hash all graphs inside P and produce the feature vectors
                    hashes = self.calculate_LSH(P, self._u[t], self._b[t]).tolist()

                    hd = dict(chain(
                            iteritems(self._hd[t]),
                            iter((j, i) for i, j in enumerate(
                                    filterfalse(lambda x: x in self._hd[t],
                                                {tuple(l) for l in hashes}),
                                    len(self._hd[t])))))

                    features = np.array([hd[tuple(l)] for l in hashes])

                    # Accumulate the results.
                    for k in range(n):
                        phi[k][t] = Counter(features[indexes[k]:indexes[k+1]])

                    # calculate the Propagation matrix if needed
                    if t < self.t_max-1:
                        for k in range(n):
                            start, end = indexes[k:k+2]
                            P[start:end, :] = np.dot(transition_matrix[k], P[start:end, :])

                return [phi[k] for k in range(n)]

    def calculate_LSH(self, X, u, b):
        """Calculate Local Sensitive Hashing needed for propagation kernels.

        See :cite:`neumann2015propagation`, p.12.

        Parameters
        ----------
        X : np.array
            A float array of shape (N, D) with N vertices and D features.

        u : np.array, shape=(D, 1)
            A projection vector.

        b : float
            An offset (times w).

        Returns
        -------
        lsh : np.array.
            The local sensitive hash coresponding to each vertex.

        """
        return np.floor((X*u+b)/self.w)


if __name__ == '__main__':
    from grakel.datasets import fetch_dataset
    import argparse
    # Create an argument parser for the installer of pynauty
    parser = argparse.ArgumentParser(
        description='Measuring classification accuracy '
                    'on multiscale_laplacian_fast')

    parser.add_argument(
        '--dataset',
        help='choose the datset you want the tests to be executed',
        type=str,
        default=None
    )

    parser.add_argument(
        '--full',
        help='fit_transform the full graph',
        action="store_true")

    parser.add_argument(
        '--attr',
        help='define if the attributed kernel will be used',
        action="store_true")

    parser.add_argument(
        '--tmax',
        help='choose the datset you want the tests to be executed',
        type=int,
        default=5
    )

    parser.add_argument(
        '--w',
        help='choose the datset you want the tests to be executed',
        type=float,
        default=0.01
    )

    # Get the dataset name
    args = parser.parse_args()
    has_attributes = bool(args.attr)
    dataset_name = args.dataset
    tmax = args.tmax
    w = args.w
    full = bool(args.full)
    if dataset_name is None:
        if has_attributes:
            dataset_name = "BZR"
        else:
            dataset_name = "MUTAG"

    # The baseline dataset for node/edge-attributes
    dataset_attr = fetch_dataset(dataset_name,
                                 with_classes=True,
                                 prefer_attr_nodes=has_attributes,
                                 verbose=True)

    from tqdm import tqdm
    from time import time

    from sklearn.metrics import accuracy_score
    from sklearn.model_selection import KFold
    from sklearn import svm

    def sec_to_time(sec):
        """Print time in a correct format."""
        dt = list()
        days = int(sec // 86400)
        if days > 0:
            sec -= 86400*days
            dt.append(str(days) + " d")

        hrs = int(sec // 3600)
        if hrs > 0:
            sec -= 3600*hrs
            dt.append(str(hrs) + " h")

        mins = int(sec // 60)
        if mins > 0:
            sec -= 60*mins
            dt.append(str(mins) + " m")

        if sec > 0:
            dt.append(str(round(sec, 2)) + " s")
        return " ".join(dt)

    # Loads the Mutag dataset from:
    # https://ls11-www.cs.tu-dortmund.de/staff/morris/graphkerneldatasets
    # the biggest collection of benchmark datasets for graph_kernels.
    G, y = dataset_attr.data, dataset_attr.target
    C_grid = (10. ** np.arange(-7, 7, 2) / len(G)).tolist()

    stats = {"acc": list(), "time": list()}

    kf = KFold(n_splits=10, random_state=42, shuffle=True)
    niter = kf.get_n_splits(y)

    for (k, (train_index, test_index)) in tqdm(enumerate(kf.split(G, y)),
                                               total=niter):
        # Train-test split of graph data
        tri = train_index.tolist()
        tei = test_index.tolist()

        G_train, G_test = list(), list()
        y_train, y_test = list(), list()
        for (i, (g, t)) in enumerate(zip(G, y)):
            if len(tri) and i == tri[0]:
                G_train.append(g)
                y_train.append(t)
                tri.pop(0)
            elif len(tei) and i == tei[0]:
                G_test.append(g)
                y_test.append(t)
                tei.pop(0)

        start = time()
        if has_attributes:
            gk = PropagationAttr(M="L1", t_max=tmax, w=w, normalize=True)
        else:
            gk = Propagation(M="H", t_max=tmax, w=w, normalize=True)

        # Calculate the kernel matrix.
        if full:
            K = gk.fit_transform(G)
            K_train = K[train_index, :][:, train_index]
            K_test = K[test_index, :][:, train_index]
        else:
            K_train = gk.fit_transform(G_train)
            K_test = gk.transform(G_test)
        end = time()

        # Cross validation on C, variable
        acc = 0
        for c in C_grid:
            # Initialise an SVM and fit.
            clf = svm.SVC(kernel='precomputed', C=c)

            # Fit on the train Kernel
            clf.fit(K_train, y_train)

            # Predict and test.
            y_pred = clf.predict(K_test)

            # Calculate accuracy of classification.
            acc = max(acc, accuracy_score(y_test, y_pred))

        stats["acc"].append(acc)
        stats["time"].append(end-start)

    print("Mean values of", niter, "iterations:")
    print("Propagation", "> Accuracy:",
          str(round(np.mean(stats["acc"])*100, 2)),
          "% | Took:", sec_to_time(np.mean(stats["time"])))
