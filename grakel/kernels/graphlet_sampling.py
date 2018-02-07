"""The graphlet sampling kernel :cite:`Shervashidze2009EfficientGK`."""
import collections
import itertools
import math
import random
import warnings

import numpy as np

import pynauty

from scipy.interpolate import interp1d
from sklearn.exceptions import NotFittedError
from sklearn.utils.validation import check_is_fitted

from grakel.graph import Graph
from grakel.tools import matrix_to_dict
from grakel.kernels import kernel

default_random_seed_value = 15487103


class graphlet_sampling(kernel):
    r"""The graphlet sampling kernel.

    See :cite:`Shervashidze2009EfficientGK`.

    If either "delta", "epsilon", "a" or "n_samples" is given calculates
    the kernel value for the given (or derived) random picked n_samples, by
    randomly sampling from k from 3 to 5.
    Otherwise calculates the kernel value drawing all possible connected
    samples of size k.

    Parameters
    ----------
    random_seed : int, default=15487103

    k : int, default=5
        The dimension of the given graphlets.

    delta : float, default=0.05
        Confidence level (typically 0.05 or 0.1).
        For calculation of the number of samples achieving the certain bound.
        n_samples argument must not be provided and for initialising the
        default value either "epsilon" or "a" must be set.

    epsilon : float, default=0.05
        Precision level (typically 0.05 or 0.1).
        For calculation of the number of samples achieving the certain bound.
        n_samples argument must not be provided and for initialising the
        default value either "delta" or "a" must be set.

    a : int
        Number of isomorphism classes of graphlets.
        If -1 the number is the maximum possible, from a database 1 until 9
        or else predicted through interpolation.
        For calculation of the number of samples achieving the certain bound.
        n_samples argument must not be provided and for initialising the
        default value either "delta" or "epsilon" must be set.

    n_samples : int
        Sets the value of randomly drawn random samples,
        from sizes between 3..k.

    Attributes
    ----------
    X : dict
        A dictionary of pairs between each input graph and a bins where the
        sampled graphlets have fallen.

    _sample_graphlets : function
        A function taking as input a binary adjacency matrix, parametrised
        to work for the certain samples, k and deterministic/propabilistic
        mode.

    _graph_bins : dict
        A dictionary of graph bins holding pynauty objects

    _nx : int
        Holds the number of sampled X graphs.

    _ny : int
        Holds the number of sampled Y graphs.

    _X_diag : np.array, shape=(_nx, 1)
        Holds the diagonal of X kernel matrix in a numpy array, if calculated
        (`fit_transform`).

    _phi_X : np.array, shape=(_nx, len(_graph_bins))
        Holds the features of X in a numpy array, if calculated.
        (`fit_transform`).

    """

    _graph_format = "adjacency"
    _graph_bins = dict()

    def __init__(self, **kargs):
        """Initialise a subtree_wl kernel."""
        self._valid_parameters |= {"random_seed",
                                   "k", "delta", "epsilon", "a", "n_samples"}
        super(graphlet_sampling, self).__init__(**kargs)

        random.seed(int(kargs.get("random_seed", default_random_seed_value)))

        k = kargs.get("k", 5)

        if k > 10:
            warnings.warn('graphlets are too big - computation may be slow')
        elif k < 3:
            raise ValueError('k must be bigger than 3')

        if "n_samples" in kargs:
            # Get the numbr of samples
            n_samples = kargs["n_samples"]

            # Display a warning if arguments ignored
            args = [arg for arg in ["delta", "epsilon", "a"] if arg in kargs]
            if len(args):
                warnings.warn('Number of samples defined as input, ' +
                              'ignoring arguments:', ', '.join(args))

            # Initialise the sample graphlets function
            self._sample_graphlets = lambda A: \
                sample_graphlets_probabilistic(A, k, n_samples)
        elif "delta" in kargs or "epsilon" in kargs or "a" in kargs:
            # Otherwise if delta exists
            delta = kargs.get("delta", 0.05)
            # or epsilon
            epsilon = kargs.get("epsilon", 0.05)
            # or a
            a = kargs.get("a", -1)

            # check the fit constraints
            if delta > 1 or delta < 0:
                raise ValueError('delta must be in the range (0,1)')

            if epsilon > 1 or epsilon < 0:
                raise ValueError('epsilon must be in the range (0,1)')

            if type(a) is not int:
                raise ValueError('a must be an integer')
            elif a == 0:
                raise ValueError('a cannot be zero')
            elif a < -1:
                raise ValueError('negative a smaller than -1 have no meaning')

            if(a == -1):
                if(k > 9):
                    warnings.warn(
                        'warning for such size number of isomorphisms is not' +
                        ' known - interpolation on know values will be used')
                    # Use interpolations
                    fallback_map = {1: 1, 2: 2, 3: 4, 4: 8, 5: 19, 6: 53, 7:
                                    209, 8: 1253, 9: 13599}

                    isomorphism_prediction = \
                        interp1d(list(fallback_map.keys()),
                                 list(fallback_map.values()), kind='cubic')
                    a = isomorphism_prediction(k)
                else:
                    a = fallback_map[k]

            # and calculate number of samples
            n_samples = math.ceil(2*(a*np.log10(2) +
                                  np.log10(1/delta))/(epsilon**2))

            self._sample_graphlets = lambda A: \
                sample_graphlets_probabilistic(A, k, n_samples)
        else:
            self._sample_graphlets = lambda A: \
                sample_graphlets_all_connected(A, k)

    def transform(self, X):
        """Calculate the kernel matrix, between given and fitted dataset.

        Parameters
        ----------
        X : iterable
            Each element must be an iterable with at most three features and at
            least one. The first that is obligatory is a valid graph structure
            (adjacency matrix or edge_dictionary) while the second is
            node_labels and the third edge_labels (that fitting the given graph
            format). If None the kernel matrix is calculated upon fit data.
            The test samples.

        Returns
        -------
        K : numpy array, shape = [n_targets, n_input_graphs]
            corresponding to the kernel matrix, a calculation between
            all pairs of graphs between target an features

        """
        self._method_calling = 3
        # Check is fit had been called
        check_is_fitted(self, ['X'])

        # Input validation and parsing
        if X is None:
            raise ValueError('transform input cannot be None')
        else:
            Y = self.parse_input(X)

        # Transform - calculate kernel matrix
        try:
            check_is_fitted(self, ['phi_X'])
            phi_x = self._phi_X
        except NotFittedError:
            phi_x = np.zeros(shape=(self._nx, len(self._graph_bins)))
            for ((i, j), v) in self.X.items():
                phi_x[i, j] = v
            self._phi_X = phi_x
        phi_y = np.zeros(shape=(self._ny, len(self._graph_bins) +
                                len(self._Y_graph_bins)))
        for ((i, j), v) in Y.items():
            phi_y[i, j] = v

        # store _phi_Y for independent (of normalization arg diagonal-calls)
        self._phi_Y = phi_y
        km = np.dot(phi_y[:, :len(self._graph_bins)], phi_x.T)
        if self._normalize:
            X_diag, Y_diag = self.diagonal()
            km /= np.sqrt(np.outer(Y_diag, X_diag))
        return km

    def fit_transform(self, X):
        """Fit and transform, on the same dataset.

        Parameters
        ----------
        X : iterable
            Each element must be an iterable with at most three features and at
            least one. The first that is obligatory is a valid graph structure
            (adjacency matrix or edge_dictionary) while the second is
            node_labels and the third edge_labels (that fitting the given graph
            format). If None the kernel matrix is calculated upon fit data.
            The test samples.

        Returns
        -------
        K : numpy array, shape = [n_targets, n_input_graphs]
            corresponding to the kernel matrix, a calculation between
            all pairs of graphs between target an features

        """
        self._method_calling = 2
        self.fit(X)

        # calculate feature matrices.
        phi_x = np.zeros(shape=(self._nx, len(self._graph_bins)))
        for ((i, j), v) in self.X.items():
            phi_x[i, j] = v

        # Transform - calculate kernel matrix
        self._phi_X = phi_x
        km = np.dot(phi_x, phi_x.T)

        self._X_diag = np.diagonal(km).reshape(km.shape[0], 1)
        if self._normalize:
            return np.divide(km, np.sqrt(np.outer(self._X_diag, self._X_diag)))
        else:
            return km

    def diagonal(self):
        """Calculate the kernel matrix diagonal for fitted data.

        A funtion called on transform on a seperate dataset to apply
        normalization on the exterior.

        Parameters
        ----------
        None.

        Returns
        -------
        X_diag : np.array
            The diagonal of the kernel matrix, of the fitted data.
            This consists of kernel calculation for each element with itself.

        Y_diag : np.array
            The diagonal of the kernel matrix, of the transformed data.
            This consists of kernel calculation for each element with itself.

        """
        # Check is fit had been called
        check_is_fitted(self, ['_phi_X', '_phi_Y'])
        try:
            check_is_fitted(self, ['X_diag'])
        except NotFittedError:
            # Calculate diagonal of X
            self._X_diag = np.sum(np.square(self._phi_X), axis=1)
            self._X_diag = np.reshape(self._X_diag, (self._X_diag.shape[0], 1))
        # Calculate diagonal of Y
        Y_diag = np.sum(np.square(self._phi_Y), axis=1)

        return self._X_diag, np.reshape(Y_diag, (Y_diag.shape[0], 1))

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
        local_values : dict
            A dictionary of pairs between each input graph and a bins where the
            sampled graphlets have fallen.

        """
        if not isinstance(X, collections.Iterable):
            raise ValueError('input must be an iterable\n')
        else:
            i = -1
            if self._method_calling == 1:
                self._graph_bins = dict()
            elif self._method_calling == 3:
                self._Y_graph_bins = dict()
            local_values = dict()
            for (idx, x) in enumerate(iter(X)):
                if type(x) is Graph:
                    A = x.get_adjacency_matrix()
                elif isinstance(x, collections.Iterable) and \
                        len(x) in [0, 1, 2, 3]:
                    if len(x) == 0:
                        warnings.warn('Ignoring empty element on ' +
                                      'index: '+str(idx))
                        continue
                    else:
                        A = Graph(x[0], {}, {},
                                  self._graph_format).get_adjacency_matrix()
                else:
                    raise ValueError('each element of X must be either a ' +
                                     'graph or an iterable with at least 1 ' +
                                     'and at most 3 elements\n')
                A = (A > 0).astype(int)
                i += 1
                # sample graphlets based on the initialized method
                samples = self._sample_graphlets(A)

                if self._method_calling == 1:
                    for (j, sg) in enumerate(samples):
                        # add the graph to an isomorphism class
                        if len(self._graph_bins) == 0:
                            self._graph_bins[0] = sg
                            local_values[(i, 0)] = 1
                        else:
                            newbin = True
                            for j in range(len(self._graph_bins)):
                                if pynauty.isomorphic(self._graph_bins[j], sg):
                                    newbin = False
                                    if (i, j) not in local_values:
                                        local_values[(i, j)] = 1
                                    local_values[(i, j)] += 1
                                    break
                            if newbin:
                                local_values[(i, len(self._graph_bins))] = 1
                                self._graph_bins[len(self._graph_bins)] = sg
                elif self._method_calling == 3:
                    for (j, sg) in enumerate(samples):
                        # add the graph to an isomorphism class
                        newbin = True
                        for j in range(len(self._graph_bins)):
                            if pynauty.isomorphic(self._graph_bins[j], sg):
                                newbin = False
                                if (i, j) not in local_values:
                                    local_values[(i, j)] = 1
                                local_values[(i, j)] += 1
                                break
                        if newbin:
                            if len(self._Y_graph_bins) == 0:
                                self._Y_graph_bins[0] = sg
                                local_values[(i,
                                              len(self._graph_bins))] = 1
                            else:
                                newbin_Y = True
                                start = len(self._graph_bins)
                                for j in range(len(self._Y_graph_bins)):
                                    if pynauty.isomorphic(
                                            self._Y_graph_bins[j], sg):
                                        newbin_Y = False
                                        bin_key = (i, j + start)
                                        if bin_key not in local_values:
                                            local_values[bin_key] = 1
                                        local_values[bin_key] += 1
                                        break
                                if newbin_Y:
                                    idx = start + len(self._Y_graph_bins)
                                    local_values[(i, idx)] = 1
                                    self._Y_graph_bins[idx] = sg

            if i == -1:
                raise ValueError('parsed input is empty')

            if self._method_calling == 1:
                self._nx = i+1
            elif self._method_calling == 3:
                self._ny = i+1
            return local_values


def sample_graphlets_probabilistic(A, k, n_samples):
    """Propabilistical sampling of n_samples of 3..k sized graphs.

    Parameters
    ----------
    A : np.array
        A binary array defining a certain graph.

    k : int
        The maximum dimension of the sampled graphlets.

    n_samples : int
        Sets the value of randomly drawn random samples,
        from sizes between 3..k

    Returns
    -------
    graphlets : generator
        Returns a generator of sampled graphlets (as pynauty graphs),
        from sizes between 3..k.

    """
    s = list(range(A.shape[0]))
    min_r, max_r = min(3, A.shape[0]), min(k, A.shape[0])
    if min_r == max_r:
        rsamp = lambda *args: min_r
    else:
        rsamp = lambda *args: random.randint(min_r, max_r)

    to_edge_dict_binary = lambda x, k: matrix_to_dict(x, '==', 1, False)
    for i in range(n_samples):
        index_rand = random.sample(s, rsamp())
        Q = A[index_rand, :][:, index_rand]
        yield pynauty.Graph(Q.shape[0], True, to_edge_dict_binary(Q, k))


def sample_graphlets_all_connected(A, k):
    """All the connected graphlets of size k of a given graph.

    Parameters
    ----------
    A : np.array
        A binary array defining a certain graph.

    k : int
        The maximum dimension of the sampled graphlets.

    Returns
    -------
    graphlets : generator
        Returns a generator of sampled graphlets (as pynauty graphs),
        of size k.

    """
    to_edge_dict_binary = lambda x: matrix_to_dict(x, '==', 1, False)
    for i in itertools.permutations(range(A.shape[0]), min(k, A.shape[0])):
        Q = A[i, :][:, i]
        if 0 not in np.sum(Q, axis=1):
            yield pynauty.Graph(Q.shape[0], True, to_edge_dict_binary(Q))
