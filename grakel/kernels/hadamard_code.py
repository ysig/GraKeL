"""The hadamard code kernel as defined in :cite:`icpram16`."""
import collections
import math
import warnings

import numpy as np

from scipy.linalg import hadamard

from sklearn.exceptions import NotFittedError
from sklearn.utils.validation import check_is_fitted

from grakel.graph import Graph
from grakel.kernels import kernel


class hadamard_code(kernel):
    """Hadamard code kernel, as proposed in :cite:`icpram16`.

    Parameters
    ----------
    base_kernel : `grakel.kernels.kernel` or tuple
        If tuple it must consist of a valid kernel object and a
        dictionary of parameters. General parameters concerning
        normalization, concurrency, .. will be ignored, and the
        ones of given on `__init__` will be passed in case it is needed.

    hc_type : str, valid_inputs={"simple", "shortened"}, default="simple"
        The hadamard code kernel type as defined in :cite:`icpram16`.

    rho : int, condition_of_appearance: hc_type=="shortened", default=-1
        The size of each single bit arrays. If -1 is chosen r is calculated as
        the biggest possible that satisfies an equal division.

    L : int, condition_of_appearance: hc_type=="shortened", default=4
        The number of bytes to store the bitarray of each label.

    niter : int, default=5
        The number of iterations.

    Attributes
    ----------
    _niter : int
        The number of iterations.

    _shortened : bool
        A flag signifying if the shortened version of the algorithm
        will be used.

    _add : function
        A function setted relevant to the version of the algorithm
        for adding hashed labels.

    _get : function
        A function setted relevant to the version of the algorithm
        for getting a label element.

    _base_kernel : function
        A void function that initializes a base kernel object.

    """

    _graph_format = "auto"

    def __init__(self, **kargs):
        """Initialise a `hadamard_code` kernel."""
        base_params = self._valid_parameters.copy()
        self._valid_parameters |= {"base_kernel",
                                   "niter",
                                   "hc_type",
                                   "rho",
                                   "L"}
        super(hadamard_code, self).__init__(**kargs)

        self._niter = kargs.get("niter", 5)
        if self._niter < 1:
            raise ValueError('niter must be an integer bigger than zero')

        hc_type = kargs.get("hc_type", 'simple')
        if hc_type == 'simple':
            self._shortened = False
            self._add = np.add
            self._get = lambda A, i: A[i, :]
        elif hc_type == 'shortened':
            # extract rho
            self._rho = kargs.get('rho', -1)
            if self._rho <= 0 and self._rho != -1:
                raise ValueError('rho must be bigger than zero or equal to 1')

            # extract L
            self._L = kargs.get('L', 4)
            if self._L <= 0:
                raise ValueError('L must be a positive integer as it ' +
                                 'corresponds to the number of bytes')

            self._L = self._L * 8

            # initialise addition labels and get for an internal use
            self._shortened = True
            self._add = lambda x, y: x + y
            self._get = lambda A, i: A[i]
        else:
            raise ValueError('unrecognised hadamard code kernel type')

        if "base_kernel" not in kargs:
            raise ValueError('User must provide a base kernel.')
        else:
            if type(kargs["base_kernel"]) is type and \
                    issubclass(kargs["base_kernel"], kernel):
                base_kernel = kargs["base_kernel"]
                params = dict()
            else:
                try:
                    base_kernel, params = kargs["base_kernel"]
                except Exception:
                    raise ValueError('Base kernel was not provided in the ' +
                                     'correct way. Check documentation.')

                if not (type(base_kernel) is type and
                        issubclass(base_kernel, kernel)):
                    raise ValueError('The first argument must be a valid ' +
                                     'grakel.kernel.kernel Object')
                if type(params) is not dict:
                    raise ValueError('If the second argument of base ' +
                                     'kernel exists, it must be a diction' +
                                     'ary between parameters names and values')
                params.pop("normalize", None)
                for p in base_params:
                        params.pop(p, None)

            params["normalize"] = False
            params["verbose"] = self._verbose
            params["executor"] = self._executor
            self._base_kernel = lambda *args: base_kernel(**params)

    def parse_input(self, X):
        """Parse input for weisfeiler lehman.

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
        base_kernel : object
            Returns base_kernel. Only if called from `fit` or `fit_transform`.

        K : np.array
            Returns the kernel matrix. Only if called from `transform` or
            `fit_transform`.

        """
        # Input validation and parsing
        if not isinstance(X, collections.Iterable):
            raise ValueError('input must be an iterable\n')
        else:
            nx = 0
            if self._method_calling in [1, 2]:
                nl = 0
                labels_enum = dict()
                base_kernel = dict()
                for kidx in range(self._niter):
                    base_kernel[kidx] = self._base_kernel()
            elif self._method_calling == 3:
                nl = len(self._labels_enum)
                labels_enum = dict(self._labels_enum)
                base_kernel = self.X
            inp = list()
            neighbors = list()
            for (idx, x) in enumerate(iter(X)):
                is_iter = False
                if isinstance(x, collections.Iterable):
                    x, is_iter = list(x), True
                if is_iter and len(x) in [0, 2, 3]:
                    if len(x) == 0:
                        warnings.warn('Ignoring empty element on index: '
                                      + str(idx))
                        continue
                    else:
                        x = Graph(x[0], x[1], {},
                                  graph_format=self._graph_format)
                elif type(x) is not Graph:
                    raise ValueError('each element of X must be either a ' +
                                     'graph object or a list with at least ' +
                                     'a graph like object and node labels ' +
                                     'dict \n')

                label = x.get_labels(purpose='any')
                inp.append((x.get_graph_object(), label))
                neighbors.append(x.get_edge_dictionary())
                for v in set(label.values()):
                    if v not in labels_enum:
                        labels_enum[v] = nl
                        nl += 1
                nx += 1
            if nx == 0:
                raise ValueError('parsed input is empty')

        # Calculate the hadamard matrix
        ord_H = int(2**(math.ceil(math.log2(nl))))
        H = hadamard(ord_H)

        if self._shortened:
            # In order to apply shortened correctly all the transform
            # elements must be considered upon fit. This is not possible
            # on the fit transform model. So we keep the fit labels for fit
            # and add new on transform.
            if self._method_calling in [1, 2]:
                if self._rho == -1:
                    rho = int(self._L/(ord_H-1))
                    if rho <= 0:
                        raise ValueError('The default calculated rho is too ' +
                                         'small, raise more L ~ L must' +
                                         ' approach the number of labels (' +
                                         str(ord_H) + ')')
                elif self._rho*(ord_H-1) > self._L:
                    raise ValueError('rho*(ord_H-1)='+str(self._rho)+'*(' +
                                     str(ord_H) + '-1) > L='+str(self._L))
                else:
                    rho = self._rho

                bit_array = dict()
                bool_H = (H > 0).astype(bool)
                # construct the bit array as string concatenation
                header_size = (self._L - rho*(ord_H-1))
                for i in range(ord_H):
                    bit_array[i] = header_size * str(int(bool_H[i, -1]))
                    for j in range(ord_H-1, -1, -1):
                        bit_array[i] = bit_array[i] + \
                            rho*str(int(bool_H[i, j]))
                    bit_array[i] = int('0b'+bit_array[i], 2)
                Labeling = bit_array

                # Store information for handling of transform
                self._fit_bit_array_ = bit_array
                self._fit_ord_H_ = ord_H
                self._fit_rho_ = rho
            elif self._method_calling == 3:
                if self._fit_ord_H_ < ord_H:
                    # Length *L* is being ommited by the new elements
                    # Length *L* stays valid for elements of fit.
                    rho = self._fit_rho_
                    bit_array = dict()
                    bool_H = (H > 0).astype(bool)
                    header_size = (self._L - rho*(self._fit_ord_H_-1))
                    for i in range(self._fit_ord_H_, ord_H):
                        # All the first elements keep fit_rho steady
                        for j in range(ord_H, self._fit_ord_H_, -1):
                            bit_array[i] = bit_array[i] + \
                                rho*str(int(bool_H[i, j]))

                        # continue normally first element of fit as in fit
                        bit_array[i] = bit_array[i] + \
                            header_size * str(int(bool_H[i, -1]))

                        # all continuing elements the same
                        for j in range(self._fit_ord_H_-1, -1, -1):
                            bit_array[i] = bit_array[i] + \
                                rho*str(int(bool_H[i, j]))
                        bit_array[i] = int('0b'+bit_array[i], 2)

                    Labeling = dict(self._fit_bit_array_, **bit_array)
                else:
                    Labeling = self._fit_bit_array_
        else:
            Labeling = H

        # Intial labeling of vertices based on their corresponding Hadamard
        # code (i-th row of the Hadamard matrix) where i is the i-th label on
        # enumeration
        new_graphs = list()
        for (obj, label) in inp:
            new_labels = dict()
            for (k, v) in label.items():
                new_labels[k] = self._get(Labeling, labels_enum[v])
            new_graphs.append((obj, label))

        # Add the zero iteration element
        if self._method_calling == 1:
            base_kernel[0].fit(new_graphs)
        elif self._method_calling == 2:
            K = base_kernel[0].fit_transform(new_graphs)
        elif self._method_calling == 3:
            K = base_kernel[0].transform(new_graphs)

        # Main
        for i in range(1, self._niter):
            graphs = new_graphs
            new_graphs = list()
            for ((obj, old_labels), neighbor) in zip(graphs, neighbors):
                # Find unique labels and sort them for both graphs
                # Keep for each node the temporary
                new_labels = dict()
                for (k, ns) in neighbor.items():
                    new_labels[k] = old_labels[k]
                    for q in ns:
                        new_labels[k] = self._add(new_labels[k], old_labels[q])
                new_graphs.append((obj, new_labels))

            # calculate kernel
            if self._method_calling == 1:
                base_kernel[i].fit(new_graphs)
            elif self._method_calling == 2:
                K += base_kernel[i].fit_transform(new_graphs)
            elif self._method_calling == 3:
                K += base_kernel[i].transform(new_graphs)

        if self._method_calling == 1:
            self._labels_enum = labels_enum
            return base_kernel
        elif self._method_calling == 2:
            self._labels_enum = labels_enum
            return K, base_kernel
        elif self._method_calling == 3:
            return K

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
            km = self.parse_input(X)

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
        if X is None:
            raise ValueError('transform input cannot be None')
        else:
            km, self.X = self.parse_input(X)

        self._X_diag = np.reshape(np.diagonal(km), (km.shape[0], 1))
        if self._normalize:
            return np.divide(km,
                             np.sqrt(np.outer(self._X_diag, self._X_diag)))
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
        # Check if fit had been called
        check_is_fitted(self, ['X'])
        try:
            check_is_fitted(self, ['_X_diag'])
            Y_diag = self.X[0].diagonal()[1]
            for i in range(1, self._niter):
                Y_diag += self.X[i].diagonal()[1]
        except NotFittedError:
            # Calculate diagonal of X
            X_diag, Y_diag = self.X[0].diagonal()
            X_diag.flags.writeable = True
            for i in range(1, self._niter):
                x, y = self.X[i].diagonal()
                X_diag += x
                Y_diag += y
            self._X_diag = X_diag

        return self._X_diag, Y_diag
