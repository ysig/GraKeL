"""The main graph kernel class, implemented as a sci-kit transformer."""
import copy
import time
import warnings

import numpy as np

from concurrent.futures import ThreadPoolExecutor
from scipy.linalg import svd
from sklearn.base import BaseEstimator
from sklearn.base import TransformerMixin
from sklearn.utils.validation import check_is_fitted

from grakel.kernels import graphlet_sampling
from grakel.kernels import random_walk
from grakel.kernels import subtree_wl
from grakel.kernels import shortest_path
from grakel.kernels import shortest_path_attr
from grakel.kernels import weisfeiler_lehman

np.random.seed(int(time.time()))

global supported_base_kernels, supported_general_kernels, default_n_components

valid_parameters = {"kernel",
                    "Nystroem",
                    "concurrency",
                    "normalize",
                    "verbose"}

supported_base_kernels = [
    "subtree_wl", "random_walk",
    "shortest_path", "subtree_rg",
    "graphlet_sampling", "subgraph_matching",
    "multiscale_laplacian",
    "lovasz_theta", "svm_theta",
    "neighborhood_hash", "neighborhood_subgraph_pairwise_distance",
    "odd_sth", "propagation",
    "pyramid_match", "jsm"
    ]

supported_general_kernels = [
    "weisfeiler_lehman",
    "hadamard_code"
    ]

default_verbose_value = True

default_random_seed_value = 5344973

default_n_components = 100


class GraphKernel(BaseEstimator, TransformerMixin):
    """A decorator for graph kernels.

    Parameters
    ----------
    kernel : list(dict(key:str, value:value))
         a list of dictionaries, or a single dicitonary that have the following
         structure:
               * "name" : [str] - with the kernel name
               * "name_of_parameter_1" : value
               * "name_of_parameter_2" : value
               *                          ...
               * "name_of_parameter_k" : value

               available "names" / "parametres" are:
                   1. base_kernels (the structure must always reach a
                   base kernel)
                        - "dirac"
                            *No arguments*
                        - "random_walk"
                            + (**o**) "lamda" : [float] < 1
                            + (**o**) method_type : [str] "fast", "baseline"

                        - "shortest_path"
                            + (**o**) "algorithm_type" : [str] "dijkstra",
                            "floyd_warshall"
                            + (**o**) "as_attributes" : [bool]
                            + (**o**) "attribute_kernel" : [function] :
                            (attribute_x, attribute_y) -> number
                            + (**o**) "with_labels" : [bool]

                        - "subtree_rg"
                            + (**o**) "h" : [int]

                        - "graphlets_sampling"
                            + (**o**) "k" : [int]
                            + (**o**) "delta" : [float]
                            + (**o**) "epsilon" : [float]
                            + (**o**) "a" : [int]
                            + (**o**) "n_samples" : [int]

                        - "multiscale_laplacian"
                            + (**o**) "L" : [int]
                            + (**o**) "gamma"

                        - "subgraph_matching"
                            + (**o**) "kv" : [function] :
                            (node_x, node_y, Lx, Ly) -> number
                            + (**o**) "ke" : [function] :
                            (edge_x, edge_y, Lx, Ly) -> number
                            + (**o**) "lw" :
                            a lambda weight function for cliques: set -> number

                        - "lovasz_theta"
                            + (**o**) "n_samples" : [int] > 1
                            + (**o**) "subsets_size_range" :
                            [touple] of two [int]
                            + (**o**) "metric" :
                            [function] (number, number) -> number

                        - "svm_theta"
                            + (**o**) "n_samples" : [int] > 1
                            + (**o**) "subsets_size_range" :
                            [touple] with 2 [int] elements
                            + (**o**) "metric" :
                            [function] (number, number) -> number

                        - "neighborhood_hash"
                            + (**o**) "nh_type" :
                            [str] "simple" or "count-sensitive"
                            + (**o**) "R" : [int] > 0
                            + (**o**) "bytes" : [int] > 0

                        - "neighborhood_subgraph_pairwise_distance"
                            + (**o**) "r" : (int) positive integer
                            + (**o**) "d" : (int) positive integer

                        - "odd_sth"
                            + (**o**) "h" : [int] > 0

                        - "propagation"
                            + (**o**) t_max: [int] > 0
                            + (**o**) T: [dict] [int]: [np.arrays]
                            + (**o**) M: [str] "H", "L1", "L2", "TV"
                            + (**o**) w: [int] > 0
                            + (**o**) base_kernel: [function]
                            x:[list[int]] , y:[list[int]] -> [number]

                        - "pyramid_match"
                            + (**o**) with_labels: [bool]
                            + (**o**) d: [int] > 0
                            + (**o**) L: [int] >= 0

                        - "jsm"
                            + (**o**) M: [int] > 0
                            + (**o**) h: [int] > 0

                   2. general_kernels (this kernel will use the next kernel
                   on the list as base kernel)
                        - "weisfeiler_lehman"
                            + (**o**) "niter" : [int]
                        - "hadamard_code"
                            + (**o**) "niter" : [int]
                            + (**o**) "hc_type" : [str] "simple", "shortened"
                            + (**o**) "rho" : [int] > 0 or -1
                            + (**o**) "L" : int, condition_of_appearance:
                            hc_type=="shortened", default=4
               where (**o**): stands for optional parameters

    Nystroem : int, optional
        Defines the number of nystroem components.
        To initialize the default (100 components), set -1 or 0.

    concurrency : int, optional
        Defines the number of workers for kernel matrix calculation using
        concurrency. To intialise the default (all possible), set -1 or 0.

    normalize : bool, optional
        Normalize the output of the graph kernel.
        Ignored when Nystroem GraphKernel object is instanciated.

    verbose : bool, optional
        Define if messages will be printed on stdout.

    random_seed : int, optional
        Initialise can provide a randomness by providing a random seed.

    Attributes
    ----------
    kernel : function
        The full kernel applied between graph objects.

    nystroem : int
        Holds the nystroem, number of components.
        If not initialised, it stands as a False
        boolean variable.

    normalize : bool
        Normalize the output of the graph kernel.

    components_ : array, shape (n_components, n_features)
        Subset of training graphs used to construct the feature map.

    nystroem_normalization_ : array, shape=(n_components, n_components)
        Normalization matrix needed for embedding.
        Square root of the kernel matrix on ``components_``.

    component_indices_ : array, shape=(n_components)
        Indices of ``components_`` in the training set.

    _verbose : bool, default=False
        Print messages in on stdout.

    _pairwise_kernel_executor : function
        The final executor destined for a single kernel calculation.

    _concurrent_executor : ThreadPoolExecutor
        An executor for applying concurrency, for the fast pairwise
        computations of the kernel matrix calculation.

    _global_random_seed : int
        The seed applying randomness.

    """

    _nystroem = False
    _normalize = False

    def __init__(self, **kargs):
        """`__init__` for `GraphKernel` object."""
        self._verbose = kargs.get("verbose", default_verbose_value)

        if "Nystroem" in kargs:
            if type(kargs["Nystroem"]) is not int:
                raise ValueError('nystroem parameter must be an int, \
                                  indicating the number of components')
            elif kargs["Nystroem"] in [0, -1]:
                # picking default number of components
                self._nystroem = default_n_components
            elif kargs["Nystroem"] <= 0:
                raise ValueError('number of nystroem components \
                                  must be positive')
            else:
                self._nystroem = kargs["Nystroem"]

        if "concurrency" in kargs:
            if type(kargs["concurrency"]) is not int:
                raise ValueError('concurrency parameter must be an int, \
                                  indicating the number of components')
            elif kargs["concurrency"] in [0, -1]:
                # Initialise an executor
                self._concurrent_executor = ThreadPoolExecutor()
            elif kargs["concurrency"] <= 0:
                raise ValueError('number of concurrency workers \
                                  must be positive')
            else:
                self._concurrent_executor = ThreadPoolExecutor(
                    max_workers=kargs["concurrency"])
            self._pairwise_kernel_executor = \
                lambda fn, *eargs, **ekargs: self._concurrent_executor.submit(
                    fn, *eargs, **ekargs).result()
        else:
            self._pairwise_kernel_executor = \
                lambda fn, *eargs, **ekargs: fn(*eargs, **ekargs)

        self._normalize = kargs.get("normalize", False)

        self._global_random_seed = kargs.get("random_seed",
                                             default_random_seed_value)

        if "kernel" in kargs:
            hidden_args = {"verbose": self._verbose,
                           "normalize": self._normalize,
                           "executor": self._pairwise_kernel_executor}
            if (type(kargs["kernel"]) is dict):
                # allow single kernel dictionary inputs
                kernel, params = self._make_kernel(
                    [kargs["kernel"]], hidden_args)

            elif (type(kargs["kernel"]) is list):
                kernel, params = self._make_kernel(
                    copy.deepcopy(kargs["kernel"]), hidden_args)
            else:
                raise ValueError('unsupported kernel format')
            self._kernel = kernel(**params)
        else:
            raise ValueError('kernel must be defined at the __init__ function \
                              of the graph kernel decorator')

        unrecognised_args = set(kargs.keys()) - valid_parameters

        if len(unrecognised_args) > 0:
            warnings.warn(
                'Ignoring unrecognised arguments:',
                ', '.join('"' + str(arg) + '"' for arg in unrecognised_args))

    def _make_kernel(self, kernel_list, hidden_args):
            """Produce the desired kernel function.

            Parameters
            ----------
            kernel_list: (list)
                List of kernel dictionaries as defined at the documentation
                of class parameters.

            Returns
            -------
            function.
            Returns the kernel, as a function of two arguments.

            """
            # If nesting type:
            kernel = kernel_list.pop(0)
            if type(kernel) is not dict:
                raise ValueError('each element of the list of kernels must' +
                                 ' be a dictionary')
            if "name" not in kernel:
                raise ValueError('each dictionary concerning a kernel must' +
                                 ' have a "name" parameter designating the' +
                                 'kernel')
            kernel_name = kernel.pop("name")
            for (keys, val) in hidden_args.items():
                kernel[keys] = val
            if kernel_name in supported_base_kernels:
                if (len(kernel_list) != 0):
                    warnings.warn('rest kernel arguments are being ignored\
                                   - reached base kernel')
                if kernel_name == "subtree_wl":
                    return (subtree_wl, kernel)
                elif kernel_name == "random_walk":
                    return (random_walk, kernel)
                elif kernel_name == "shortest_path":
                    if kernel.pop("as_attributes", False):
                        return (shortest_path_attr, kernel)
                    else:
                        return (shortest_path, kernel)
                elif kernel_name == "subtree_rg":
                    raise ValueError('still developing')
                elif kernel_name == "graphlet_sampling":
                    if "random_seed" not in kernel and \
                        self._global_random_seed is not \
                            default_random_seed_value:
                            kernel["random_seed"] = self._global_random_seed
                    return (graphlet_sampling, kernel)
                elif kernel_name == "multiscale_laplacian":
                    raise ValueError('still developing')
                elif kernel_name == "subgraph_matching":
                    raise ValueError('still developing')
                elif kernel_name == "lovasz_theta":
                    raise ValueError('still developing')
                elif kernel_name == "svm_theta":
                    raise ValueError('still developing')
                elif kernel_name == "neighborhood_hash":
                    raise ValueError('still developing')
                elif kernel_name == "neighborhood_subgraph_pairwise_distance":
                    raise ValueError('still developing')
                elif kernel_name == "odd_sth":
                    raise ValueError('still developing')
                elif kernel_name == "propagation":
                    raise ValueError('still developing')
                elif kernel_name == "pyramid_match":
                    raise ValueError('still developing')
                elif kernel_name == "jsm":
                    raise ValueError('still developing')
            elif kernel_name in supported_general_kernels:
                if (len(kernel_list) == 0):
                    raise ValueError(str(kernel_name)+' is not a base kernel')
                else:
                    kernel["base_kernel"] = self._make_kernel(kernel_list, {})
                if kernel_name == "weisfeiler_lehman":
                    return (weisfeiler_lehman, kernel)
                if kernel_name == "hadamard_code":
                    raise ValueError('still developing')
            else:
                raise ValueError('unsupported kernel: ' + str(kernel_name))

    def fit(self, X, y=None):
        """Fit a dataset, for a transformer.

        Parameters
        ----------
        X : iterable
            Each element must be an iterable with at most three features and at
            least one. The first that is obligatory is a valid graph structure
            (adjacency matrix or edge_dictionary) while the second is
            node_labels and the third edge_labels (that fitting the given grap
            format). The train samples.

        y : None
            There is no need of a target in a transformer, yet the pipeline API
            requires this parameter.

        Returns
        -------
        self : object
        Returns self.

        """
        # Input validation and parsing
        if bool(self._nystroem):
            nx = len(X)
            # get basis vectors
            if self._nystroem > nx:
                n_components = nx
                warnings.warn("n_components > n_samples. This is not \
                              possible.\nn_components was set to n_samples, \
                              which results in inefficient evaluation of the \
                              full kernel.")
            else:
                n_components = self._nystroem

            n_components = min(nx, n_components)
            inds = np.random.permutation(nx)
            basis_inds = inds[:n_components]
            basis = [X[i] for i in range(basis_inds)]

            self.components_ = basis
            self._nystroem = n_components

            # sqrt of kernel matrix on basis vectors
            U, S, V = svd(self._kernel.fit_transform(basis))
            S = np.maximum(S, 1e-12)
            self._nystroem = n_components
            self.nystroem_normalization_ = np.dot(U / np.sqrt(S), V)
            self.components_ = basis
            self.component_indices_ = inds
        else:
            self._kernel.fit(X)

        # Return the transformer
        return self

    def transform(self, X):
        """Calculate the kernel matrix, between given and fitted dataset.

        Paramaters
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
        # Check if nystroem has been initialised had been called
        if bool(self._nystroem):
            check_is_fitted(self, 'components_')

        # Transform - calculate kernel matrix
        if bool(self._nystroem):
            return np.dot(self._kernel.transform(X),
                          self.nystroem_normalization_.T)
        else:
            return self._kernel.transform(X)

    def fit_transform(self, X):
        """Fit and transform, on the same dataset.

        Paramaters
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
        # Transform - calculate kernel matrix
        if bool(self._nystroem):
            self.fit(X)
            return np.dot(self._kernel.transform(X),
                          self.nystroem_normalization_.T)
        else:
            return self._kernel.fit_transform(X)
