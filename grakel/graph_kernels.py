"""The main graph kernel class, implemented as a sci-kit transformer."""
import copy
import time
import warnings

import numpy as np

# from concurrent.futures import ThreadPoolExecutor
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
from grakel.kernels import pyramid_match
from grakel.kernels import neighborhood_hash
from grakel.kernels import subgraph_matching
from grakel.kernels import neighborhood_subgraph_pairwise_distance
from grakel.kernels import lovasz_theta
from grakel.kernels import svm_theta
from grakel.kernels import jsm
from grakel.kernels import odd_sth
from grakel.kernels import propagation
from grakel.kernels import hadamard_code
from grakel.kernels import multiscale_laplacian
from grakel.kernels import multiscale_laplacian_fast
from grakel.kernels import vertex_histogram
from grakel.kernels import edge_histogram

# Python 2/3 cross-compatibility import
from future.utils import iteritems

np.random.seed(int(time.time()))

valid_parameters = {"kernel",
                    "Nystroem",
                    "n_jobs",
                    "normalize",
                    "verbose",
                    "random_seed"}

supported_base_kernels = [
    "subtree_wl", "random_walk",
    "shortest_path", "subtree_rg",
    "graphlet_sampling", "subgraph_matching",
    "multiscale_laplacian",
    "lovasz_theta", "svm_theta",
    "neighborhood_hash", "neighborhood_subgraph_pairwise_distance",
    "odd_sth", "propagation",
    "pyramid_match", "jsm",
    "propagation", "vertex_histogram", "edge_histogram"
    ]

supported_general_kernels = [
    "weisfeiler_lehman",
    "hadamard_code"
    ]

default_verbose_value = True

default_random_seed_value = 5344973

default_n_components = 100


class GraphKernel(BaseEstimator, TransformerMixin):
    r"""A decorator for graph kernels.

    Parameters
    ----------
    kernel : list(dict(key:str, value:value))
        A list of dictionaries, or a single dicitonary that have the following
        structure:
            * "name" : [str] - with the kernel name

            * "name_of_parameter_1" : value

            * "name_of_parameter_2" : value

            * :math:`\;\cdots\;`

            * "name_of_parameter_k" : value

        available "names" / "parametres" are:
            1. base_kernels (the structure must always reach a base kernel)
                - "dirac"
                    *No arguments*

                - "random_walk"
                    + (**o**) "lambda" : float

                    + (**o**) "method_type" : [str], "baseline", "fast"

                    + (**o**) "kernel_type" : [str], "geometric",
                      "exponential"

                    + (**o**) "p" : [int] > 0


                - "shortest_path"
                    + (**o**) "algorithm_type" : [str] "dijkstra",
                      "floyd_warshall"

                    + (**o**) "as_attributes" : [bool]

                    + (**o**) "attribute_kernel" : [function] :
                      (attribute_x, attribute_y) -> number

                    + (**o**) "with_labels" : [bool]

                - "subtree_rg"
                    + (**o**) "h" : [int]

                - "graphlet_sampling"
                    + (**o**) "k" : [int]

                    + (**o**) "delta" : [float]

                    + (**o**) "epsilon" : [float]

                    + (**o**) "a" : [int]

                    + (**o**) "n_samples" : [int]

                - "multiscale_laplacian"
                    + (**o**) "which" : [str] "slow", "fast"

                    + (**o**) "L" : [int] > 0

                    + (**o**) "gamma" : [float] > .0

                    + (**o**) "heta" : [float] > .0

                    + (**o**) "N" : [int] > 0, if "which": "fast"

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

                - "edge"
                    *No arguments*

                - "vertex"
                    *No arguments*

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

    Nystroem : int or bool, optional
        Defines the number of nystroem components.
        To initialize the default (100 components), set -1 or 0.

    n_jobs : int, optional
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
    kernel_ : function
        The full kernel applied between graph objects.

    nystroem_ : int
        Holds the nystroem, number of components.
        If not initialised, it stands as a False
        boolean variable.

    components_ : array, shape=(n_components, n_features)
        Subset of training graphs used to construct the feature map.

    nystroem_normalization_ : array, shape=(n_components, n_components)
        Normalization matrix needed for embedding.
        Square root of the kernel matrix on ``components_``.

    component_indices_ : array, shape=(n_components)
        Indices of ``components_`` in the training set.

    pairwise_kernel_executor_ : function
        The final executor destined for a single kernel calculation.

    concurrent_executor_ : ThreadPoolExecutor
        An executor for applying concurrency, for the fast pairwise
        computations of the kernel matrix calculation.

    """

    initialised_ = None

    def __init__(self,
                 kernel=None,
                 normalize=False,
                 verbose=False,
                 n_jobs=0,
                 random_seed=default_random_seed_value,
                 Nystroem=False):
        """`__init__` for `GraphKernel` object."""
        self.kernel = kernel
        self.normalize = normalize
        self.verbose = verbose
        self.n_jobs = n_jobs
        self.random_seed = random_seed
        self.Nystroem = Nystroem
        self.initialised_ = {"kernel": False,
                             "Nystroem": False,
                             "n_jobs": False}

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
        # Initialise the Graph Kernel.
        self.initialise_()

        # Input validation and parsing
        if bool(self.nystroem_):
            X = list(X)
            nx = len(X)
            # get basis vectors
            if self.nystroem_ > nx:
                n_components = nx
                warnings.warn("n_components > n_samples. This is not "
                              "possible.\nn_components was set to n_samples"
                              ", which results in inefficient evaluation of"
                              " the full kernel.")
            else:
                n_components = self.nystroem_

            n_components = min(nx, n_components)
            inds = np.random.permutation(nx)
            basis_inds = inds[:n_components]
            basis = [X[i] for i in basis_inds]

            # sqrt of kernel matrix on basis vectors
            U, S, V = svd(self.kernel_.fit_transform(basis))
            S = np.maximum(S, 1e-12)
            self.nystroem_ = n_components
            self.nystroem_normalization_ = np.dot(U / np.sqrt(S), V)
            self.components_ = basis
            self.component_indices_ = inds
        else:
            self.kernel.fit(X)

        # Return the transformer
        return self

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
        # Check if nystroem has been initialised had been called
        if bool(self.nystroem_):
            check_is_fitted(self, 'components_')

        # Transform - calculate kernel matrix
        if bool(self.nystroem_):
            K = self.kernel_.transform(X).dot(self.nystroem_normalization_.T)
        else:
            K = self.kernel_.transform(X)

        if K.shape == (1, 1):
            return K[0, 0]
        else:
            return K

    def fit_transform(self, X, y=None):
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

        y : None
            There is no need of a target in a transformer, yet the pipeline API
            requires this parameter.

        Returns
        -------
        K : numpy array, shape = [n_targets, n_input_graphs]
            corresponding to the kernel matrix, a calculation between
            all pairs of graphs between target an features

        """
        # Initialise the Graph Kernel
        self.initialise_()

        # Transform - calculate kernel matrix
        if bool(self.nystroem_):
            self.fit(X)
            K = self.kernel_.transform(X).dot(self.nystroem_normalization_.T)
        else:
            K = self.kernel_.fit_transform(X)

        if K.shape == (1, 1):
            return K[0, 0]
        else:
            return K

    def initialise_(self):
        """Initialise all transformer arguments, needing initialisation."""
        if not self.initialised_["Nystroem"]:
            if type(self.Nystroem) not in [int, bool]:
                raise ValueError('Nystroem parameter must be an int, '
                                 'indicating the number of components'
                                 'or a boolean')
            elif self.Nystroem is False:
                self.nystroem_ = False
            elif self.Nystroem in [0, -1] or self.Nystroem is True:
                # picking default number of components
                self.nystroem_ = default_n_components
            elif self.Nystroem <= 0:
                    raise ValueError('number of nystroem components '
                                     'must be positive')
            else:
                self.nystroem_ = self.Nystroem
            self.initialised_["Nystroem"] = True

        if not self.initialised_["n_jobs"]:
            if self.n_jobs == 0:
                pass
                # self._pairwise_kernel_executor =
                # lambda fn, *eargs, **ekargs:\
                #    fn(*eargs, **ekargs)
            else:
                warnings.warn('feature is currently not implemented')
                if type(self.n_jobs) is not int:
                    raise ValueError('n_jobs parameter must be an int, '
                                     'indicating the number of workers')
                elif self.n_jobs == -1:
                    pass
                    # Initialise an executor
                    # self._concurrent_executor = ThreadPoolExecutor()
                elif self.n_jobs <= 0:
                    raise ValueError('number of jobs (concurrent workers) '
                                     'must be positive')
                else:
                    pass
                    # self._concurrent_executor = ThreadPoolExecutor(
                    #    max_workers=kargs["n_jobs"])
                # self._pairwise_kernel_executor = lambda fn, *eargs, **ekargs:
                # \
                # self._concurrent_executor.submit(fn, *eargs, **ekargs).
                # result()
            self.initialised_["Nystroem"] = True

            self.pairwise_kernel_executor_ = lambda fn, *eargs, **ekargs: \
                fn(*eargs, **ekargs)

        if not self.initialised_["kernel"]:
            if self.kernel is None:
                raise ValueError('kernel must be defined at the __init__ '
                                 'function of the graph kernel decorator ')
            else:
                hidden_args = {"verbose": self.verbose,
                               "normalize": self.normalize,
                               "executor": self.pairwise_kernel_executor_}

                k = self.kernel
                if type(k) is dict:
                    # allow single kernel dictionary inputs
                    k = [self.kernel]
                elif type(k) is not list:
                    raise ValueError('unsupported kernel format')

                kernel, params = self.make_kernel_(
                    copy.deepcopy(k), hidden_args)

                self.kernel_ = kernel(**params)
            self.initialised_["kernel"] = True

    def make_kernel_(self, kernel_list, hidden_args):
        """Produce the desired kernel function.

        Parameters
        ----------
        kernel_list: (list)
            List of kernel dictionaries as defined at the documentation
            of class parameters.

        Returns
        -------
        kernel : kernel (class).
            Returns an instance of a kernel type object corresponding to the
            certain kernel.

        """
        kernel = kernel_list.pop(0)
        if type(kernel) is not dict:
            raise ValueError('each element of the list of kernels must'
                             ' be a dictionary')
        if "name" not in kernel:
            raise ValueError('each dictionary concerning a kernel must'
                             ' have a "name" parameter designating the'
                             'kernel')
        kernel_name = kernel.pop("name")
        for (keys, val) in iteritems(hidden_args):
            kernel[keys] = val
        if kernel_name in supported_base_kernels:
            if len(kernel_list) != 0:
                warnings.warn('rest kernel arguments are being ignored\
                               - reached base kernel')
            if kernel_name == "subtree_wl":
                return subtree_wl, kernel
            elif kernel_name == "random_walk":
                return random_walk, kernel
            elif kernel_name == "shortest_path":
                if kernel.pop("as_attributes", False):
                    return shortest_path_attr, kernel
                else:
                    return (shortest_path, kernel)
            elif kernel_name == "subtree_rg":
                raise ValueError('still developing')
            elif kernel_name == "graphlet_sampling":
                if ("random_seed" not in kernel and
                    self.random_seed is not
                        default_random_seed_value):
                        kernel["random_seed"] = self.random_seed
                return graphlet_sampling, kernel
            elif kernel_name == "multiscale_laplacian":
                if kernel.pop("which", None) == "simple":
                    kernel.pop("N", None)
                    return (multiscale_laplacian, kernel)
                else:
                    if ("random_seed" not in kernel and
                        self.random_seed is not
                            default_random_seed_value):
                        kernel["random_seed"] = self.random_seed
                    return (multiscale_laplacian_fast, kernel)
            elif kernel_name == "subgraph_matching":
                return subgraph_matching, kernel
            elif kernel_name == "lovasz_theta":
                if ("random_seed" not in kernel and
                        self.random_seed is not
                        default_random_seed_value):
                    kernel["random_seed"] = self.random_seed
                return lovasz_theta, kernel
            elif kernel_name == "svm_theta":
                if ("random_seed" not in kernel and
                    self.random_seed is not
                        default_random_seed_value):
                    kernel["random_seed"] = self.random_seed
                return svm_theta, kernel
            elif kernel_name == "neighborhood_hash":
                return neighborhood_hash, kernel
            elif kernel_name == "neighborhood_subgraph_pairwise_distance":
                return neighborhood_subgraph_pairwise_distance, kernel
            elif kernel_name == "odd_sth":
                return odd_sth, kernel
            elif kernel_name == "propagation":
                if ("random_seed" not in kernel and
                    self.random_seed is not
                        default_random_seed_value):
                    kernel["random_seed"] = self.random_seed
                return propagation, kernel
            elif kernel_name == "pyramid_match":
                return pyramid_match, kernel
            elif kernel_name == "jsm":
                return jsm, kernel
            elif kernel_name == "vertex_histogram":
                return vertex_histogram, kernel
            elif kernel_name == "edge_histogram":
                return edge_histogram, kernel
        elif kernel_name in supported_general_kernels:
            if (len(kernel_list) == 0):
                raise ValueError(str(kernel_name)+' is not a base kernel')
            else:
                kernel["base_kernel"] = self.make_kernel_(kernel_list, {})
            if kernel_name == "weisfeiler_lehman":
                return (weisfeiler_lehman, kernel)
            if kernel_name == "hadamard_code":
                return (hadamard_code, kernel)
        else:
            raise ValueError("unsupported kernel: " + str(kernel_name))

    def set_params(self, **params):
        """Call the parent method."""
        # Copy the parameters
        params = copy.deepcopy(params)

        # Iterate over the parameters
        for key, value in iteritems(params):
            key, delim, sub_key = key.partition('__')
            if delim:
                if sub_key in self.initialised_:
                    self.initialised_[sub_key] = False
            elif key in self.initialised_:
                self.initialised_[key] = False

        # Set parameters
        super(GraphKernel, self).set_params(**params)
