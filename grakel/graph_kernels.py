"""The main graph kernel class, implemented as a sci-kit transformer."""
import collections
import itertools
import time
import warnings

import numpy as np

from concurrent.futures import ThreadPoolExecutor
from scipy.linalg import svd
from sklearn.base import BaseEstimator
from sklearn.base import TransformerMixin
from sklearn.utils.validation import check_is_fitted

from grakel.graph import graph
from grakel.kernels.dirac import dirac_pair
from grakel.kernels.graphlet_sampling import graphlet_sampling_matrix
from grakel.kernels.hadamard_code import hadamard_code_matrix
from grakel.kernels.jsm import jsm_pair
from grakel.kernels.lovasz_theta import lovasz_theta_pair
from grakel.kernels.multiscale_laplacian import multiscale_laplacian_pair
from grakel.kernels.neighborhood_hash import neighborhood_hash_matrix
from grakel.kernels.neighborhood_subgraph_pairwise_distance\
    import neighborhood_subgraph_pairwise_distance_pair
from grakel.kernels.odd_sth import odd_sth_matrix
from grakel.kernels.propagation import propagation_matrix
from grakel.kernels.pyramid_match import pyramid_match_matrix
from grakel.kernels.random_walk import random_walk_pair
from grakel.kernels.shortest_path import shortest_path_matrix
from grakel.kernels.shortest_path import shortest_path_pair_attributes
from grakel.kernels.subgraph_matching import subgraph_matching_pair
from grakel.kernels.subtree_rg import subtree_rg_pair
from grakel.kernels.svm_theta import svm_theta_pair
from grakel.kernels.weisfeiler_lehman import weisfeiler_lehman_matrix

np.random.seed(int(time.time()))

global supported_base_kernels, supported_general_kernels, default_n_components

valid_parameters = {"kernel",
                    "Nystroem",
                    "concurrency",
                    "normalize",
                    "verbose"}

supported_base_kernels = [
    "dirac", "random_walk",
    "shortest_path", "subtree_rg",
    "graphlets_sampling", "subgraph_matching",
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

default_n_components = 100


class GraphKernel(BaseEstimator, TransformerMixin):
    """A general class that describes all kernels.

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
                            + (**o**) method_type : [str] "sylvester", "simple"

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
    Attributes
    ----------
    X_graph : dict
        Stores input graphs indexing from 0 to num_of_graphs-1.

    num_of_graphs : int
        The number of graphs given in fit input.

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

    """

    num_of_graphs = 0
    X_graph = None
    nystroem = False
    normalize = False

    def __init__(self, **kargs):
        """`__init__` for `GraphKernel` object."""
        self._verbose = kargs.get("verbose", default_verbose_value)

        if "kernel" in kargs:
            if (type(kargs["kernel"]) is dict):
                # allow single kernel dictionary inputs
                kernel, flag = self._make_kernel([kargs["kernel"]])
            elif (type(kargs["kernel"]) is list):
                kernel, flag = self._make_kernel(kargs["kernel"])
            else:
                raise ValueError('unsupported kernel format')
            if flag:
                self.kernel = lambda y=None, diagonal_only=False: \
                    self.calculate_kernel_matrix(
                        kernel,
                        target_graph=y,
                        kernel_type="matrix",
                        diagonal_only=diagonal_only)
            else:
                self.kernel = lambda y=None, diagonal_only=False:\
                    self.calculate_kernel_matrix(
                        kernel,
                        target_graph=y,
                        kernel_type="pairwise",
                        diagonal_only=diagonal_only)
        else:
            raise ValueError('kernel must be defined at the __init__ function \
                              of a graph kernel')

        if "Nystroem" in kargs:
            if type(kargs["Nystroem"]) is not int:
                raise ValueError('nystroem parameter must be an int, \
                                  indicating the number of components')
            elif kargs["Nystroem"] in [0, -1]:
                # picking default number of components
                self.nystroem = default_n_components
            elif kargs["Nystroem"] <= 0:
                raise ValueError('number of nystroem components \
                                  must be positive')
            else:
                self.nystroem = kargs["Nystroem"]

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

        self.normalize = kargs.get("normalize", False)

        unrecognised_args = set(kargs.keys()) - valid_parameters

        if len(unrecognised_args) > 0:
            warnings.warn(
                'Ignoring unrecognised arguments:',
                ', '.join('"' + str(arg) + '"' for arg in unrecognised_args))

    def _make_kernel(self, kernel_list):
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
            kernel_name = kernel.pop("name")
            if kernel_name in supported_base_kernels:
                if (len(kernel_list) != 0):
                    warnings.warn('rest kernel arguments are being ignored\
                                   - reached base kernel')
                if kernel_name == "dirac":
                    return (lambda x, y: dirac_pair(x, y), False)
                elif kernel_name == "random_walk":
                    return (lambda x, y: random_walk_pair(
                        x, y, **kernel), False)
                elif kernel_name == "shortest_path":
                    at = kernel.get("algorithm_type", "dijkstra")
                    if kernel.get("as_attributes", False):
                        ak = kernel.get("attribute_kernel",
                                        lambda x, y: np.dot(x, y))
                        return (lambda x, y:
                                shortest_path_pair_attributes(
                                    x,
                                    y,
                                    algorithm_type=at,
                                    attribute_kernel=ak), False)
                    else:
                        wl = kernel.get("with_labels", True)
                        return (lambda x, y:
                                shortest_path_matrix(x, y,
                                                     algorithm_type=at,
                                                     with_labels=wl),
                                True)
                elif kernel_name == "subtree_rg":
                    return (lambda x, y:
                            subtree_rg_pair(x, y, **kernel), False)
                elif kernel_name == "graphlets_sampling":
                    return (lambda x, y:
                            graphlet_sampling_matrix(x, y, **kernel), True)
                elif kernel_name == "multiscale_laplacian":
                    return (lambda x, y:
                            multiscale_laplacian_pair(x, y, **kernel), False)
                elif kernel_name == "subgraph_matching":
                    return (lambda x, y:
                            subgraph_matching_pair(x, y, **kernel), False)
                elif kernel_name == "lovasz_theta":
                    return (lambda x, y:
                            lovasz_theta_pair(x, y, **kernel), False)
                elif kernel_name == "svm_theta":
                    return (lambda x, y:
                            svm_theta_pair(x, y, **kernel), False)
                elif kernel_name == "neighborhood_hash":
                    return (lambda x, y:
                            neighborhood_hash_matrix(x, y, **kernel), True)
                elif kernel_name == "neighborhood_subgraph_pairwise_distance":
                    return (lambda x, y:
                            neighborhood_subgraph_pairwise_distance_pair(
                                x, y, **kernel), False)
                elif kernel_name == "odd_sth":
                    return (lambda x, y:
                            odd_sth_matrix(x, y, **kernel), True)
                elif kernel_name == "propagation":
                    return (lambda x, y:
                            propagation_matrix(x, y, **kernel), True)
                elif kernel_name == "pyramid_match":
                    return (lambda x, y:
                            pyramid_match_matrix(x, y, **kernel), True)
                elif kernel_name == "jsm":
                    return (lambda x, y:
                            jsm_pair(x, y, **kernel), False)
            elif kernel_name in supported_general_kernels:
                if (len(kernel_list) == 0):
                    raise ValueError(str(kernel_name)+' is not a base kernel')

                (bkp, matrix_flag) = self._make_kernel(kernel_list)
                if not matrix_flag:
                    bk = lambda x, y=None: \
                        pairwise_to_matrix_kernel(
                            x=x,
                            pairwise_kernel=bkp,
                            kernel_executor=self._pairwise_kernel_executor,
                            y=y)
                else:
                    bk = bkp

                if kernel_name == "weisfeiler_lehman":
                    return (lambda x, y:
                            weisfeiler_lehman_matrix(x, bk, y, **kernel), True)
                if kernel_name == "hadamard_code":
                    return (lambda x, y:
                            hadamard_code_matrix(x, bk, y, **kernel), True)
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
        if X is None:
            raise ValueError('fit input cannot be None')
        else:
            X_graph = parse_input(X)
            num_of_graphs = len(X_graph)

        if bool(self.nystroem):
            # get basis vectors
            if self.nystroem > num_of_graphs:
                n_components = num_of_graphs
                warnings.warn("n_components > n_samples. This is not \
                              possible.\nn_components was set to n_samples, \
                              which results in inefficient evaluation of the \
                              full kernel.")
            else:
                n_components = self.nystroem

            n_components = min(num_of_graphs, n_components)
            inds = np.random.permutation(num_of_graphs)
            basis_inds = inds[:n_components]
            basis = {i: X_graph[idx] for (i, idx) in enumerate(basis_inds)}

            self.components_ = basis
            self.nystroem = n_components
            basis_kernel = self.kernel()

            # sqrt of kernel matrix on basis vectors
            U, S, V = svd(basis_kernel)
            S = np.maximum(S, 1e-12)
            self.nystroem = n_components
            self.nystroem_normalization_ = np.dot(U / np.sqrt(S), V)
            self.components_ = basis
            self.component_indices_ = inds

        self.X_graph = X_graph
        self.num_of_graphs = num_of_graphs

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
        # Check is fit had been called
        check_is_fitted(self, ['X_graph'])
        if bool(self.nystroem):
            check_is_fitted(self, 'components_')

        # Input validation and parsing
        if X is None:
            raise ValueError('transform input cannot be None')
        else:
            target_graph = parse_input(X)

        # Transform - calculate kernel matrix
        if bool(self.nystroem):
            return np.dot(
                self.kernel(target_graph), self.nystroem_normalization_.T)
        else:
            km = self.kernel(target_graph)
            if self.normalize:
                dX, dT = self.kernel(y=target_graph, diagonal_only=True)
                km /= np.sqrt(np.dot(dT, dX.T))
            return km

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
        self.fit(X)

        # Transform - calculate kernel matrix
        if bool(self.nystroem):
            return np.dot(self.kernel(self.X_graph),
                          self.nystroem_normalization_.T)
        else:
            km = self.kernel()
            if self.normalize:
                km_diag = np.diagonal(km)
                km_diag = km_diag.reshape(km_diag.shape[0], 1)
                km /= np.sqrt(np.dot(km_diag.T, km_diag))
            return km

    def calculate_kernel_matrix(
            self,
            kernel,
            target_graph=None,
            diagonal_only=False,
            kernel_type="pairwise"):
        """Calculate the kernel matrix given a target_graph and a kernel.

        Parameters
        ----------
        target_graph : dict
            A dictionary from 0 to the number of graphs of
            target "graph objects".

        kernel: function
            A pairwise graph kernel (between "graph" type objects)
            or a matrix graph kernel (between enumerative dictionaries -
            starting from 0 - for graph type objects)

        kernel_type: str, valid_values={"matrix", "pairwise"}
            Distinguishes between the two valid types for graphs.

        diagonal_only : bool, default=False
            Calculate only the diagonal.
        Returns
        -------
        K : numpy array, shape = [n_targets, n_input_graphs], case_of_existence
        =(diagonal_only==False)
            The kernel matrix: a calculation between all pairs of graphs
            between target an features.

        diag_X : numpy array, shape = [n_input_graphs], case_of_existence=
        (diagonal_only==True)
            Returns the kernel values of each element of the input_graph
            with itself. Useful for normalization.

        diag_T : numpy array, shape = [n_targets], case_of_existence=
        (diagonal_only==True)
            Returns the kernel values of each element of the input_graph
            with itself. Useful for normalization.

        """
        if kernel_type not in ["matrix", "pairwise"]:
            raise ValueError('unsupported "kernel_type"')

        if bool(self.nystroem):
            X_graph = self.components_
            num_of_graphs = self.nystroem
        else:
            X_graph = self.X_graph
            num_of_graphs = self.num_of_graphs
        if kernel_type == "pairwise":
            if diagonal_only:
                diag_X = np.empty(shape=(num_of_graphs, 1))
                for i in range(num_of_graphs):
                    diag_X[i] = self._pairwise_kernel_executor(kernel,
                                                               X_graph[i],
                                                               X_graph[i])
                if target_graph is None:
                    return diag_X
                else:
                    diag_T = np.empty(shape=(len(target_graph), 1))
                    for i in range(len(target_graph)):
                        diag_T[i] = self._pairwise_kernel_executor(
                            kernel,
                            target_graph[i],
                            target_graph[i])
                    return diag_X, diag_T
            else:
                if target_graph is None:
                    is_symmetric = True
                    target_graph = X_graph
                    num_of_targets = num_of_graphs
                    pairs = [(i, j) for i in range(0, num_of_targets)
                             for j in range(i, num_of_graphs)]
                else:
                    is_symmetric = False
                    num_of_targets = len(target_graph)
                    pairs = list(itertools.product(range(0, num_of_targets),
                                                   range(0, num_of_graphs)))

                K = np.zeros(shape=(num_of_targets, num_of_graphs))

                for (i, j) in pairs:
                    K[i, j] = self._pairwise_kernel_executor(kernel,
                                                             target_graph[i],
                                                             X_graph[j])

                if is_symmetric:
                    K = np.triu(K) + np.triu(K, 1).T

            return K
        else:
            if diagonal_only:
                diag_X = np.empty(shape=(num_of_graphs, 1))
                for i in range(num_of_graphs):
                    diag_X[i] = self._pairwise_kernel_executor(
                        kernel,
                        {0: X_graph[i]},
                        {0: X_graph[i]})[0, 0]
                if target_graph is None:
                    return diag_X
                else:
                    diag_T = np.empty(shape=(len(target_graph), 1))
                    for i in range(len(target_graph)):
                        diag_T[i] = self._pairwise_kernel_executor(
                            kernel,
                            {0: target_graph[i]},
                            {0: target_graph[i]})[0, 0]
                    return diag_X, diag_T
            else:
                return kernel(X_graph, target_graph)


def pairwise_to_matrix_kernel(x, pairwise_kernel, kernel_executor, y=None):
    """Convert a pairwise-kernel to a matrix-kernel.

    Parameters
    ----------
    x,y : dict, default_y=None
        Enumerative dictionary of graph type objects with keys from 0 to the
        number of values. If value of Graphs_y is None the kernel matrix is
        computed between all pairs of Graphs_x, where in another case the
        kernel_matrix rows correspond to elements of Graphs_y, and columns
        to the elements of Graphs_x.

    pairwise_kernel : function
        A pairwise graph kernel (between "graph" type objects)

    kernel_excutor : Executor
        An executor for functions to allow concurrency speedups.

    Returns
    -------
    kernel_matrix : np.array
        The kernel matrix. If the Graphs_y is not None the rows correspond to
        Graphs_y and the cols to the Graphs_x (based on the given order).

    """
    nx = len(x.keys())

    if y is None:
        ny = nx
        ng = nx
        pairs = [(i, j) for i in range(0, nx) for j in range(i, ny)]
        offset = 0
        Gs = x
    else:
        ny = len(y.keys())
        ng = nx + ny
        pairs = list(itertools.product(range(nx, ng), range(0, nx)))
        offset = nx
        Gs = {i: g for (i, g) in enumerate(itertools.chain(x.values(),
                                           y.values()))}

    kernel_mat = np.zeros(shape=(ny, nx))

    for (i, j) in pairs:
        kernel_mat[i-offset, j] = \
            kernel_executor(pairwise_kernel, Gs[i], Gs[j])

    if y is None:
        kernel_mat = np.triu(kernel_mat) + np.triu(kernel_mat, 1).T

    return kernel_mat


def parse_input(X):
    """Parse the given input and raise errors if it is invalid.

    Parameters
    ----------
    X : object
        For the input to pass the test, we must have:
        Each element must be an iterable with at most three features and at
        least one. The first that is obligatory is a valid graph structure
        (adjacency matrix or edge_dictionary) while the second is node_labels
        and the third edge_labels (that fitting the given graph format).
        If None the kernel matrix is calculated upon fit data.
        The test samples.

    Returns
    -------
    X_graph : dict
        Enumerative dictionary of graph type objects with keys from 0 to the
        number of values.

    """
    if not isinstance(X, collections.Iterable):
        raise ValueError('input must be an iterable\n')
        # Not a dictionary
    else:
        X_graph = dict()
        for (i, x) in enumerate(iter(X)):
            if len(x) == 0:
                warnings.warn('Ignoring empty element on index: '+str(i))
            if len(x) == 1:
                X_graph[len(X_graph)] = graph(x[0], {}, {}, "all")
            elif len(x) == 2:
                X_graph[len(X_graph)] = graph(x[0], x[1], {}, "all")
            elif len(x) == 3:
                X_graph[len(X_graph)] = graph(x[0], x[1], x[2], "all")
            else:
                raise ValueError('each element of X must have at least one \
                                  and at most 3 elements\n')
        if len(X_graph) == 0:
            raise ValueError('parsed input is empty')
        return X_graph
