"""The main graph kernel class, implemented as a sci-kit transformer."""
import copy
import warnings

import numpy as np

from scipy.linalg import svd
from sklearn.base import BaseEstimator
from sklearn.base import TransformerMixin
from sklearn.utils import check_random_state
from sklearn.utils.validation import check_is_fitted

from grakel.kernels import GraphletSampling
from grakel.kernels import RandomWalk
from grakel.kernels import RandomWalkLabeled
from grakel.kernels import ShortestPath
from grakel.kernels import ShortestPathAttr
from grakel.kernels import WeisfeilerLehman
from grakel.kernels import NeighborhoodHash
from grakel.kernels import PyramidMatch
from grakel.kernels import SubgraphMatching
from grakel.kernels import NeighborhoodSubgraphPairwiseDistance
from grakel.kernels import LovaszTheta
from grakel.kernels import SvmTheta
from grakel.kernels import OddSth
from grakel.kernels import Propagation
from grakel.kernels import PropagationAttr
from grakel.kernels import HadamardCode
from grakel.kernels import MultiscaleLaplacian
from grakel.kernels import VertexHistogram
from grakel.kernels import EdgeHistogram
from grakel.kernels import GraphHopper
from grakel.kernels import CoreFramework
from grakel.kernels import WeisfeilerLehmanOptimalAssignment

# Python 2/3 cross-compatibility import
from future.utils import iteritems

# Supported base kernels
sbk = [
    ["vertex_histogram", "subtree_wl", "VH", "ST-WL"],
    ["edge_histogram", "EH"],
    ["random_walk", "RW"],
    ["shortest_path", "SP"],
    ["graphlet_sampling", "GR"],
    ["subgraph_matching", "SM"],
    ["multiscale_laplacian", "ML"],
    ["lovasz_theta", "LOVT"],
    ["svm_theta", "SVMT"],
    ["neighborhood_hash", "NH"],
    ["neighborhood_subgraph_pairwise_distance", "NSPD"],
    ["odd_sth", "ODD"],
    ["propagation", "PR"],
    ["pyramid_match", "PM"],
    ["graph_hopper", "GH"],
    ["weisfeiler_lehman_optimal_assignment", "WL-OA"]
    ]

sbks = set(e for ls in sbk for e in ls)

# Supported frameworks
sf = [
    ["weisfeiler_lehman", "WL"],
    ["hadamard_code", "HC"],
    ["core_framework", "CORE"]
    ]

sfs = set(e for ls in sf for e in ls)

# Supported kernels message
sep = u"\n \u27E1 "
sk_msg = ("Base-Kernels\n" + 12*"-" + sep + sep.join(','.join(synonyms) for synonyms in sbk) +
          "\n\nFrameworks\n" + 10*"-" + sep + sep.join(','.join(synonyms) for synonyms in sf))

# Defaults
default_n_components = 100


class GraphKernel(BaseEstimator, TransformerMixin):
    r"""A generic wrapper for graph kernels.

    Parameters
    ----------
    kernel : list*(dict or str)
        A single element or a list of :code:`dict` with:

            * "name" : [str] - with the kernel name

            * "name_of_parameter_1" : value

            * "name_of_parameter_2" : value

            * :math:`\;\cdots\;`

            * "name_of_parameter_k" : value

        or of :code:`str`, designating a kernel name.

        available "name" or "name-alias" / "parametres" are:
            1. base_graph_kernels (the structure must always reach a base kernel)

                - "vertex_histogram" or "subtree_wl" or "VH" or "ST-WL"
                    + (**o**) "sparse" : bool or 'auto'

                - "edge_histogram" or "EH"
                    + (**o**) "sparse" : bool or 'auto'

                - "random_walk" or "RW"
                    + (**o**) "with_labels" : bool

                    + (**o**) "lamda" : float

                    + (**o**) "method_type" : [str], "baseline", "fast"

                    + (**o**) "kernel_type" : [str], "geometric", "exponential"

                    + (**o**) "p" : [int] > 0


                - "shortest_path" or "SP"
                    + (**o**) "algorithm_type" : [str] "dijkstra", "floyd_warshall"

                    + (**o**) "as_attributes" : [bool]

                    + (**o**) "metric" : [function] : (attribute_x, attribute_y) -> number

                    + (**o**) "with_labels" : [bool]

                - "graphlet_sampling" or "GR"
                    + (**o**) "k" : [int]

                    + (**o**) "sampling" : [dict] or **None**

                - "multiscale_laplacian" or "ML"
                    + (**o**) "L" : [int] > 0

                    + (**o**) "gamma" : [float] > .0

                    + (**o**) "heta" : [float] > .0

                    + (**o**) "n_samples" : [int] > 0, if "which": "fast"

                    + (**o**) "P" : [int] > 0, if "which": "fast"

                - "subgraph_matching" or "SM"
                    + (**o**) "kv" : [function] : (node_x, node_y, Lx, Ly) -> number

                    + (**o**) "ke" : [function] : (edge_x, edge_y, Lx, Ly) -> number

                    + (**o**) "lw" : a lambda weight function for cliques: set -> number

                - "lovasz_theta" or "LOVT"
                    + (**o**) "n_samples" : [int] > 1

                    + (**o**) "subsets_size_range" : [tuple] of two [int]

                    + (**o**) "metric" : [function] (number, number) -> number

                - "svm_theta" or "SVMT"
                    + (**o**) "n_samples" : [int] > 1

                    + (**o**) "subsets_size_range" : [tuple] with 2 [int] elements

                    + (**o**) "metric" : [function] (number, number) -> number

                - "neighborhood_hash" or "NH"
                    + (**o**) "nh_type" : [str] "simple" or "count-sensitive"

                    + (**o**) "R" : [int] > 0

                    + (**o**) "bits" : [int] > 0

                - "neighborhood_subgraph_pairwise_distance" or "NSPD"
                    + (**o**) "r" : (int) positive integer

                    + (**o**) "d" : (int) positive integer

                - "odd_sth" or "ODD"
                    + (**o**) "h" : [int] > 0

                - "propagation" or "PR"
                    + (**o**) t_max: [int] > 0

                    + (**o**) T: [dict] [int]: [np.arrays]

                    + (**o**) with_attributes: [bool], default=False

                    + (**o**) M: [str] {"H", "TV"} if `with_attributes=True` else {"L1", "L2"}

                    + (**o**) w: [int] > 0

                    + (**o**) metric: [function] x:[Counter] , y:[Counter] -> [number]

                - "pyramid_match" or "PM"
                    + (**o**) with_labels: [bool]

                    + (**o**) d: [int] > 0

                    + (**o**) L: [int] >= 0

                - "graph_hopper" or "GH"
                    + (**o**) kernel_type: [str: {'linear', 'gaussian'}] or [tuple: {('gaussian', mu)}]
                      or [function] x:[(np.array, np.array)] , y:[(np.array, np.array)] -> [number]

                - "weisfeiler_lehman_optimal_assignment" or "WL-OA"
                    + (**o**) "n_iter" : [int] >= 0

            2. frameworks (if a next kernel in the list it asssigned as a base-kernel, else see default)
                - "weisfeiler_lehman" or "WL" / default="VH"
                    + (**o**) "n_iter" : [int] >= 0

                - "hadamard_code" or "HC" / default="VH"
                    + (**o**) "n_iter" : [int] > 0

                - "core_framework" or "CORE" / default="SP"
                    + (**o**) "min_core" : [int] >= -1

        where (**o**): stands for optional parameters

    Nystroem : int or bool, optional
        Defines the number of nystroem components.
        To initialize the default (100 components), set -1 or 0.

    n_jobs : int or None, optional
        Defines the number of jobs of a joblib.Parallel objects needed for parallelization
        or None for direct execution. The use or not of this function depends on each kernel.

    normalize : bool, optional
        Normalize the output of the graph kernel.
        Ignored when Nystroem GraphKernel object is instanciated.

    verbose : bool, optional
        Define if messages will be printed on stdout.

    random_state :  RandomState or int, default=None
        A random number generator instance or an int to initialize a RandomState as a seed.

    Attributes
    ----------
    _initialized : dict
        Monitors which parameter derived object should be _initialized.

    kernel_ : function
        The full kernel applied between graph objects.

    nystroem_ : int
        Holds the nystroem, number of components.
        If not _initialized, it stands as a False
        boolean variable.

    components_ : array, shape=(n_components, n_features)
        Subset of training graphs used to construct the feature map.

    nystroem_normalization_ : array, shape=(n_components, n_components)
        Normalization matrix needed for embedding.
        Square root of the kernel matrix on ``components_``.

    component_indices_ : array, shape=(n_components)
        Indices of ``components_`` in the training set.

    random_state_ : RandomState
        A RandomState object handling all randomness of the class.

    """

    def __init__(self,
                 kernel="shortest_path",
                 normalize=False,
                 verbose=False,
                 n_jobs=None,
                 random_state=None,
                 Nystroem=False):
        """`__init__` for `GraphKernel` object."""
        self.kernel = kernel
        self.normalize = normalize
        self.verbose = verbose
        self.n_jobs = n_jobs
        self.random_state = random_state
        self.Nystroem = Nystroem
        self._initialized = {"kernel": False,
                             "Nystroem": False,
                             "random_state": False,
                             "normalize": False,
                             "verbose": False,
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
        # Initialize the Graph Kernel.
        self.initialize()

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
            inds = self.random_state_.permutation(nx)
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
            self.kernel_.fit(X)

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
        # Transform - calculate kernel matrix
        check_is_fitted(self, 'kernel_')
        if hasattr(self, 'nystroem_') and bool(self.nystroem_):
            # Check if nystroem has been initialized had been called
            check_is_fitted(self, 'components_')
            K = self.kernel_.transform(X).dot(self.nystroem_normalization_.T)
        else:
            K = self.kernel_.transform(X)

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
        # Initialize the Graph Kernel
        self.initialize()

        # Transform - calculate kernel matrix
        if bool(self.nystroem_):
            self.fit(X)
            K = self.kernel_.transform(X).dot(self.nystroem_normalization_.T)
        else:
            K = self.kernel_.fit_transform(X)

        return K

    def initialize(self):
        """Initialize all transformer arguments, needing initialisation."""
        if not self._initialized["Nystroem"]:
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
            self._initialized["Nystroem"] = True

        if any(not self._initialized[param]
               for param in ["random_state", "normalize", "verbose", "n_jobs", "kernel"]):
            # Intialise random_state_
            if not self._initialized["random_state"]:
                self.random_state_ = check_random_state(self.random_state)

            k = self.kernel
            if type(k) is dict or type(k) is str:
                # allow single kernel dictionary inputs
                k = [self.kernel]
            elif type(k) is not list:
                raise ValueError('A "kernel" must be defined at the __init__ '
                                 'function of the graph kernel generic wrapper.'
                                 'Valid kernel types are dict, str, and list of dict or str.')

            hidden_args = {"verbose": self.verbose, "normalize": self.normalize,
                           "n_jobs": self.n_jobs}

            # Initialize a new kernel each time a new fit is being called
            kernel, params = self.make_kernel_(copy.deepcopy(k), hidden_args)
            self.kernel_ = kernel(**params)
            for param in ["random_state", "normalize", "verbose", "n_jobs", "kernel"]:
                self._initialized[param] = True

    def make_kernel_(self, kernel_list, hidden_args):
        """Produce the desired kernel function.

        Parameters
        ----------
        kernel_list : (list)
            List of kernel dictionaries as defined at the documentation
            of class parameters.

        Returns
        -------
        kernel : kernel (class).
            Returns an instance of a kernel type object corresponding to the
            certain kernel.

        """
        kernel = kernel_list.pop(0)
        if type(kernel) is str:
            kernel_name, kernel = str(kernel), dict()
        elif type(kernel) is not dict:
            raise ValueError('each element of the list of kernels must be a dictionary or a string')
        else:
            if "name" not in kernel:
                raise ValueError('each dictionary concerning a kernel must '
                                 'have a "name" parameter designating the '
                                 'kernel')
            kernel_name = kernel.pop("name")

        for (keys, val) in iteritems(hidden_args):
            if keys in kernel:
                warnings.warn('Overriding global kernel attribute ' + str(keys) + ' with ' + str(val) +
                              '. Please set this attribute as an argument of GraphKernel.')
            kernel[keys] = val

        def get_random_state_(kernel):
            return kernel.pop(
                "random_state",
                (self.random_state_ if self.random_state is not None else None))

        if kernel_name in sbks:
            if len(kernel_list) != 0:
                warnings.warn('Kernel List not empty while reaching a base-kernel - '
                              'the rest kernel names will be ignored')

            if kernel_name in sbk[0]:
                return VertexHistogram, kernel
            elif kernel_name in sbk[1]:
                return EdgeHistogram, kernel
            elif kernel_name in sbk[2]:
                if kernel.pop("with_labels", False):
                    return RandomWalkLabeled, kernel
                else:
                    return RandomWalk, kernel
            elif kernel_name in sbk[3]:
                if kernel.pop("as_attributes", False):
                    return ShortestPathAttr, kernel
                else:
                    return ShortestPath, kernel
            elif kernel_name in sbk[4]:
                kernel["random_state"] = get_random_state_(kernel)
                return GraphletSampling, kernel
            elif kernel_name in sbk[5]:
                return SubgraphMatching, kernel
            elif kernel_name in sbk[6]:
                kernel["random_state"] = get_random_state_(kernel)
                return (MultiscaleLaplacian, kernel)
            elif kernel_name in sbk[7]:
                kernel["random_state"] = get_random_state_(kernel)
                return LovaszTheta, kernel
            elif kernel_name in sbk[8]:
                kernel["random_state"] = get_random_state_(kernel)
                return SvmTheta, kernel
            elif kernel_name in sbk[9]:
                return NeighborhoodHash, kernel
            elif kernel_name in sbk[10]:
                return NeighborhoodSubgraphPairwiseDistance, kernel
            elif kernel_name in sbk[11]:
                return OddSth, kernel
            elif kernel_name in sbk[12]:
                kernel["random_state"] = get_random_state_(kernel)
                if kernel.pop("with_attributes", False):
                    return PropagationAttr, kernel
                else:
                    return Propagation, kernel
            elif kernel_name in sbk[13]:
                return PyramidMatch, kernel
            elif kernel_name in sbk[14]:
                return GraphHopper, kernel
            elif kernel_name in sbk[15]:
                return WeisfeilerLehmanOptimalAssignment, kernel

        elif kernel_name in sfs:
            if len(kernel_list):
                kernel["base_graph_kernel"] = self.make_kernel_(kernel_list, {})
            if kernel_name in sf[0]:
                return (WeisfeilerLehman, kernel)
            elif kernel_name in sf[1]:
                return (HadamardCode, kernel)
            elif kernel_name in sf[2]:
                return (CoreFramework, kernel)
        else:
            raise ValueError("Unsupported kernel: " + str(kernel_name) + "\n"
                             "Supported kernels are:\n\n" + sk_msg)

    def set_params(self, **params):
        """Call the parent method."""
        # Copy the parameters
        params = copy.deepcopy(params)

        # Iterate over the parameters
        for key, value in iteritems(params):
            key, delim, sub_key = key.partition('__')
            if delim:
                if sub_key in self._initialized:
                    self._initialized[sub_key] = False
            elif key in self._initialized:
                self._initialized[key] = False

        # Set parameters
        super(GraphKernel, self).set_params(**params)
