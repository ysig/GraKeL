"""The core kernel framework as defined in :cite:`nikolentzos2018degeneracy`."""
# Author: Ioannis Siglidis <y.siglidis@gmail.com>
# License: BSD 3 clause
import collections
import warnings

import numpy as np

from sklearn.exceptions import NotFittedError
from sklearn.utils.validation import check_is_fitted

from grakel.graph import Graph
from grakel.kernels import Kernel
from grakel.kernels.shortest_path import ShortestPath

# Python 2/3 cross-compatibility import
from six import iteritems


class CoreFramework(Kernel):
    """The core kernel framework, as proposed in :cite:`nikolentzos2018degeneracy`.

    Parameters
    ----------
    base_graph_kernel : `grakel.kernels.kernel` or tuple, default=None
        If tuple it must consist of a valid kernel object and a
        dictionary of parameters. General parameters concerning
        normalization, concurrency, .. will be ignored, and the
        ones of given on `__init__` will be passed in case it is needed.
        Default `base_graph_kernel` is `VertexHistogram`.

    min_core : int, default=-1
        Core numbers bigger than min_core will only be considered.

    Attributes
    ----------
    base_graph_kernel_ : function
        A void function that initializes a base kernel object.

    """

    _graph_format = "dictionary"

    def __init__(self, n_jobs=None, verbose=False,
                 normalize=False, min_core=-1, base_graph_kernel=None):
        """Initialise a `hadamard_code` kernel."""
        super(CoreFramework, self).__init__(
            n_jobs=n_jobs, verbose=verbose, normalize=normalize)

        self.min_core = -1
        self.base_graph_kernel = base_graph_kernel
        self._initialized.update({"min_core": False, "base_graph_kernel": False})

    def initialize(self):
        """Initialize all transformer arguments, needing initialization."""
        if not self._initialized["n_jobs"]:
            if self.n_jobs is not None:
                warnings.warn('no implemented parallelization for CoreFramework')
            self._initialized["n_jobs"] = True

        if not self._initialized["base_graph_kernel"]:
            base_graph_kernel = self.base_graph_kernel
            if base_graph_kernel is None:
                base_graph_kernel, params = ShortestPath, dict()
            elif type(base_graph_kernel) is type and issubclass(base_graph_kernel, Kernel):
                params = dict()
            else:
                try:
                    base_graph_kernel, params = base_graph_kernel
                except Exception:
                    raise TypeError('Base kernel was not formulated in '
                                    'the correct way. '
                                    'Check documentation.')

                if not (type(base_graph_kernel) is type and
                        issubclass(base_graph_kernel, Kernel)):
                    raise TypeError('The first argument must be a valid '
                                    'grakel.kernel.kernel Object')
                if type(params) is not dict:
                    raise ValueError('If the second argument of base '
                                     'kernel exists, it must be a diction'
                                     'ary between parameters names and '
                                     'values')
                params.pop("normalize", None)

            params["normalize"] = False
            params["verbose"] = self.verbose
            params["n_jobs"] = None
            self.base_graph_kernel_ = base_graph_kernel
            self.params_ = params
            self._initialized["base_graph_kernel"] = True

        if not self._initialized["min_core"]:
            if type(self.min_core) is not int or self.min_core < -1:
                raise TypeError("'min_core' must be an integer bigger than -1")
            self._initialized["min_core"] = True

    def parse_input(self, X):
        """Parse input and create features, while initializing and/or calculating sub-kernels.

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
        base_graph_kernel : object
            Returns base_graph_kernel. Only if called from `fit` or `fit_transform`.

        K : np.array
            Returns the kernel matrix. Only if called from `transform` or
            `fit_transform`.

        """
        # Input validation and parsing
        if not isinstance(X, collections.Iterable):
            raise TypeError('input must be an iterable\n')
        else:
            nx, max_core_number, core_numbers, graphs = 0, 0, [], []
            for (idx, x) in enumerate(iter(X)):
                is_iter = False
                extra = tuple()
                if isinstance(x, collections.Iterable):
                    x, is_iter = list(x), True
                if is_iter and len(x) >= 0:
                    if len(x) == 0:
                        warnings.warn('Ignoring empty element on index: '
                                      + str(idx))
                        continue
                    elif len(x) == 1:
                        x = Graph(x[0], {}, {}, graph_format="adjacency")
                    elif len(x) == 2:
                        x = Graph(x[0], x[1], {}, graph_format="adjacency")
                    elif len(x) >= 3:
                        if len(x) > 3:
                            extra += tuple(x[3:])
                        x = Graph(x[0], x[1], x[2], graph_format="adjacency")
                elif type(x) is Graph:
                    x.desired_format("adjacency")
                    x = Graph(x.get_adjacency_matrix(),
                              x.get_labels(purpose="adjacency", label_type="vertex", return_none=True),
                              x.get_labels(purpose="adjacency", label_type="edge", return_none=True))
                else:
                    raise TypeError('each element of X must be either a '
                                    'graph object or a list with at least '
                                    'a graph like object and node labels '
                                    'dict \n')
                # workaround for leaving a sparse representation for x
                x.change_format(self._graph_format)
                c = core_number(x)
                max_core_number = max(max_core_number, max(c.values()))
                core_numbers.append(c)
                graphs.append((x, extra))

                nx += 1
            if nx == 0:
                raise ValueError('parsed input is empty')

        if max_core_number <= self.min_core:
            raise ValueError('The maximum core equals the min_core boundary set in init.')

        # Add the zero iteration element
        if self._method_calling == 2:
            K = np.zeros(shape=(nx, nx))
        elif self._method_calling == 3:
            self._dummy_kernel = dict()
            K = np.zeros(shape=(nx, self._nx))

        # Main
        base_graph_kernel, indexes_list = dict(), dict()
        for i in range(max_core_number, self.min_core, -1):
            subgraphs, indexes = list(), list()
            for (idx, (cn, (g, extra))) in enumerate(zip(core_numbers, graphs)):
                vertices = [k for k, v in iteritems(cn) if v >= i]
                if len(vertices) > 0:
                    # Calculate subgraph and store the index of the non-empty vertices
                    sg = g.get_subgraph(vertices)
                    sub_extra = list()
                    indexes.append(idx)
                    if len(extra) > 0:
                        vs = np.array(sg.get_vertices(purpose='any'))
                        for e in extra:
                            # This case will only be reached by now if the user add the propagation
                            # kernel as subkernel with a custom propagation matrix. This is a workaround!
                            if type(e) is np.array and len(e.shape) == 2:
                                e = e[vs, :][:, vs]
                            sub_extra.append(e)
                        subgraphs.append((sg, ) + tuple(sub_extra))
                    else:
                        subgraphs.append(sg)
            indexes = np.array(indexes)
            indexes_list[i] = indexes

            # calculate kernel
            if self._method_calling == 1 and indexes.shape[0] > 0:
                base_graph_kernel[i] = self.base_graph_kernel_(**self.params_)
                base_graph_kernel[i].fit(subgraphs)
            elif self._method_calling == 2 and indexes.shape[0] > 0:
                base_graph_kernel[i] = self.base_graph_kernel_(**self.params_)
                ft_subgraph_mat = base_graph_kernel[i].fit_transform(subgraphs)
                for j in range(indexes.shape[0]):
                    K[indexes[j], indexes] += ft_subgraph_mat[j, :]
            elif self._method_calling == 3:
                if self._max_core_number < i or self._fit_indexes[i].shape[0] == 0:
                    if len(indexes) > 0:
                        # add a dummy kernel for calculating the diagonal
                        self._dummy_kernel[i] = self.base_graph_kernel_(**self.params_)
                        self._dummy_kernel[i].fit(subgraphs)
                else:
                    if indexes.shape[0] > 0:
                        subgraph_tmat = self.X[i].transform(subgraphs)
                        for j in range(indexes.shape[0]):
                            K[indexes[j], self._fit_indexes[i]] += subgraph_tmat[j, :]

        if self._method_calling == 1:
            self._nx = nx
            self._max_core_number = max_core_number
            self._fit_indexes = indexes_list
            return base_graph_kernel
        elif self._method_calling == 2:
            self._nx = nx
            self._max_core_number = max_core_number
            self._fit_indexes = indexes_list
            return K, base_graph_kernel
        elif self._method_calling == 3:
            self._t_nx = nx
            self._max_core_number_trans = max_core_number
            self._transform_indexes = indexes_list
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

        self._is_transformed = True
        if self.normalize:
            X_diag, Y_diag = self.diagonal()
            old_settings = np.seterr(divide='ignore')
            km /= np.sqrt(np.outer(Y_diag, X_diag))
            km = np.nan_to_num(km)
            np.seterr(**old_settings)

        return km

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
        self._method_calling = 2
        self._is_transformed = False
        self.initialize()
        if X is None:
            raise ValueError('transform input cannot be None')
        else:
            km, self.X = self.parse_input(X)

        self._X_diag = np.diagonal(km)
        if self.normalize:
            old_settings = np.seterr(divide='ignore')
            km = np.nan_to_num(np.divide(km, np.sqrt(np.outer(self._X_diag, self._X_diag))))
            np.seterr(**old_settings)
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
            Y_diag = np.zeros(shape=(self._t_nx,))
            if self._is_transformed:
                max_core_number = min(self._max_core_number_trans, self._max_core_number)
                for i in range(max_core_number, self.min_core, -1):
                    tidx = self._transform_indexes[i]
                    if tidx.shape[0] > 0:
                        Y_diag[self._transform_indexes[i]] += self.X[i].diagonal()[1]
        except NotFittedError:
            # Calculate diagonal of X
            X_diag = np.zeros(shape=(self._nx,))
            if self._is_transformed:
                max_core_number = min(self._max_core_number_trans, self._max_core_number)
                Y_diag = np.zeros(shape=(self._t_nx,))
                for i in range(max_core_number, self.min_core, -1):
                    fidx = self._fit_indexes[i]
                    tidx = self._transform_indexes[i]
                    if tidx.shape[0] > 0 and fidx.shape[0] > 0:
                        x, y = self.X[i].diagonal()
                        X_diag[fidx] += x
                        Y_diag[tidx] += y
                if max_core_number < self._max_core_number:
                    for i in range(self._max_core_number, self._max_core_number_trans, -1):
                        fidx = self._fit_indexes[i]
                        if fidx.shape[0] > 0:
                            X_diag[fidx] += self.X[i].diagonal()
            else:
                for i in range(self._max_core_number, self.min_core, -1):
                    fidx = self._fit_indexes[i]
                    if fidx.shape[0] > 0:
                        X_diag[fidx] += self.X[i].diagonal()
            self._X_diag = X_diag
        if self._is_transformed:
            if len(self._dummy_kernel):
                for (idx, bk) in iteritems(self._dummy_kernel):
                    Y_diag[self._transform_indexes[idx]] += bk.diagonal()
            return self._X_diag, Y_diag
        else:
            return self._X_diag


def core_number(G):
    """Calculate the core number for each vertex.

    Parameters
    ----------
    G : grakel.Graph
        A graph type object corresponding to the input.

    Returns
    -------
    core : dict
        A dictionary containing for each node its core number.

    """
    nbrs, degrees = dict(), dict()
    for v in G.get_vertices(purpose='any'):
        ns = G.neighbors(v)
        nbrs[v] = ns
        degrees[v] = len(ns)
    nodes = sorted(degrees, key=degrees.get)
    bin_boundaries = [0]
    curr_degree = 0
    for i, v in enumerate(nodes):
        if degrees[v] > curr_degree:
            bin_boundaries.extend([i]*(degrees[v]-curr_degree))
            curr_degree = degrees[v]
    node_pos = dict((v, pos) for pos, v in enumerate(nodes))
    core = degrees
    for v in nodes:
        for u in nbrs[v]:
            if core[u] > core[v]:
                nbrs[u].remove(v)
                pos = node_pos[u]
                bin_start = bin_boundaries[core[u]]
                node_pos[u] = bin_start
                node_pos[nodes[bin_start]] = pos
                nodes[bin_start], nodes[pos] = nodes[pos], nodes[bin_start]
                bin_boundaries[core[u]] += 1
                core[u] -= 1
    return core


if __name__ == '__main__':
    from grakel.datasets import fetch_dataset
    import argparse
    # Create an argument parser for the installer of pynauty
    parser = argparse.ArgumentParser(
        description='Measuring classification accuracy '
                    ' on multiscale_laplacian_fast')

    parser.add_argument(
        '--dataset',
        help='choose the dataset you want the tests to be executed',
        type=str,
        default="MUTAG"
    )

    parser.add_argument(
        '--full',
        help='fit_transform the full graph',
        action="store_true")

    parser.add_argument(
        '--mc',
        help='the min_core kernel parameter',
        type=int,
        default=-1)

    # Get the dataset name
    args = parser.parse_args()
    dataset_name = args.dataset
    full = bool(args.full)
    mc = int(args.mc)
    # The baseline dataset for node/edge-attributes
    dataset_attr = fetch_dataset(dataset_name,
                                 with_classes=True,
                                 produce_labels_nodes=True,
                                 prefer_attr_nodes=False,
                                 verbose=True)

    from tqdm import tqdm
    from time import time

    from sklearn.metrics import accuracy_score
    from sklearn.model_selection import KFold
    from sklearn import svm
    from grakel.kernels import WeisfeilerLehman
    from grakel.kernels import VertexHistogram
    # from grakel.kernels import ShortestPath

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
        bk = (WeisfeilerLehman, dict(base_graph_kernel=VertexHistogram))
        # bk = (ShortestPath, dict(with_labels=False))
        # gk = WeisfeilerLehman(normalize=True, base_graph_kernel=VertexHistogram)
        gk = CoreFramework(normalize=True, base_graph_kernel=bk, min_core=mc)

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
    print("Core-Framework/WL/Subtree", "> Accuracy:",
          str(round(np.mean(stats["acc"])*100, 2)),
          "% | Took:", sec_to_time(np.mean(stats["time"])))
