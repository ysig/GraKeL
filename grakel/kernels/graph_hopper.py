"""The Graph Hopper kernel as defined in :cite:`feragen2013scalable`."""
# Author: Ioannis Siglidis <y.siglidis@gmail.com>
# License: BSD 3 clause
import numpy as np

from collections import defaultdict
from collections import Iterable
from numbers import Real
from warnings import warn
from numpy.matlib import repmat

from grakel.kernels import Kernel
from grakel.graph import Graph
from grakel.graph import dijkstra

# Python 2/3 cross-compatibility import
from six.moves import filterfalse


class GraphHopper(Kernel):
    """Graph Hopper Histogram kernel as found in :cite:`feragen2013scalable`.

    Parameters
    ----------
    kernel_type : str, tuple or function
        For `kernel_type` of **type**:
            + **str** can either be 'linear', 'gaussian', 'bridge'.
            + **tuple** can be of the form ('gaussian', mu) where mu is a number.
            + **function** can be a function that takes two tuples of np.arrays for each graph
              corresponding to the M matrix and the attribute matrix and returns a number.


    Attributes
    ----------
    metric_ : function
        The base metric applied between features.

    calculate_norm_ : bool
        Defines if the norm of the attributes will be calculated
        (in order to avoid recalculation when using it with e.g. gaussian).

    """

    _graph_format = "all"

    def __init__(self, n_jobs=None, normalize=False, verbose=False, kernel_type='linear'):
        """Initialize an Graph Hopper kernel."""
        super(GraphHopper, self).__init__(n_jobs=n_jobs,
                                          normalize=normalize,
                                          verbose=verbose)
        self.kernel_type = kernel_type
        self._initialized.update({"kernel_type": False})

    def initialize(self):
        """Initialize all transformer arguments, needing initialization."""
        super(GraphHopper, self).initialize()
        if not self._initialized["kernel_type"]:
            if type(self.kernel_type) is str:
                if self.kernel_type == "linear":
                    self.metric_ = linear_kernel
                    self.calculate_norm_ = False
                elif self.kernel_type == "gaussian":
                    self.metric_ = lambda x, y: gaussian_kernel(x, y, 1)
                    self.calculate_norm_ = True
                elif self.kernel_type == "bridge":
                    self.metric_ = bridge_kernel
                    self.calculate_norm_ = False
                else:
                    raise ValueError('Unsupported kernel with name "' + str(self.kernel_type) + '"')
            elif (type(self.kernel_type) is tuple and len(self.kernel_type) == 2 and
                    self.kernel_type[0] == "gaussian" and isinstance(self.kernel_type[1], Real)):
                self.metric_ = lambda x, y: gaussian_kernel(x, y, self.kernel_type[1])
                self.calculate_norm_ = True
            elif callable(self.kernel_type):
                self.metric_ = self._kernel_type
                self.calculate_norm_ = False
            else:
                raise TypeError('Unrecognized "kernel_type": can either be a str '
                                'from the supported: "linear", "gaussian", "bridge" '
                                'or tuple ("gaussian", mu) or a callable.')

    def parse_input(self, X):
        """Parse and check the given input for the Graph Hopper kernel.

        Parameters
        ----------
        X : iterable
            For the input to pass the test, we must have:
            Each element must be an iterable with at most three features and at
            least one. The first that is obligatory is a valid graph structure
            (adjacency matrix or edge_dictionary) while the second is
            node_labels and the third edge_labels (that fitting the given graph
            format).

        Returns
        -------
        out : np.array, shape=(len(X), n_labels)
            A np array for frequency (cols) histograms for all Graphs (rows).

        """
        if not isinstance(X, Iterable):
            raise TypeError('input must be an iterable\n')
        else:
            ni = 0
            diam = list()
            graphs = list()
            for (i, x) in enumerate(iter(X)):
                is_iter = False
                if isinstance(x, Iterable):
                    is_iter = True
                    x = list(x)

                if type(x) is Graph:
                    g = Graph(x.get_adjacency_matrix(),
                              x.get_labels(purpose="adjacency"),
                              {},
                              self._graph_format)
                elif is_iter and len(x) == 0 or len(x) >= 2:
                    if len(x) == 0:
                        warn('Ignoring empty element on index: '+str(i))
                        continue
                    elif len(x) >= 2:
                        g = Graph(x[0], x[1], {}, "adjacency")
                        g.change_format(self._graph_format)
                else:
                    raise TypeError('each element of X must be either a '
                                    'graph object or a list with at least '
                                    'a graph like object and node, ')

                spm, attr = g.build_shortest_path_matrix(labels="vertex")
                nv = g.nv()
                try:
                    attributes = np.array([attr[j] for j in range(nv)])
                except TypeError:
                    raise TypeError('All attributes of a single graph should have the same dimension.')
                diam.append(int(np.max(spm[spm < float("Inf")])))
                graphs.append((g.get_adjacency_matrix(), nv, attributes))
                ni += 1

        if self._method_calling == 1:
            max_diam = self._max_diam = max(diam) + 1
        else:
            max_diam = max(self._max_diam, max(diam) + 1)

        out = list()
        for i in range(ni):
            AM, node_nr, attributes = graphs[i]
            des = np.zeros(shape=(node_nr, node_nr, max_diam), dtype=int)
            occ = np.zeros(shape=(node_nr, node_nr, max_diam), dtype=int)

            # Convert adjacency matrix to dictionary
            idx_i, idx_j = np.where(AM > 0)
            ed = defaultdict(dict)
            for (a, b) in filterfalse(lambda a: a[0] == a[1], zip(idx_i, idx_j)):
                ed[a][b] = AM[a, b]

            for j in range(node_nr):
                A = np.zeros(shape=AM.shape)

                # Single-source shortest path from node j
                D, p = dijkstra(ed, j)

                D = np.array(list(D.get(k, float("Inf")) for k in range(node_nr)))
                p[j] = -1

                # Restrict to the connected component of node j
                conn_comp = np.where(D < float("Inf"))[0]

                # To-be DAG adjacency matrix of connected component of node j
                A_cc = A[conn_comp, :][:, conn_comp]

                # Adjacency matrix of connected component of node j
                AM_cc = AM[conn_comp, :][:, conn_comp]
                D_cc = D[conn_comp]
                conn_comp_converter = np.zeros(shape=(A.shape[0], 1), dtype=int)
                for k in range(conn_comp.shape[0]):
                    conn_comp_converter[conn_comp[k]] = k
                conn_comp_converter = np.vstack([0, conn_comp_converter])
                p_cc = conn_comp_converter[np.array(list(p[k] for k in conn_comp)) + 1]

                # Number of nodes in connected component of node j
                conncomp_node_nr = A_cc.shape[0]
                for v in range(conncomp_node_nr):
                    if p_cc[v] > 0:
                        # Generate A_cc by adding directed edges of form (parent(v), v)
                        A_cc[p_cc[v], v] = 1

                    # Distance from v to j
                    v_dist = D_cc[v]

                    # All neighbors of v in the undirected graph
                    v_nbs = np.where(AM_cc[v, :] > 0)[0]

                    # Distances of neighbors of v to j
                    v_nbs_dists = D_cc[v_nbs]

                    # All neighbors of v in undirected graph who are
                    # one step closer to j than v is; i.e. SP-DAG parents
                    v_parents = v_nbs[v_nbs_dists == (v_dist - 1)]

                    # Add SP-DAG parents to A_cc
                    A_cc[v_parents, v] = 1

                # Computes the descendants & occurence vectors o_j(v), d_j(v)
                # for all v in the connected component
                occ_p, des_p = od_vectors_dag(A_cc, D_cc)

                if des_p.shape[0] == 1 and j == 0:
                    des[j, 0, 0] = des_p
                    occ[j, 0, 0] = occ_p
                else:
                    # Convert back to the indices of the original graph
                    for v in range(des_p.shape[0]):
                        for l in range(des_p.shape[1]):
                            des[j, conn_comp[v], l] = des_p[v, l]
                    # Convert back to the indices of the original graph
                    for v in range(occ_p.shape[0]):
                        for l in range(occ_p.shape[1]):
                            occ[j, conn_comp[v], l] = occ_p[v, l]

            M = np.zeros(shape=(node_nr, max_diam, max_diam))
            # j loops through choices of root
            for j in range(node_nr):
                des_mat_j_root = np.squeeze(des[j, :, :])
                occ_mat_j_root = np.squeeze(occ[j, :, :])
                # v loops through nodes
                for v in range(node_nr):
                    for a in range(max_diam):
                        for b in range(a, max_diam):
                            # M[v,:,:] is M[v]; a = node coordinate in path, b = path length
                            M[v, a, b] += des_mat_j_root[v, b - a]*occ_mat_j_root[v, a]

            if self.calculate_norm_:
                out.append((M, attributes, np.sum(attributes ** 2, axis=1)))
            else:
                out.append((M, attributes))
        return out

    def pairwise_operation(self, x, y):
        """Graph Hopper kernel as proposed in :cite:`feragen2013scalable`.

        Parameters
        ----------
        x, y : tuple
            Extracted features from `parse_input`.

        Returns
        -------
        kernel : number
            The kernel value.

        """
        xp, yp = x[0], y[0]
        m = min(xp.shape[1], yp.shape[1])
        m_sq = m**2
        if x[0].shape[1] > m:
            xp = xp[:, :m, :][:, :, :m]
        elif y[0].shape[1] > m:
            yp = yp[:, :m, :][:, :, :m]

        return self.metric_((xp.reshape(xp.shape[0], m_sq),) + x[1:],
                            (yp.reshape(yp.shape[0], m_sq),) + y[1:])


def linear_kernel(x, y):
    """Graph Hopper linear pairwise kernel as proposed in :cite:`feragen2013scalable`.

    Parameters
    ----------
    x, y : tuple
        Extracted features from `parse_input`.

    Returns
    -------
    kernel : number
        The kernel value.

    """
    M_i, NA_i = x
    M_j, NA_j = y
    weight_matrix = np.dot(M_i, M_j.T)
    NA_linear_kernel = np.dot(NA_i, NA_j.T)
    return np.dot(weight_matrix.flat, NA_linear_kernel.flat)


def gaussian_kernel(x, y, mu):
    """Graph Hopper gaussian pairwise kernel as proposed in :cite:`feragen2013scalable`.

    Parameters
    ----------
    x, y : tuple
        Extracted features from `parse_input`.

    mu : Number
        The mean value of the gaussian.

    Returns
    -------
    kernel : number
        The kernel value.

    """
    M_i, NA_i, norm2_i = x
    M_j, NA_j, norm2_j = y
    weight_matrix = np.dot(M_i, M_j.T)
    NA_linear_kernel = np.dot(NA_i, NA_j.T)
    NA_squared_distmatrix = ((-2*NA_linear_kernel.T + norm2_i).T + norm2_j)
    nodepair = np.exp(-mu*NA_squared_distmatrix)
    return np.dot(weight_matrix.flat, nodepair.flat)


def bridge_kernel(x, y):
    """Graph Hopper bridge kernel as proposed in :cite:`feragen2013scalable`.

    Parameters
    ----------
    x, y : tuple
        Extracted features from `parse_input`.

    Returns
    -------
    kernel : number
        The kernel value.

    """
    M_i, NA_i = x
    M_j, NA_j = y
    weight_matrix = np.dot(M_i, M_j.T)
    NAs = np.vstack([NA_i, NA_j])
    NAs_linear_kernel = np.dot(NAs, NAs.T)
    NAs_distances = kernelmatrix2distmatrix(NAs_linear_kernel)
    NA_i_NA_j_distances = NAs_distances[:NA_i.shape[0], NA_i.shape[0]:]
    nodepair = (4-NA_i_NA_j_distances)/4
    nodepair[nodepair < 0] = 0
    return np.dot(weight_matrix.flat, nodepair.flat)


def kernelmatrix2distmatrix(K):
    """Convert a Kernel Matrix to a Distance Matrix.

    Parameters
    ----------
    K : np.array, n_dim=2
        The kernel matrix.

    Returns
    -------
    D : np.array, n_dim=2
        The distance matrix.

    """
    diag_K = K.diagonal().reshape(K.shape[0], 1)
    return np.sqrt(diag_K + diag_K.T - 2*K)


def od_vectors_dag(G, shortestpath_dists):
    """Compute the set of occurrence and distance vectors for G.

    Defined in :cite:`feragen2013scalable`.

    Parameters
    ----------
    G : np.array, n_dim=2
        DAG induced from a gappy tree where the indexing of nodes gives a
        breadth first order of the corresponding original graph

    shortestpath_dists : np.array, n_dim=1
        Shortest path distances from the source node.

    Returns
    -------
    occ : np.array, n_dim=2
        n x d descendant matrix occ, where n: `G.shape[0]` loops through the
        nodes of G, and d: 'diameter of G'. The rows of the occ matrix will be
        padded with zeros on the right.

    des : np.array, n_dim=2
        n x d descendant matrix des, where n: `G.shape[0]` loops through the
        nodes of G, and d: 'diameter of G'. The rows of the des matrix will be
        padded with zeros on the right.

    """
    dag_size = G.shape[0]
    DAG_gen_vector = shortestpath_dists + 1

    # This only works when the DAG is a shortest path DAG on an unweighted graph
    gen_sorted = DAG_gen_vector.argsort()
    re_sorted = gen_sorted.argsort()
    sortedG = G[gen_sorted, :][:, gen_sorted]
    delta = int(np.max(DAG_gen_vector))

    # Initialize:
    # For a node v at generation i in the tree, give it the vector
    # [0 0 ... 1 ... 0] of length h_tree with the 1 at the ith place.
    occ = np.zeros(shape=(dag_size, delta), dtype=int)
    occ[0, 0] = 1

    # Initialize:
    # For a node v at generation i in the tree, give it the vector
    # [0 0 ... 1 ... 0] of length delta with the 1 at the ith place.
    des = np.zeros(shape=(dag_size, delta), dtype=int)
    des[:, 0] = np.ones(shape=(1, dag_size))

    for i in range(dag_size):
        edges_starting_at_ith = np.where(np.squeeze(sortedG[i, :]) == 1)[0]
        occ[edges_starting_at_ith, :] = occ[edges_starting_at_ith, :] + \
            repmat(np.hstack([0, occ[i, :-1]]), edges_starting_at_ith.shape[0], 1)

        # Now use message-passing from the bottom of the DAG to add up the
        # edges from each node. This is easy because the vertices in the DAG
        # are depth-first ordered in the original tree; thus, we can just start
        # from the end of the DAG matrix.
        edges_ending_at_ith_from_end = np.where(np.squeeze(sortedG[:, dag_size - i - 1]) == 1)[0]
        des[edges_ending_at_ith_from_end, :] = (
            des[edges_ending_at_ith_from_end, :] +
            repmat(np.hstack([0, des[dag_size - i - 1, :-1]]),
                   edges_ending_at_ith_from_end.shape[0], 1))

    return occ[re_sorted, :], des[re_sorted, :]


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
        default="BZR"
    )

    parser.add_argument(
        '--full',
        help='fit_transform the full graph',
        action="store_true")

    mec = parser.add_mutually_exclusive_group()

    mec.add_argument(
        '--linear',
        help='choose a linear kernel',
        action="store_true")

    mec.add_argument(
        '--gaussian',
        help='choose a gaussian kernel (optionaly add a mu: default=1)',
        nargs='?',
        type=str,
        const='1',
        default=None)

    mec.add_argument(
        '--bridge',
        help='choose a bridge kernel',
        action="store_true")

    # Get the dataset name
    args = parser.parse_args()
    dataset_name = args.dataset

    if args.gaussian is not None:
        kernel_type = ('gaussian', float(args.gaussian))
    elif bool(args.bridge):
        kernel_type = 'bridge'
    else:
        kernel_type = 'linear'

    full = bool(args.full)
    # The baseline dataset for node/edge-attributes
    dataset_attr = fetch_dataset(dataset_name,
                                 with_classes=True,
                                 prefer_attr_nodes=True,
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
        gk = GraphHopper(normalize=True, kernel_type=kernel_type)

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
    print("GraphHopper", "> Accuracy:",
          str(round(np.mean(stats["acc"])*100, 2)),
          "% | Took:", sec_to_time(np.mean(stats["time"])))
