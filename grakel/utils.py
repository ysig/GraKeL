"""A file containing useful external functions"""
# Author: Ioannis Siglidis <y.siglidis@gmail.com>
# License: BSD 3 clause
import os
import numpy as np

from collections import defaultdict
from collections import Iterable

from sklearn.base import TransformerMixin
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import KFold
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import ShuffleSplit
from sklearn.base import BaseEstimator
from sklearn.svm import SVC
from sklearn.utils import Bunch
from sklearn.utils import check_random_state
from sklearn.utils.validation import check_is_fitted

from grakel import Graph
from grakel.graph import is_adjacency as valid_matrix


class KMTransformer(BaseEstimator, TransformerMixin):
    """A Kernel Matrix Transformer.

    Usefull for using precalculated Kernel Matrices inside scikit-learn pipeline.

    Parameters
    ----------
    K : array-like, shape=[n, n]
        If given an array the input can be as follows:

            + array-like lists of lists

            + np.array

            + sparse matrix (scipy.sparse)

        It can also be embedded in an sklearn Bunch object as a mat (argument)

    Attributes
    ----------
    K_ : numpy.array, shape=[n, n]

    """
    def __init__(self, K=None):
        """Initialise the Kernel Matrix Transformer"""

        self.K = K
        self._initialized = {"K": False}

    def initialize(self):
        """Initialize all transformer arguments, needing initialisation."""
        if not self._initialized["K"]:
            if self.K is None:
                M = np.array([[1.0]])
            else:
                K = self.K
                if isinstance(K, Bunch):
                    try:
                        K = K.mat
                    except Exception:
                        raise ValueError('If in an sklearn Bunch K must be under mat')
                flag, M = valid_matrix(K, transform=True)
                if not flag:
                    raise ValueError('The provided K cannot be converted to a '
                                     'two dimensional np.array.')
            self.K_ = M
            self._initialized["K"] = True

    def fit(self, X, y=None):
        """Fit a list of indeces.

        Parameters
        ----------
        X : Iterable of int:
            Indexes for the X of the kernel matrix.

        Returns
        -------
        self : object
            Returns self.

        """
        self.initialize()
        if any(x < 0 or x > self.K_.shape[0] for x in X):
            raise ValueError('')
        else:
            self.X = np.array(X)

        return self

    def fit_transform(self, X, y=None):
        """Fit and transform, on the same dataset.

        Parameters
        ----------
        X : Iterable of int
            Indexes for the first dimension of the kernel matrix.

        y : None
            There is no need of a target in a transformer, yet the pipeline API
            requires this parameter.

        Returns
        -------
        K : numpy array, shape = [len(X), len(X)]
            Corresping to the values of the X indexes with themselfs.

        """
        # Initialize the Graph Kernel
        self.initialize()
        if any(x < 0 or x > self.K_.shape[0] for x in X):
            raise ValueError('')
        else:
            self.X = np.array(X)

        return self.K_[self.X, :][:, self.X]

    def transform(self, X):
        """Calculate the kernel matrix, between given and fitted dataset.

        Parameters
        ----------
        X : Iterable of int
            Indexes for the second dimension of the kernel matrix.

        Returns
        -------
        K : numpy array, shape = [len(Y), len(X)]
            Corresping to the values of the Y indexes with X.

        """
        check_is_fitted(self, 'X')
        if any(x < 0 or x > self.K_.shape[0] for x in X):
            raise ValueError('')

        return self.K_[X, :][:, self.X]


def cross_validate_Kfold_SVM(K, y,
                             n_iter=10, n_splits=10, C_grid=None,
                             random_state=None, scoring="accuracy", fold_reduce=None):
    """Cross Validate a list of precomputed kernels with an SVM.

    Parameters
    ----------
    K : list
        A list that must contain either numpy arrays or iterables of numpy arrays.

    y : list
        List of lists that for every element of K contains numbers of score for all iterations.

    n_iter : int
        Number of iteration for the K-Fold.

    n_splits : int
        Number of splits for the K-Fold.

    random_state :  RandomState or int, default=None
        A random number generator instance or an int to initialize a RandomState as a seed.

    fold_reduce : callable or None
        A function that summarizes information between all folds.
        Input must be a list of n_splits elements corresponding to scoring.
        If None default is np.mean.

    scoring : string, callable, list/tuple, dict or None, default: None
        As in ``scoring`` in :xref:`lgscv`.

    Returns
    -------
    out : list
       A list that contains a list of the reduced folds for each iteration, for each primary
       element of K.
    """
    # Initialise C_grid
    if C_grid is None:
        C_grid = ((10. ** np.arange(-7, 7, 2)) / len(y)).tolist()
    elif type(C_grid) is np.array:
        C_grid = np.squeeze(C_grid)
        if len(C_grid.shape) != 1:
            raise ValueError('C_grid should either be None or a squeezable to 1 dimension np.array')
        else:
            C_grid = list(C_grid)

    # Initialise fold_reduce:
    if fold_reduce is None:
        fold_reduce = np.mean
    elif not isinstance(callable, fold_reduce):
        raise ValueError('fold_reduce should be a callable')

    # Initialise and check random state
    random_state = check_random_state(random_state)

    # Initialise sklearn pipeline objects
    kfolder = KFold(n_splits=n_splits, random_state=random_state, shuffle=True)
    estimator = make_pipeline(KMTransformer(), SVC(kernel='precomputed'))

    # Make all the requested folds
    nfolds = tuple(tuple(kfolder.split(y)) for _ in range(n_iter))

    out = list()
    for ks in K:
        mid = list()
        if valid_matrix(ks):
            pg = {"svc__C": C_grid, "kmtransformer__K": [Bunch(mat=ks)]}
        elif isinstance(ks, Iterable) and all(valid_matrix(k) for k in ks):
            pg = [{"svc__C": C_grid, "kmtransformer__K": [Bunch(mat=k)]} for k in ks]
        else:
            raise ValueError('Not a valid object for kernel matrix/ces')

        for kfolds in nfolds:
            fold_info = list()
            for train, test in kfolds:
                gs = GridSearchCV(estimator, param_grid=pg, scoring=scoring,
                                  cv=ShuffleSplit(n_splits=1,
                                                  test_size=0.1,
                                                  random_state=random_state)).fit(train, y[train])
                fold_info.append(gs.score(test, y[test]))
            mid.append(fold_reduce(fold_info))
        out.append(mid)
    return out


def graph_from_networkx(X, node_labels_tag=None, edge_labels_tag=None, edge_weight_tag=None,
                        as_Graph=False):
    """Transform an iterable of networkx objects to an iterable of Graphs.

    A function for helping a user that has a collection of graphs in networkx to use grakel.

    Parameters
    ----------
        node_labels_tag : str or None
            Define where to search for labels of nodes, inside the `node` attribute of each graph.
            If None no labels are assigned.

        edge_labels_tag : tuple
            Define where to search for labels of edges, inside the `edge` attribute of each graph.
            If None no labels are assigned.

        edge_weight_tag : int
            Define where to search for weights inside the `edge` attribute of each graph.
            If None 1.0 weights are assigned.

        as_Graph : bool, default=False
            Return each output as a grakel.Graph object.

    Returns
    -------
        grakel_graphs : generator
            Returns a generator containing the same collections to be given as an input for
            any grakel kernel.

    """
    import networkx as nx
    v2 = False
    if nx.__version__ > '2.':
        v2 = True

    if node_labels_tag is None:
        def nodel_init():
            return None

        def nodel_put(*args):
            pass
    elif type(node_labels_tag) is str:
        def nodel_init():
            return dict()

        def nodel_put(nl, u, d):
            nl[u] = d[u][node_labels_tag]
    else:
        raise ValueError('node_labels_tag must be a str indicating the '
                         'tag of the labels inside nodes or None')

    if edge_labels_tag is None:
        def edgel_init():
            return None

        def edgel_put(*args):
            pass
    elif type(edge_labels_tag) is str:
        def edgel_init():
            return dict()

        if v2:
            def edgel_put(el, u, d):
                el[u] = d[u][edge_labels_tag]
        else:
            def edgel_put(el, u, d):
                el[u] = d[u[0]][u[1]][edge_labels_tag]
    else:
        raise ValueError('edge_labels_tag must be a str indicating the '
                         'tag of the labels inside edges or None')

    if edge_weight_tag is None:
        def get_weight(*args):
            return 1.0
    elif type(edge_weight_tag) is str:
        if v2:
            def get_weight(d, e):
                return d[e][edge_weight_tag]
        else:
            def get_weight(d, e):
                return d[e[0]][e[1]][edge_weight_tag]

    else:
        raise ValueError('weight_labels_tag must be a str indicating  '
                         'tag of the labels inside edges or None (1.0)')

    if not isinstance(X, Iterable):
        raise ValueError('X must be an iterable')

    if v2:
        def take_ne(graph):
            return graph.nodes, graph.edges
    else:
        def take_ne(graph):
            return graph.node, graph.edge

    for G in X:
        graph_object = dict()
        nl = nodel_init()
        el = edgel_init()
        nodes, edges = take_ne(G)
        for u in G.nodes():
            graph_object[u] = dict()
            nodel_put(nl, u, nodes)
            for v in G.neighbors(u):
                graph_object[u][v] = get_weight(edges, (u, v))
                edgel_put(el, (u, v), edges)

        if as_Graph:
            yield Graph(graph_object, nl, el)
        else:
            yield [graph_object, nl, el]


def graph_from_pandas(edge_df, node_df=None, directed=False, as_Graph=False):
    """Produces a collection of Graph Objects from pandas dataframes.

    A function for helping a user that has a bunch of graph csv
    files make a collection for grakel input.

    Parameters
    ----------
        edge_df : tuple(pandas.DataFrame, pct, tuple(pct, pct), pct or None, pct or list(pct) or None)
            A tuple of many elements:

                1. A pandas dataframe containing all the edges of a collection of graphs
                2. The column name of the graph index
                3. A tuple containing source and destination column names
                4. The column name of weights column (None if non-existent)
                5. The column name pointing the column containing edge-labels or a list pointing
                   the columns of attributes (None if edge-labels are non-existent).

            If node_df exists, correspondance of graph indexes and node indexes must be exact.

        node_df : tuple(pandas.DataFrame, pct, pct or list(pct) or None) or None, default=None
            A tuple of many elements:

                1. A pandas dataframe containing all the nodes of a collection of graphs
                2. The column name of the graph index
                3. The column name pointing the column containing node-labels or a list pointing
                   the columns of attributes (None if node-labels are non-existent).

            Node id must correspond to node number.

        directed : bool, default=False
            A variable indicating with True if both directions of the graph have been inserted inside
            the dataframe. Otherwise the graph is considered as undirected.

        as_Graph : bool, default=False
            Return each output as a grakel.Graph object.

    Returns
    -------
        grakel_graphs : dict
            Dictionary with keys te keys of the graphs inside the pandas dataframe and values the
            corresponding graph object

    """
    from pandas import DataFrame

    if node_df is None:
        graphs = dict()
        throw_error = False
    elif (type(node_df) is tuple and
          len(node_df) == 3 and
          type(node_df[0]) is DataFrame and
          (node_df[1] is None or node_df[1] in node_df[0].columns) and
          (node_df[2] is None or
           (type(node_df[2]) is not list and node_df[2] in node_df[0].columns) or
           (type(node_df[2]) is list and all(i in node_df[0].columns for i in node_df[2])))):
        df, gtag, labs = node_df
        graphs = defaultdict(lambda: defaultdict(dict))
        if labs is None:
            for index, row in df.iterrows():
                gidx = row[gtag]
                graphs[gidx]["graph"][index] = dict()
                graphs[gidx]["node_label"] = None
        elif type(labs) is list:
            for index, row in df.iterrows():
                gidx = row[gtag]
                graphs[gidx]["graph"][index] = dict()
                graphs[gidx]["node_label"][index] = np.array(row[c] for c in labs)
        else:
            for index, row in df.iterrows():
                gidx = row[gtag]
                graphs[gidx]["graph"][index] = dict()
                graphs[gidx]["node_label"][index] = row[labs]
        throw_error = True
    else:
        raise ValueError('node_df must be a tuple containing a pandas.DataFrame object '
                         'a column name corresponding to the valid column index for graphs indexes '
                         'and a None, column name or list of column names corresponding to no-labels, '
                         'labels or attributes.')

    if (type(edge_df) is tuple and
        len(edge_df) == 5 and
        type(edge_df[0]) is DataFrame and
        (edge_df[1] in edge_df[0].columns) and
        (type(edge_df[2]) is tuple and len(edge_df[2]) == 2 and
         all(c in edge_df[0].columns for c in edge_df[2])) and
        (edge_df[3] is None or edge_df[3] in edge_df[0].columns) and
        (edge_df[4] is None or
         (type(edge_df[4]) is not list and edge_df[4] in edge_df[0].columns) or
         (type(edge_df[4]) is list and all(c in edge_df[0].columns for c in edge_df[4])))):
        df, gtag, (src_c, dst_c), w_c, labs = edge_df
        if w_c is not None:
            def get_weight(row):
                return row[w_c]
        else:
            def get_weight(row):
                return 1.

        if labs is None:
            def set_label(*args):
                pass
        elif type(labs) is list:
            def set_label(row, d, e):
                if d["edge_label"] is None:
                    d["edge_label"] = dict()
                el = np.array([row[l] for l in labs])
                d["edge_label"][e] = el
                if not directed:
                    d["edge_label"][(e[1], e[0])] = el
        else:
            def set_label(row, d, e):
                if d["edge_label"] is None:
                    d["edge_label"] = dict()
                el = row[labs]
                d["edge_label"][e] = el
                if not directed:
                    d["edge_label"][(e[1], e[0])] = el

        for index, row in df.iterrows():
            src, dst = row[src_c], row[dst_c]
            gidx = row[gtag]
            if gidx not in graphs:
                if throw_error:
                    raise ValueError('This graph didn\'t appear on node labels dataframe')
                else:
                    graphs[gidx] = dict()
                    graphs[gidx]["graph"] = dict()
                    if "node_label" not in graphs[gidx]:
                        graphs[gidx]["node_label"] = None
                    graphs[gidx]["edge_label"] = None
            if src not in graphs[gidx]["graph"]:
                if throw_error:
                    raise ValueError('This node didn\'t appear on node labels dataframe for graph '
                                     'with id ' + str(gidx))
                else:
                    graphs[gidx]["graph"][src] = dict()
            graphs[gidx]["graph"][src][dst] = get_weight(row)
            if not directed:
                if dst not in graphs[gidx]["graph"]:
                    if throw_error:
                        raise ValueError('This node didn\'t appear on node labels dataframe '
                                         'for graph with id ' + str(gidx))
                    else:
                        graphs[gidx]["graph"][dst] = dict()
                graphs[gidx]["graph"][dst][src] = get_weight(row)
            set_label(row, graphs[gidx], (src, dst))
    else:
        raise ValueError('edge_df must be a tuple containing a pandas.DataFrame object '
                         'a column name corresponding to a valid column index of the dataframe '
                         'that contain the index of the graph that an edge belongs inside the '
                         'dataframe, a column name corresponding to weights inside the dataframe'
                         '(or None if weights do not exist) and a list of column names '
                         'corresponding or column name corresponding to labels or None '
                         'if none of this exists')

    return {k: Graph(graphs[k]["graph"], graphs[k]["node_label"], graphs[k]["edge_label"])
            if as_Graph else [graphs[k]["graph"], graphs[k]["node_label"], graphs[k]["edge_label"]]
            for k in graphs.keys()}


def graph_from_csv(edge_files, node_files=None, index_type=str,
                   directed=False, sep=",", as_Graph=False):
    """Produces a collection of Graph Objects from a collection of csv files.

    A function for helping a user that has a bunch of graph csv
    files make a collection for grakel input.

    Parameters
    ----------
        edge_files : tuple(iter(str), bool, bool or None)
            edge_files = (iter(<edge_file_address>), <weight_flag>, <attributes_flag>).
            Each line of all the edge files must have the following structure:

            .. code::

                <v_i><sep><v_j>[weight][labels/attr]

                weights: "<sep><weight>", <weight_flag> is True
                         "", <weight_flag> is False
                labels/attr: "<sep><label>", <attributes_flag> is False
                             "<sep><a_1>, .. , <sep><a_n>", <attributes_flag> is True
                             "", <attributes_flag> is None

        node_files : tuple(iter(str), bool) or None, default=None
            node_files = (iter(<node_file_address>), index_type, <attributes_flag>).
            Each line of the node file must have the following structure (if exists):

            .. code::

                <v_i><sep>[labels/attr]

                labels/attr: "<sep><label>", <attributes_flag> is False
                             "<sep><a_1>, .. , <sep><a_n>", <attributes_flag> is True
                             "", <attributes_flag> is None

        directed : bool, default=False
            Defines if the graph should be considered as directed.

        sep : str, default=","
            The separator for the csv files.

        as_Graph : bool, default=False
            Return each output as a grakel.Graph object.

    Returns
    -------
        grakel_graphs : generator
            Returns a generator containing the same collections to be given as an input for
            any grakel kernel.

    """
    def edge_files_error():
        raise ValueError('edge_file argument must contain an iterable of strings of edge files, '
                         'a bool weight_flag and attributes_flag bool or None')

    def node_files_error():
        raise ValueError('node_files argument can be None or contain an iterable of strings'
                         'of edge files and attributes_flag bool or None')

    if type(index_type) is not type:
        raise ValueError('index_type must be a class `type` object')

    if (type(edge_files) is not tuple or len(edge_files) != 3 or
            type(edge_files[1]) is not bool or
            type(edge_files[2]) not in [bool, None]):
        edge_files_error()
    else:
        if edge_files[1]:
            def get_weight(l):
                return float(l.pop(0))
        else:
            def get_weight(l):
                return 1.0
        efs = list()
        for e in edge_files[0]:
            if type(e) is not str:
                edge_files_error()
            if not os.path.isfile(e):
                raise ValueError('Each edge file address must be a valid address')
            efs.append(e)

    if type(node_files) is None:
        nfs = None
    elif (type(node_files) is tuple and len(node_files) == 2 and
            type(node_files[1]) in [bool, None]):
        nfs = list()
        for n in node_files[0]:
            if type(n) is not str:
                edge_files_error()
            if not os.path.isfile(n):
                raise ValueError('Each edge file address must be a valid address')
            nfs.append(n)
    else:
        node_files_error()

    if edge_files[2] is None:
        def get_el(ef, g):
            el = None
            for line in ef:
                q = line.split(sep)
                m = q[2:]
                ea, eb = [index_type(e) for e in q[:2]]
                g[ea][eb] = get_weight(m)
                if not directed:
                    g[eb][ea] = g[ea][eb]
            return el
    elif edge_files[2]:
        def get_el(ef, g):
            el = dict()
            for line in ef:
                q = line.strip('\n').split(sep)
                ea, eb = [index_type(e) for e in q[:2]]
                m = q[2:]
                g[ea][eb] = get_weight(m)
                el[(ea, eb)] = np.array([float(num) for num in m])
                if not directed:
                    el[(eb, ea)] = el[(ea, eb)]
                    g[eb][ea] = g[ea][eb]
            return el
    else:
        def get_el(ef, g):
            el = dict()
            for line in ef:
                q = line.strip('\n').split(sep)
                ea, eb = [index_type(e) for e in q[:2]]
                m = q[2:]
                g[ea][eb] = get_weight(m)
                el[(ea, eb)] = str(m[0])
                if not directed:
                    el[(eb, ea)] = el[(ea, eb)]
                    g[eb][ea] = g[ea][eb]
            return el

    if nfs is None:
        for efa in efs:
            nl = None
            graph_object = defaultdict(dict)
            with open(efa, 'r') as ef:
                el = get_el(ef, graph_object)
            graph_object = dict(graph_object)

            if as_Graph:
                yield Graph(graph_object, nl, el)
            else:
                yield [graph_object, nl, el]
    else:
        if len(nfs) != len(efs):
            raise ValueError('The number of edge files and the number of node files must be the same')
        for (nfa, efa) in zip(nfs, efs):
            graph_object = dict()
            with open(nfa, 'r') as nf:
                if node_files[1] is None:
                    nl = None
                    for line in nf:
                        graph_object[index_type(line.strip('\n').split(sep)[0].strip())] = dict()
                elif node_files[1]:
                    nl = dict()
                    for line in nf:
                        q = line.strip('\n').split(sep)
                        graph_object[index_type(q[0].strip())] = dict()
                        nl[index_type(q[0].strip())] = np.array([float(i.strip()) for i in q[1:]])
                else:
                    nl = dict()
                    for line in nf:
                        q = line.strip('\n').split(sep)
                        graph_object[index_type(q[0].strip())] = dict()
                        nl[index_type(q[0].strip())] = q[1]

            with open(efa, 'r') as ef:
                el = get_el(ef, graph_object)

            if as_Graph:
                yield Graph(graph_object, nl, el)
            else:
                yield [graph_object, nl, el]
