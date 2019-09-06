"""The ODD-Sth kernel as defined in :cite:`Martino2012ATK`."""
# Author: Ioannis Siglidis <y.siglidis@gmail.com>
# License: BSD 3 clause
import copy
import warnings

import numpy as np

from collections import Iterable
from collections import defaultdict

from sklearn.utils.validation import check_is_fitted
from sklearn.exceptions import NotFittedError

from grakel.kernels import Kernel
from grakel.graph import Graph

# Python 2/3 cross-compatibility import
from six import iteritems


class OddSth(Kernel):
    """ODD-Sth kernel as proposed in :cite:`Martino2012ATK`.

    Parameters
    ----------
    h : int, default=None
        Maximum (single) dag height.
        If None there is no restriction.

    Attributes
    ----------
    _nx : int
        The number of parsed inputs on fit.

    _ny : int
        The number of parsed inputs on transform.

    _phi_x : np.array, n_dim=2
        A numpy array corresponding all the frequency values for
        each vertex, coresponding to the fitted data, in the
        resulting bigDAG of the fitted and transformed data.

    _phi_y : np.array, n_dim=2
        A numpy array corresponding all the frequency values for
        each vertex, corresponding to the transformed data, in the
        resulting bigDAG of the fitted and transformed data.

    _C : np.array, n_dim=1
        A numpy array corresponding to the vertex depth of each node, in
        the resulting bigDAG of the fitted and transformed data.

    """

    def __init__(self, n_jobs=None, normalize=False, verbose=False, h=None):
        """Initialise an `odd_sth` kernel."""
        super(OddSth, self).__init__(n_jobs=n_jobs,
                                     normalize=normalize,
                                     verbose=verbose)
        self.h = h
        self._initialized.update({"h": False})

    def initialize(self):
        """Initialize all transformer arguments, needing initialization."""
        if not self._initialized["n_jobs"]:
            if self.n_jobs is not None:
                warnings.warn('no implemented parallelization for OddSth')
            self._initialized["n_jobs"] = True

        if not self._initialized["h"]:
            if self.h is not None and (type(self.h) is not int or self.h <= 0):
                raise ValueError('h must be an integer bigger than zero')

            self.h_ = (-1 if self.h is None else self.h)
            self._initialized["h"] = True

    def parse_input(self, X):
        """Parse and create features for the propagation kernel.

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
        out : tuple
            A tuple corresponding to the calculated bigDAG.

        """
        if not isinstance(X, Iterable):
            raise TypeError('input must be an iterable\n')
        else:
            i = 0
            out = None
            if self._method_calling == 3:
                out = copy.deepcopy(self.X)
            for (idx, x) in enumerate(iter(X)):
                is_iter = isinstance(x, Iterable)
                if is_iter:
                    x = list(x)
                if is_iter and (len(x) == 0 or len(x) >= 2):
                    if len(x) == 0:
                        warnings.warn('Ignoring empty element' +
                                      ' on index: '+str(idx))
                        continue
                    elif len(x) >= 2:
                        x = Graph(x[0], x[1], {}, self._graph_format)
                elif type(x) is not Graph:
                    raise TypeError('each element of X must have either ' +
                                    'a graph with labels for node and edge ' +
                                    'or 3 elements consisting of a graph ' +
                                    'type object, labels for vertices and ' +
                                    'labels for edges.')
                out = big_dag_append(make_big_dag(x, self.h_), out, merge_features=False)
                i += 1

            if self._method_calling == 1:
                self._nx = i
            elif self._method_calling == 3:
                self._ny = i
            if i == 0:
                raise ValueError
                ('parsed input is empty')
            return out

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

        y : Object, default=None
            Ignored argument, added for the pipeline.

        Returns
        -------
        K : numpy array, shape = [_nx, _nx]
            corresponding to the kernel matrix, a calculation between
            all pairs of graphs between target an features

        """
        self._method_calling = 2
        self.fit(X)

        # calculate feature matrices.
        ref = dict(enumerate(self.X[0].keys()))

        C = np.empty(shape=(len(ref), 1))
        phi = np.empty(shape=(len(ref), self._nx))
        for (i, v) in iteritems(ref):
            # number of identical subtrees
            # equal the D element
            C[i] = self.X[0][v][0]
            for j in range(self._nx):
                phi[i, j] = self.X[0][v][1][j]

        km = np.dot(phi.T, np.multiply(phi, C))

        self._X_diag = np.diagonal(km)
        if self.normalize:
            return np.divide(km, np.sqrt(np.outer(self._X_diag, self._X_diag)))
        else:
            return km
        return km

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
            full_dag = self.parse_input(X)

        ref = dict(enumerate(full_dag[0].keys()))

        C = np.empty(shape=(len(ref), 1))
        phi_x = np.empty(shape=(len(ref), self._nx))
        phi_y = np.empty(shape=(len(ref), self._ny))
        for (i, v) in enumerate(ref.keys()):
            # number of identical subtrees equal the D element
            C[i] = full_dag[0][v][0]
            for j in range(self._nx):
                phi_x[i, j] = full_dag[0][v][1][j]
            for j in range(self._ny):
                phi_y[i, j] = full_dag[0][v][1][j + self._nx]

        self._phi_x, self._phi_y, self._C = phi_x, phi_y, C
        km = np.dot(phi_y.T, np.multiply(phi_x, C))
        self._is_transformed = True
        if self.normalize:
            X_diag, Y_diag = self.diagonal()
            km /= np.sqrt(np.outer(Y_diag, X_diag))
        return km

    def diagonal(self):
        """Calculate the kernel matrix diagonal of the fitted data.

        Parameters
        ----------
        None.

        Returns
        -------
        X_diag : np.array
            The diagonal of the kernel matrix, of the fitted. This consists
            of each element calculated with itself.


        """
        # Check is fit had been called
        check_is_fitted(self, ['_phi_x', '_C'])
        try:
            check_is_fitted(self, ['_X_diag'])
        except NotFittedError:
            # Calculate diagonal of X
            self._X_diag = np.dot(np.square(self._phi_X),
                                  self._C).reshape((self._nx, 1))

        try:
            check_is_fitted(self, ['_phi_y'])
            Y_diag = np.dot(np.square(self._phi_y).T,
                            self._C).reshape((self._ny, 1))
            return self._X_diag, Y_diag
        except NotFittedError:
            return self._X_diag


def make_big_dag(g, h):
    """Compose a big dag out of all dags of a graph.

    Parameters
    ----------
    g : graph
        A graph type object.

    h : int, default=None
        Maximum (single) dag height.
        If None there is no restriction.

    Returns
    -------
    big_dag : tuple
        The big dag tuple consisting of:
            + A dictionary for vertices containing at each value the frequency
              the depth and the ID
            + A hash map from each ID to a list of vertices.
            + A list of sorted vertices based on ordering
            + A dictionary of edges.
            + A dictionary of labels

    """
    big_dag = None
    for v in g.get_vertices(purpose='any'):
        dag_odd = make_dag_odd(v, g, h)
        dag = tuple(hash_trees(dag_odd))+tuple([])+(dag_odd[1], dag_odd[3])
        big_dag = big_dag_append(dag, big_dag)

    _, D_edges, D_ordering, _ = odd(big_dag[0], big_dag[2], big_dag[3])

    big_dag = (big_dag[0],
               big_dag[1],
               sorted(list(big_dag[0].keys()),
               key=lambda x:
                      (D_ordering[x], big_dag[3][x])), D_edges, big_dag[3])

    return big_dag


def make_dag_odd(v, g, h):
    """Calculate the vertex rooted dag and apply inverse topological sorting.

    Parameters
    ----------
    v : hashable
        The dag vertex root.

    g : graph
        A graph type object from where v is coming from.

    h : int, default=None
        Maximum depth of the exploration.

    Returns
    -------
    odd_dag : tuple
        A a tuple representing the topologically sorted dag:
            + A set of vertices
            + A dictionary of sorted edges for each node based on ordering
              and labels
            + A dictionary of ordering for each node
            + A dictionary of labels for each node

    """
    vertices, edges = dag(v, g, h)
    return odd(vertices, edges, g.get_labels(purpose='any'))


def dag(v, g, h):
    """BFS exploration that returns a dag.

    Parameters
    ----------
    v : hashable
        The dag vertex root.

    g : graph
        A graph type object from where v is coming from.

    h : int, default=None
        Maximum depth of the exploration.

    Returns
    -------
    vertices : set
        The dag vertices.

    edges : dict
        The dictionary of edges. For each node a list of vertices.

    """
    q = [(v, 0)]
    vertices = dict()
    edges = defaultdict(list)

    vertices[v] = 0
    while len(q) > 0:
        u, level = q.pop(0)

        if level == h:
            break

        for n in g.neighbors(u, purpose='any'):
            if n not in vertices:
                edges[u].append(n)
                q.append((n, level+1))
                vertices[n] = level+1
            elif vertices[n] >= level+1:
                edges[u].append(n)

    vertices = set(vertices.keys())
    return vertices, edges


def odd(vertices, edges, labels):
    """Calculate the inverse topological order of a DAG and sorts it's edges.

    Parameters
    ----------
    vertices : dict or set
        A set of vertices.

    edges : dict
        Edges between vertices.

    labels : dict
        Labels for each vertex.

    Returns
    -------
    vertices : dict or set
        A dictionary or a set of vertices depending on the given input.

    edges : dict
        An edges dictionary, where for each edge the list of adjacent vertex
        is sorted corresponding to the given ordering or labels.

    ordering : dict
        The inverse topological ordering for each vertex.

    labels : dict
        The labels dictionary for each vertex.

    """
    # Kahn's algorithm for topological sorting of dags

    # Step-1: Compute in-degree (number of incoming edges)
    # for each of the vertex present in the DAG and initialize
    # the count of visited nodes as 0.
    indegrees = dict()

    if type(vertices) is set:
        zero_indegrees = vertices.copy()
        visited_nodes = len(vertices)
    elif type(vertices) is dict:
        zero_indegrees = set(vertices.keys())
        visited_nodes = len(vertices.keys())
    else:
        raise TypeError('unsupported vertices type')

    for (k, e) in iteritems(edges):
        for v in e:
            if v not in indegrees:
                indegrees[v] = 1
            else:
                indegrees[v] += 1
            zero_indegrees.discard(v)

    # Step-2: Pick all the vertices with in-degree as
    # 0 and add them into a queue
    q = list(zero_indegrees)
    ordering = dict()
    while len(q) > 0:
        # Step-3: Remove a vertex from the queue and then
        # Increment count of visited nodes by 1.
        # Decrease in-degree by 1 for all its neighboring nodes.
        # If in-degree of a neighboring nodes is reduced to zero,
        # then add it to the queue.

        q.sort(key=lambda x: labels[x])
        e = q.pop(0)
        ordering[e] = visited_nodes
        for k in edges[e]:
            if k in indegrees:
                if indegrees[k] == 1:
                    indegrees.pop(k)
                    q.append(k)
                else:
                    indegrees[k] -= 1
        visited_nodes -= 1

    # apply the ordering
    for k in edges.keys():
        edges[k].sort(key=lambda x: (ordering[x], labels[x]))

    return vertices, edges, ordering, labels


def hash_trees(tree):
    """Hashes trees and adds frequencies and a hash map.

    Parameters
    ----------
    tree : tuple
        A tuple of elements corresponding to a tree:
            + A set of vertices
            + A dictionary of edges
            + A dictionary that corresponds to an ordering of vertices
            + A dictionary of labels

    Returns
    -------
    hash_tree : tuple
        A hashed version of the tree:
            + A dictionary from vertices to tuples containing the subtree size,
              frequencies and node ID
            + A dictionary between hashes and vertices
              (representing a valid hash map)
            + A list of ordered vertices, based on an the inverse
              topological ordering

    """
    (vertex, edge, ordering, labels) = tree
    v_ordered = sorted(list(vertex), key=lambda x: (ordering[x], labels[x]))
    vertices_hash_map = dict()
    vertices = dict()
    for v in v_ordered:
        if v not in edge or len(edge[v]) == 0:
            if labels[v] not in vertices_hash_map:
                vertices_hash_map[labels[v]] = list()
            vertices_hash_map[labels[v]].append(v)
            vertices[v] = [0, 1, str(labels[v])]
        else:
            neighbors_ids = []
            d = 0
            for n in edge[v]:
                d += 1 + vertices[n][0]
                neighbors_ids.append(vertices[n][2])

            ID = str(labels[v])+'('+(','.join(neighbors_ids))+')'

            vertices[v] = [d, 1, ID]
            if ID not in vertices_hash_map:
                vertices_hash_map[ID] = list()
            vertices_hash_map[ID].append(v)

    return (vertices, vertices_hash_map, v_ordered)


def big_dag_append(dag, big_dag=None, merge_features=True):
    """Calculate the *minimal DAG* or *BigDAG*.

    See :cite:`Martino2006` and notated in :cite:`Martino2012ATK`.

    Parameters
    ----------
    dag : tuple
        A tuple representing a single dag containing:
            + A dictionary from vertices to tuples containing the subtree size,
              frequencies and node ID
            + A dictionary between hashes and vertices
              (representing a valid hash map)
            + A list of ordered vertices, based on inverse topological ordering
            + A dictionary corresponding to edges from vertex to a list of
              adjacent vertices
            + A dictionary of labels

    big_dag : tuple, default=None
        The dag on which the dag will be appended.
        If None: builds it from dag, else the tuple is the same as the format
        of this function output.

    merge_features : bool, default=True
        If True increments frequencies when a same element is found, else keeps
        them as vectors.

    Returns
    -------
    big_dag : tuple
        A tuple representing Big_Dag:
            + A dictionary from vertices to tuples containing the subtree size,
              frequencies and node ID
            + A dictionary between hashes and vertices
              (representing a valid hash map)
            + A dictionary corresponding to edges from vertex to a list of
              adjacent vertices.
            + A dictionary of labels.

    """
    if big_dag is None:
        D_labels = dict()
        D_hash_map = dict()
        D_vertices = dict()
        D_edges = dict()
        nodes_idx = 0
        nf = 1
    else:
        (D_vertices, D_hash_map, D_edges, D_labels) = big_dag
        nodes_idx = len(D_vertices.keys())
        if not merge_features:
            f = True
            for v in D_vertices.keys():
                D_vertices[v][1].append(0)
                if f:
                    nf = len(D_vertices[v][1])
                    f = False
            if f:
                nf = 1

    (vertices, vertices_hash_map, v_ordered, edges, labels) = dag
    for q in v_ordered:
        key = vertices[q][2]
        if key in D_hash_map:
            node = D_hash_map[key][0]
            # update frequency
            if merge_features:
                D_vertices[node][1] += vertices[q][1]
            else:
                D_vertices[node][1][-1] += vertices[q][1]
        else:
            D_labels[nodes_idx] = labels[q]
            d_edges = list()
            # in order to avoid collisions
            v_nodes = set()
            for c in edges[q]:
                ck = vertices[c][2]
                if ck in D_hash_map:
                    node = D_hash_map[ck][0]
                    # add edges with their new indexes
                    if node not in v_nodes:
                        d_edges.append(node)
                        v_nodes.add(node)

            D_edges[nodes_idx] = d_edges
            D_hash_map[key] = [nodes_idx]

            if merge_features:
                freq = vertices[q][1]
            else:
                freq = (nf-1)*[0]+[vertices[q][1]]

            D_vertices[nodes_idx] = [vertices[q][1], freq, key]
            nodes_idx += 1

    return (D_vertices, D_hash_map, D_edges, D_labels)


if __name__ == "__main__":
    # Dag Example from original paper
    tree_a = ({0, 1, 2, 3},
              {0: [1, 2], 1: [3], 2: [], 3: []},
              {0: 'a', 1: 'b', 2: 'd', 3: 'c'})
    tree_b = ({0, 1, 2, 3, 4}, {0: [1, 2], 1: [3], 2: [4], 3: [], 4: []},
              {0: 'a', 1: 'b', 2: 'c', 3: 'c', 4: 'd'})

    todd_a = odd(tree_a[0], tree_a[1], tree_a[2])
    todd_b = odd(tree_b[0], tree_b[1], tree_b[2])

    Atree = tuple(hash_trees(todd_a))+tuple([])+(todd_a[1], todd_a[3])
    Btree = tuple(hash_trees(todd_b))+tuple([])+(todd_b[1], todd_b[3])
    big_dag = None
    big_dag = big_dag_append(Atree, big_dag)
    big_dag = big_dag_append(Btree, big_dag)

    print("Tree A:\n", Atree, "\nTree B:\n", Btree, "\nBig Dag:\n", big_dag)
