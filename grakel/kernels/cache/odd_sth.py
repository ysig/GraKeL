"""The ODD-Sth kernel as defined in :cite:`Martino2012ATK`."""

import itertools

import numpy as np

from grakel.graph import graph


def odd_sth(X, Y, Lx, Ly, h=None):
    """ODD-Sth kernel as proposed in :cite:`Martino2012ATK`.

    Parameters
    ----------
    X,Y : *valid-graph-format*
        The pair of graphs on which the kernel is applied.

    L{x,y} : dict
        Coresponding graph labels for vertices.

    h : int, default=None
        Maximum (single) dag height.
        If None there is no restriction.

    Returns
    -------
    kernel : number
        The kernel value.

    """
    Gx = graph(X, Lx)
    Gy = graph(Y, Ly)
    return float(odd_sth_matrix({0: Gx}, {0: Gy}, h=None)[0, 0])


def odd_sth_matrix(Graphs_x, Graphs_y=None, h=None):
    """ODD-Sth kernel, as proposed in :cite:`Martino2012ATK`.

    Parameters
    ----------
    Graphs_{x,y} : dict, default_y=None
        Enumerative dictionary of graph type objects with keys from 0
        to the number of values. If value of Graphs_y is None the kernel
        matrix is computed between all pairs of Graphs_x where in another
        case the kernel_matrix rows correspond to elements of Graphs_y
        and columns to the elements of Graphs_x.

    h : int, default=None
        Maximum (single) dag height.
        If None there is no restriction.

    Returns
    -------
    kernel_matrix : np.array
        The kernel matrix.

    """
    h_x = len(Graphs_x.keys())

    if Graphs_y is None:
        h_y = h_x
        g_iter = Graphs_x.values()
        ng = h_x
        pairs = [(i, j) for i in range(0, h_x) for j in range(i, h_x)]
        offset = 0
    else:
        h_y = len(Graphs_y.keys())
        g_iter = itertools.chain(Graphs_x.values(), Graphs_y.values())
        ng = h_x+h_y
        pairs = list(itertools.product(range(h_x, ng), range(0, h_x)))
        offset = h_x

    big2Dag = None
    for g in g_iter:
        big_dag = make_big_dag(g, h=h)
        big2Dag = big_dag_append(big_dag, big2Dag, merge_features=False)

    C = dict()
    for v in big2Dag[0].keys():
        # number of identical subtrees
        # equal the D element
        C[v] = big2Dag[0][v][0]
    K = np.zeros(shape=(h_y, h_x))
    for (i, j) in pairs:
        k = 0
        for v in C.keys():
            # TODO(write K calculation in terms of matrix multiplication)
            k += big2Dag[0][v][1][i]*big2Dag[0][v][1][j]*C[v]
        K[i-offset, j] = k

    if Graphs_y is None:
        K = np.triu(K) + np.triu(K, 1).T

    return K


def make_big_dag(g, h=None):
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
        dag_odd = make_dag_odd(v, g, h=h)
        dag = tuple(hash_trees(dag_odd))+tuple([])+(dag_odd[1], dag_odd[3])
        big_dag = big_dag_append(dag, big_dag)

    _, D_edges, D_ordering, _ = odd(big_dag[0], big_dag[2], big_dag[3])

    big_dag = (big_dag[0],
               big_dag[1],
               sorted(list(big_dag[0].keys()),
               key=lambda x:
                      (D_ordering[x], big_dag[3][x])), D_edges, big_dag[3])

    return big_dag


def make_dag_odd(v, g, h=None):
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
    vertices, edges = dag(v, g, h=None)
    return odd(vertices, edges, g.get_labels(purpose='any'))


def dag(v, g, h=None):
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
    if h is None:
        h = -1
    else:
        h = int(h)
        if h <= 0:
            raise ValueError('depth of visits must be bigger than zero')

    q = [(v, 0)]
    vertices = dict()
    edges = dict()

    vertices[v] = 0
    while len(q) > 0:
        u, level = q.pop(0)

        if level == h:
            break

        if u not in edges:
            edges[u] = list()
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
        raise ValueError('unsupported vertices type')

    for (k, e) in edges.items():
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

    See :cite:`Aiolli2006FastOK` and notated in :cite:`Martino2012ATK`.

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


"""
# Dag Example from original paper
tree_a = \
    ({0,1,2,3}, {0:[1,2], 1:[3], 2:[], 3:[]}, {0:'a', 1:'b', 2:'d', 3:'c'})
tree_b = \
    ({0,1,2,3,4}, {0:[1,2], 1:[3], 2:[4], 3:[], 4:[]},
        {0:'a', 1:'b', 2:'c', 3:'c', 4:'d'})

todd_a = odd(tree_a[0], tree_a[1], tree_a[2])
todd_b = odd(tree_b[0], tree_b[1], tree_b[2])

Atree = tuple(hash_trees(todd_a))+tuple([])+(todd_a[1],todd_a[3])
Btree = tuple(hash_trees(todd_b))+tuple([])+(todd_b[1],todd_b[3])
big_dag = None
big_dag = big_dag_append(Atree, big_dag)
big_dag = big_dag_append(Btree, big_dag)

print("Tree A:\n",Atree,"\nTree B:\n",Btree,"\nBig Dag:\n", big_dag)
"""
