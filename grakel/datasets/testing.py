"""Dataset generation fo rtesting purposes."""
import numpy as np
from numbers import Real
from sklearn.utils import check_random_state


def generate_dataset(
        n_graphs=100,
        r_vertices=(1, 20),
        r_connectivity=(0.2, 0.8),
        r_weight_edges=(1, 1),
        n_graphs_test=20,
        random_state=None,
        features="nl"):
    """
    Generate a random graph datset for testing.

    Labels there is always one distinct label for a graph inside test.

    Parameters
    ----------
    n_graphs : int
        The total dataset size.

    r_vertices : tuple, len=2
        A range of vertices that can appear inside the random generated graphs.

    r_connectivity : tuple, len=2
        A tuple of Real numbers corresponding to the amount of connectivity inside the extracted graph.

    r_weight_edges : tuple, len=2
        A tuple of Real numbers corresponding to the range of weight values between the graphs edges.

    n_graphs_test : int
        Number of graphs inside test.

    random_state :  RandomState or int, default=None
        A random number generator instance or an int to initialize a RandomState as a seed.

    features : str or tuple or None
        - str : can be either "nl" or "na", "el" or "ea"

        - tuple : cartesian product of ["nl", "na"] and ["el", "ea"].
            If "na" or "ea" appears it can be followed by an integer bigger
            than zero concerning the attribute size, default is 4, whereas if
            "nl" or "el" appears it can be followed by an integer bigger
            than zero concerning the number of distinct labels.

    Returns
    -------
    train : list
        A list of iterables, where its contains an adjacency matrix, a dict of node-labels/attributes
        (optional) and a dict of edge-labels/attributes (optional). Corresponds to the fit/fit_transform
        input.

    test : list
        A list of iterables, where its contains an adjacency matrix, a dict of node-labels/attributes
        (optional) and a dict of edge-labels/attributes (optional). Corresponds to the transform
        input, where there exists a new label from fit, if node-labels is desired by the user and
        similar with transform.

    """
    if type(n_graphs) is not int or n_graphs < 1:
        raise TypeError('Number of graphs must be an integer bigger than 1.')

    if type(n_graphs_test) is not int or n_graphs_test < 1 or n_graphs_test >= n_graphs:
        raise TypeError('Number of graphs inside test must be smaller than the total '
                        'number of graphs and positive.')

    rs = check_random_state(random_state)

    if (type(r_vertices) is not tuple or len(r_vertices) != 2
            or r_vertices[0] > r_vertices[1] or r_vertices[0] < 0 or
            type(r_vertices[0]) is not int or type(r_vertices[1]) is not int):
        raise TypeError('r_vertices must be a tuple containing the minimum '
                        'and the maximum numebr of vertices graph ')

    if (type(r_connectivity) is not tuple or len(r_connectivity) != 2
            or r_connectivity[0] > r_connectivity[1] or r_connectivity[0] < 0 or r_connectivity[1] > 1
            or not isinstance(r_connectivity[0], Real) or not isinstance(r_connectivity[1], Real)):
        raise TypeError('r_connectivity must be a Real tuple containing the minimum '
                        'and the maximum percent of allowed connectivity in the range of [0, 1]')

    if (type(r_weight_edges) is not tuple or len(r_weight_edges) != 2
            or r_weight_edges[0] > r_weight_edges[1] or r_weight_edges[0] <= 0
            or not isinstance(r_weight_edges[0], Real) or not isinstance(r_weight_edges[1], Real)):
        raise TypeError('r_weight_edges must be a tuple containing the minimum '
                        'and the maximum values for an edge weight ')

    if type(features) is str or features is None:
        features = (features, )
    pnl, pel, pna, pea = False, False, False, False
    dna, dea, nnl, nel = 4, 4, 5, 5

    if type(features) is not tuple or len(features) > 4 or len(features) == 0:
        raise TypeError('features can either be a tuple of at most four elements and at least one'
                        'or a string')
    elif features[0] is not None:
        fiter = iter(features)
        flag = True
        try:
            f = next(fiter)
        except StopIteration:
            flag = False
        while flag:
            if type(f) is str:
                if f == "nl":
                    if pna:
                        raise ValueError('The dataset can either have node-labels or node-attributes')
                    pnl = True
                    try:
                        f = next(fiter)
                        if type(f) is int:
                            if f >= 2:
                                nel = f
                                f = next(fiter)
                            else:
                                raise TypeError('The number of distinct node labels must '
                                                'bigger equal to 2')
                    except StopIteration:
                        break
                elif f == "na":
                    try:
                        if pnl:
                            raise ValueError('The dataset can either have node-labels '
                                             'or node-attributes')
                        pna = True
                        f = next(fiter)
                        if type(f) is int:
                            if f >= 1:
                                dna = f
                                f = next(fiter)
                            else:
                                raise TypeError('The dimension of node-attributes must be positive')
                    except StopIteration:
                        break
                elif f == "el":
                    if pea:
                        raise ValueError('The dataset can either have edge-labels or edge-attributes')
                    pel = True
                    try:
                        f = next(fiter)
                        if type(f) is int:
                            if f >= 2:
                                nel = f
                                f = next(fiter)
                            else:
                                raise TypeError('The number of distinct edge labels must bigger '
                                                'equal to 2')
                    except StopIteration:
                        break
                elif f == "ea":
                    if pel:
                        raise ValueError('The dataset can either have edge-labels or edge-attributes')
                    pea = True
                    try:
                        f = next(fiter)
                        if type(f) is int:
                            if f >= 1:
                                dea = f
                                f = next(fiter)
                            else:
                                raise TypeError('The dimension of node-attributes must be positive')
                    except StopIteration:
                        break
                else:
                    raise ValueError('feature markers can either be "nl", "na", "el", "ea"')
            else:
                raise TypeError('feature markers must be a string from either "nl", "na", "el", "ea"')

    # Lets produce some graphs
    def randint(*args, **kargs):
        if kargs["low"] == kargs["high"]:
            return np.full(kargs["size"], kargs["low"])
        else:
            return rs.randint(low=kargs["low"], high=kargs["high"], size=kargs["size"])

    def rand(*args, **kargs):
        low, high = kargs["low"], kargs["high"]
        if low == high:
            return np.full(kargs["size"], low)
        else:
            return rs.rand(*kargs["size"])*(high - low) + low

    graphs = list()
    for i in range(n_graphs):
        # Pick a size
        size = randint(low=r_vertices[0], high=r_vertices[1], size=(1,))[0]
        # Draw all possible tuples
        max_ne = int((size-1)*size/2)
        nedges = randint(low=int(max_ne * r_connectivity[0]), high=int(max_ne * r_connectivity[1]),
                         size=(1,))
        if nedges == max_ne:
            A = rand(low=r_weight_edges[0], high=r_weight_edges[1], size=(size, size)).astype(float)
            upper_triangular = np.triu(A, 1)
            A = upper_triangular + upper_triangular.T
        elif nedges == 0:
            A = np.zeros(shape=(size, size), dtype=float)
        else:
            if nedges > max_ne/2:
                A = rand(low=r_weight_edges[0], high=r_weight_edges[1], size=(size, size)).astype(float)
                grey_idxs = set(rs.choice(max_ne, max_ne - nedges, replace=False))
                c = 0
                for i in range(size):
                    A[i, i] = 0
                    for j in range(i+1, size):
                        if c in grey_idxs:
                            A[i, j] = 0
                        A[j, i] = A[i, j]
                        c += 1
            else:
                A = np.zeros(shape=(size, size), dtype=float)
                grey_idxs = set(rs.choice(max_ne, max_ne - nedges, replace=False))
                c = 0
                for i in range(size):
                    A[i, i] = 0
                    for j in range(i+1, size):
                        if c in grey_idxs:
                            A[i, j] = rand(low=r_weight_edges[0],
                                           high=r_weight_edges[1], size=(1,))[0].astype(float)

                        A[j, i] = A[i, j]
                        c += 1
        graphs.append(A)

    graph_node_labels = list()
    if pnl:
        # Node Labels
        for i in range(n_graphs - n_graphs_test):
            node_labels = dict(enumerate(rs.choice(nnl-1, graphs[i].shape[0])))
            graph_node_labels.append(node_labels)
        for i in range(n_graphs - n_graphs_test, n_graphs):
            node_labels = dict(enumerate(rs.choice(nnl, graphs[i].shape[0])))
            graph_node_labels.append(node_labels)
        # existance gurantee for the new label inside the graphs of test
        graph_node_labels[-1][0] = nnl
    elif pna:
        # Node Attributes
        for i in range(n_graphs):
            node_labels = dict(enumerate(rs.rand(graphs[i].shape[0], dna)))
            graph_node_labels.append(node_labels)

    graph_edge_labels = list()
    if pel:
        # Edge Labels
        for i in range(n_graphs - n_graphs_test):
            g = graphs[i]
            idx_i, idx_j = np.where(g > 0)
            edge_labels = dict(zip(zip(idx_i, idx_j), rs.choice(nel-1, idx_i.shape[0])))
            graph_edge_labels.append(edge_labels)
        has_not_the_new_label = True
        for i in range(n_graphs - n_graphs_test, n_graphs):
            g = graphs[i]
            idx_i, idx_j = np.where(g > 0)
            edge_labels = dict(zip(zip(idx_i, idx_j), rs.choice(nel, idx_i.shape[0])))
            if has_not_the_new_label and idx_i.shape[0] > 0:
                # existance gurantee for the new label inside the graphs of test
                # if not the new set of new graphs has no edges
                has_not_the_new_label = False
                edge_labels[idx_i[0], idx_j[0]] = nnl
            graph_edge_labels.append(edge_labels)
    elif pea:
        # Edge Attributes
        for i in range(n_graphs):
            g = graphs[i]
            idx_i, idx_j = np.where(g > 0)
            edge_labels = dict(zip(zip(idx_i, idx_j), rs.rand(idx_i.shape[0], dea)))
            graph_edge_labels.append(edge_labels)

    # Zip all together base
    out = (graphs,)
    if pna or pnl:
        out += (graph_node_labels,)
        if pel or pea:
            out += (graph_edge_labels,)
    elif pel or pea:
        graph_node_labels = iter(None for i in range(len(graphs)))
        out += (graph_node_labels, graph_edge_labels)

    out = list(zip(*out))
    return out[:-n_graphs_test], out[-n_graphs_test:]
