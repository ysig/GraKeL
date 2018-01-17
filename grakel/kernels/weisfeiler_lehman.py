"""The weisfeiler lehman kernel :cite:`Shervashidze2011WeisfeilerLehmanGK`."""

import itertools

import numpy as np

from grakel.graph import graph


def weisfeiler_lehman(X, Y, Lx, Ly, base_kernel, niter=5):
    """Compute the Weisfeiler Lehman Kernel.

     See :cite:`Shervashidze2011WeisfeilerLehmanGK`.

    Parameters
    ----------
    X,Y : *valid-graph-format*
        The pair of graphs on which the kernel is applied.

    L{x,y} : dict
        Corresponding graph labels for vertices.

    base_kernel : function (graph, graph -> number)
        A valid pairwise base_kernel for graphs.

    niter : int, default=5
        The number of iterations.

    Returns
    -------
    kernel : number
        The kernel value.

    """
    Ga = graph(X, Lx)
    Gb = graph(Y, Ly)
    bkm = lambda x, y: np.array([[base_kernel(x[0], y[0])]])
    return weisfeiler_lehman_matrix({0: Ga}, bkm, {0: Gb}, niter)[0, 0]


def weisfeiler_lehman_matrix(
        Graphs_a, base_kernel_matrix, Graphs_b=None,  niter=5):
    """Compute the Weisfeler Lehman Kernel Matrix.

    See :cite:`Shervashidze2011WeisfeilerLehmanGK`.

    Parameters
    ----------
    Graphs_{x,y} : dict, default_y=None
        Enumerative dictionary of graph type objects with keys from 0 to the
        number of values. If value of Graphs_y is None the kernel matrix is
        computed between all pairs of Graphs_x where in another case the
        kernel_matrix rows correspond to elements of Graphs_y, and columns
        to the elements of Graphs_x.

    base_kernel_matrix : function (dict(graph), dict(graph) -> np.array)
        Rows of the np.array should correspond to the second dictionary of
        graphs and cols to the first. If the second graph is None the kernel
        should be computed upon itself.

    niter : int, default=5
        The number of iterations.

    Returns
    -------
    kernel_mat : np.array
        The kernel matrix.

    """
    if niter < 1:
        raise ValueError('n_iter must be an integer bigger than zero')

    ng_a = len(Graphs_a.keys())

    if Graphs_b is None:
        g_iter = [Graphs_a[i] for i in range(0, ng_a)]
        ng = ng_a
        ng_b = ng_a
    else:
        ng_b = len(Graphs_b.keys())
        g_iter = itertools.chain([Graphs_a[i] for i in range(0, ng_a)],
                                 [Graphs_b[i] for i in range(0, ng_b)])
        ng = ng_a + ng_b

    Gs = dict()
    G_ed = dict()
    L_orig = dict()

    # get all the distinct values of current labels
    WL_labels = dict()
    WL_labels_inverse = dict()
    distinct_values = set()

    for (i, g) in enumerate(g_iter):
        g.desired_format("dictionary")
        Gs[i] = g
        G_ed[i] = g.edge_dictionary
        L_orig[i] = g.node_labels
        # calculate all the distinct values
        distinct_values |= set(L_orig[i].values())
    distinct_values = sorted(list(distinct_values))

    # assign a number to each label
    label_count = 0
    for dv in distinct_values:
        WL_labels[label_count] = dv
        WL_labels_inverse[dv] = label_count
        label_count += 1

    for i in range(ng):
        L = dict()
        for k in L_orig[i].keys():
            L[k] = WL_labels_inverse[L_orig[i][k]]

        # add new labels
        Gs[i].relabel(L)

    kernel = np.zeros(shape=(ng_b, ng_a))
    kernel = np.add(kernel, base_kernel_matrix(Graphs_a, Graphs_b))

    for i in range(niter):
        label_set = set()
        L_temp = dict()
        for j in range(ng):
            # Find unique labels and sort
            # them for both graphs
            # Keep for each node the temporary
            L_temp[j] = dict()
            for v in G_ed[j].keys():
                nlist = list()
                for neighbor in G_ed[j][v].keys():
                    nlist.append(Gs[j].node_labels[neighbor])
                credential = str(Gs[j].node_labels[v])+","+str(sorted(nlist))
                L_temp[j][v] = credential
                label_set.add(credential)

        label_list = sorted(list(label_set))
        for dv in label_list:
            WL_labels[label_count] = dv
            WL_labels_inverse[dv] = label_count
            label_count += 1

        # Recalculate labels
        for j in range(ng):
            L = dict()
            for k in L_temp[j].keys():
                L[k] = WL_labels_inverse[L_temp[j][k]]
            # relabel
            Gs[j].relabel(L)

        # calculate kernel
        kernel = np.add(kernel, base_kernel_matrix(Graphs_a, Graphs_b))

    # Restore original labels
    for i in range(ng):
        Gs[i].relabel(L_orig[i])

    if Graphs_b is None:
        kernel = np.triu(kernel) + np.triu(kernel, 1).T

    return kernel
