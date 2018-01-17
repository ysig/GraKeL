"""The subtree kernel as defined in :cite:`Ramon2003ExpressivityVE`."""

import itertools

from grakel.graph import graph
from grakel.tools import nested_dict_add
from grakel.tools import nested_dict_get


def subtree_rg(X, Y, Lx, Ly, h=5):
    """Calculate Ramon Gartner subtree kernel.

    See :cite:`Ramon2003ExpressivityVE`, s. 5.

    Parameters
    ----------
    X,Y : *valid-graph-format*
        The pair of graphs on which the kernel is applied.

    L{x,y} : dict
        Corresponding graph labels for vertices.

    h : int, default=5
        The sub-tree depth.

    Returns
    -------
    kernel : number
        The kernel value.

    """
    Gx = graph(X, Lx)
    Gy = graph(Y, Ly)
    return subtree_rg_pair(Gx, Gy, h=5)


def subtree_rg_pair(Gx, Gy, h=5):
    """Calculate The Ramon Gartner subtree kernel.

    See :cite:`Ramon2003ExpressivityVE`, s. 5.

    Parameters
    ----------
    G_{x,y} : graph
        The pair graphs on which the kernel is applied.

    h : int, default=5
        The sub-tree depth.

    Returns
    -------
    kernel : number
        The kernel value.

    """
    Gx.desired_format("dictionary")
    Gy.desired_format("dictionary")
    kernel = 0
    dynamic_dictionary = dict()
    for u in Gx.edge_dictionary.keys():
        for v in Gy.edge_dictionary.keys():
            kernel += \
                subtree_rg_core_dynamic(u, v, Gx, Gy, h, dynamic_dictionary)
    return kernel


def subtree_rg_core_dynamic(
        u, v, g_x, g_y, h, dynamic_dict, p_u=None, p_v=None):
    """Calculate the inner kernel of the Ramon Gartner subtree kernels sum.

    See :cite:`Ramon2003ExpressivityVE`, s. 5.

    Parameters
    ----------
    u,v : hashable(s)
        Vertices that correspond to graphs g_x, g_y.

    g_{x,y} : graph
        Graph formats for the corresponding graphs.

    h : int
        The height of the subtree exploration.

    dynamic_dict : dict
        A dictionary that is used for applying dynamic programming.
    p_{u,v} : hashable(s) or None, default=None
        Parent nodes, to avoid repeating explorations.

    Returns
    -------
    kernel : number
        The kernel value.

    Note
    ----
    This is an efficiency proposal using various data structures and dynamic
    programming assuming that for each node you take all of its neighbors and
    not its predecesor. A dictionary is being formed that holds level, u, v and
    pred of u and pred of v.

    """
    if h == 1:
        return int(g_x.label(u) == g_y.label(v))

    elif h > 1:
        # Calculate R-maximal: First group all nodes with the same label make
        # a list of lists of all pairs for each label.
        lbx, lby = g_x.get_label_group(), g_y.get_label_group()

        # Calculate neighbors and remove previous.
        # Avoid backtracking to father.
        nx = g_x.neighbors(u)
        if p_u is not None and p_u in nx:
            nx.remove(p_u)

        # Do the same for y.
        ny = g_y.neighbors(v)
        if p_v is not None and p_v in ny:
            ny.remove(p_v)

        # What happens when one list has no neighbors?
        xy_and, xy_or = len(nx)*len(ny), len(nx)+len(ny)
        if(xy_or == 0):
            # both trees finish at the same point
            return int(g_x.label(u) == g_y.label(v))
        elif(xy_and == 0):
            # else trees are different so output zero (?)
            return 0

        # Calculate the set of neighbours with common labels.
        Rset = []
        snx = set(nx)
        sny = set(ny)
        for kx in lbx:
            if kx in lby:
                # Substract from the common label only the valid neighbors.
                snxk = list(set(lbx[kx]).intersection(snx))
                snyk = list(set(lby[kx]).intersection(sny))
                # If both sets are bigger than zero calculate all the possible
                # pairs and add to the Rset
                if len(snxk) > 0 and len(snyk) > 0:
                    pair = [lbx[kx], lby[kx]]
                    Rset += list(itertools.product(*pair))

        # Designate all nodes with the same start and all with the same finish.
        right, left = dict(), dict()

        # Dictionary to store for every pair the index.
        Rset_dict, Rset_inv_dict = dict(), dict()

        # Enumerate all index pairs.
        for (w, z) in Rset:
            if w not in right:
                right[w] = list()
            if z not in left:
                left[z] = list()
            right[w].append(len(Rset_dict))
            left[z].append(len(Rset_dict))
            Rset_inv_dict[len(Rset_dict)] = (w, z)
            Rset_dict[(w, z)] = len(Rset_dict)

        # For each tuple of indexes store to left and to right bins according
        # with the fact that they have the same left or right indexes.
        Kbins_flat_r, Kbins_flat_l = list(), list()
        for k in right.keys():
            Kbins_flat_r.append(list(right[k]))
        for k in left.keys():
            Kbins_flat_l.append(list(left[k]))

        # Calculate all possible combinations between all the different bins.
        # Product equals combinations because the lists are disjoint.
        setA, setB = set(itertools.product(*Kbins_flat_r)),\
            set(itertools.product(*Kbins_flat_l))

        # Produce R-maximal by taking the intersection between
        # all valid subsets.
        Rmaximal = setA.intersection(setB)

        # Calculate for all the valid tuples, their kernel values by using
        # dynamic programming on level, nodes, predecessors.
        Rv = dict()
        for (w, wp) in {Rset_inv_dict[i] for s in Rmaximal for i in s}:
            r = nested_dict_get(dynamic_dict, h-1, u, v, w, wp)
            if (r is None):
                kh = subtree_rg_core_dynamic(
                    w, wp, g_x, g_y, h-1, dynamic_dict, u, v)
                nested_dict_add(dynamic_dict, kh, h-1, u, v, w, wp)
                Rv[Rset_dict[(w, wp)]] = kh
            else:
                Rv[Rset_dict[(w, wp)]] = r

        # Holds the values of all sets inside R, that are not zero
        M_values = dict()
        for s in Rmaximal:
            # Keep only non zero elements
            # (all subsets containing zero-elements will be zero)
            non_zero_elements = set()
            for i in s:
                if(Rv[i] != 0):
                    non_zero_elements.add(i)

            if bool(non_zero_elements):
                # Calculate all the subset values faster
                # using dynamic programming
                plough_subsets(non_zero_elements, Rv, M_values)

        return sum(M_values.values())

    else:
        raise ValueError('h must pe a positive integer')


def plough_subsets(initial_set, Rv, value):
    """Calculate all the subset kernel values, without repeating operations.

    Parameters
    ----------
    initial_set : set
        The set that will be ploughed further.

    Rv : dict
        Values for single elements (dictionary on integers).

    value : dict
        Stores all values for the relevant sets.

    Returns
    -------
    None.

    """
    # If only one element return its value
    if (len(initial_set) == 1):
        frozen_initial_set = frozenset(initial_set)
        p = initial_set.pop()
        value[frozen_initial_set] = Rv[p]

    elif (len(initial_set) > 1):
        # Explore all subsets for value dictionary
        # to fill and calculate current.
        flag = True
        frozen_initial_set = frozenset(initial_set)
        for s in initial_set:
            # copy set, delete element and if not calculated: plough
            temp_set = initial_set.copy()
            temp_set.discard(s)
            frozen_temp_set = frozenset(temp_set)
            if frozen_temp_set not in value:
                plough_subsets(temp_set, Rv, value)

            if(flag):
                value[frozen_initial_set] = Rv[s]*value[frozen_temp_set]
                flag = False
