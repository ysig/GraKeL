"""The sugraph mathing kernel as defined by :cite:`Kriege2012SubgraphMK`."""

from grakel.graph import graph

global kv_default, ke_default, lw_default

# Define default vertex, edge and lambda weight functions
kv_default = lambda x, y, Lx, Ly: 1 if Lx[x] == Ly[y] else 0


def ke_default(x, y, Ex, Ey, Lex, Ley):
    """Calculate the default edge kernel :eq:`sm_ke_def`."""
    cond_a, cond_b = x in Ex, x in Ey
    if (cond_a and cond_b and Lex[x] == Ley[y]) or \
            ((not cond_a) and (not cond_b)):
        return 1
    else:
        return 0


lw_default = lambda x: int(bool(x))


def subgraph_matching(
        X, Y, Lx, Ly, Lex, Ley, kv=kv_default, ke=ke_default, lw=lw_default):
    """Calculate the subgraph matching kernel.

    See :cite:`Kriege2012SubgraphMK`.

    Parameters
    ----------
    X,Y : *valid-graph-format*
        The pair of graphs on which the kernel is applied.

    L{x,y} : dict
        Corresponding graph labels for vertices.

    kv : function (symbol, symbol, dict, dict -> number),
    default=:math:`k_{v}^{default}`
        The kernel function for two nodes based on label dictionaries.

    ke : function (tuple, tuple, dict, dict, dict, dict -> number),
    default=:math:`k_{e}^{default}`
        The kernel function for edges (tuples) of two graphs based on their
        edges (dict 1,2) and edge labels (dict 3,4).

    lw : function (set -> number), default=:math:`l_{w}^{default}`
        A lambda weight function for cliques.

    Returns
    -------
    kernel : number
        The kernel value.

    """
    Gx = graph(X, Lx, Ley)
    Gy = graph(Y, Ly, Ley)
    return subgraph_matching_pair(Gx, Gy, kv=kv, ke=ke, lw=lw)


def subgraph_matching_pair(
        Gx, Gy, kv=kv_default, ke=ke_default, lw=lw_default):
    """Calculate the subgraph matching kernel.

    See :cite:`Kriege2012SubgraphMK`.

    Parameters
    ----------
    G{x,y} : graph
        The pair graphs on which the kernel is applied.

    kv : function (symbol, symbol, dict, dict -> number),
    default=:math:`k_{v}^{default}`
        The kernel function for two nodes based on label dictionaries.

    ke : function (tuple, tuple, dict, dict, dict, dict -> number),
    default=:math:`k_{e}^{default}`
        The kernel function for edges (tuples) of two graphs based on their
        edges (dict 1,2) and edge labels (dict 3,4).

    lw : function (set -> number), default=:math:`l_{w}^{default}`
        A lambda weight function for cliques

    Returns
    -------
    kernel : number
        The kernel value.

    """
    # Calculate product graph
    Vp, Ep, c = weighted_product_graph(Gx, Gy, kv=kv, ke=ke)

    # Initialize values
    C = set()
    value = [.0]

    subgraph_matching_core(1, C, Vp, Ep, c, lw, value)
    return value[0]


def subgraph_matching_core(w, C, P, Ep, c, lw, value):
    """Calculate the core algorithm of the subgraph matching kernel.

    See :cite:`Kriege2012SubgraphMK` (p.5, Algorithm 1).

    Parameters
    ----------
    G{x,y} : graph
        The pair graphs on which the kernel is applied.

    w : number
        The current weight.

    C : set
        The current clique.

    P : set
        Current vertices.

    c : dict
        Costs of vertices, edges.

    lw : function: (set -> number)
        A lambda weight function for cliques.

    value : list
        One element list used as mutable global variable for
        the recursion of the algorithm.

    Returns
    -------
    None.

    """
    while len(P) > 0:
        v = P.pop()
        w = w*c[v]

        for u in list(C):
            w = w*c[(u, v)]
        C.add(v)
        value[0] += w*lw(C)

        # Create the new P
        Pp = P.copy()
        Pp.add(v)
        Pp &= Ep[v]

        # apply subgraph matching for the new clique
        subgraph_matching_core(w, C, Pp, Ep, c, lw, value)


def weighted_product_graph(Gx, Gy, kv, ke):
    """Calculate the weighted product graph.

    For a definition see :cite:`Kriege2012SubgraphMK` (p.5, Definition 5).

    Parameters
    ----------
    G{x,y} : graph
        The pair graphs on which the kernel is applied.

    kv : function (symbol, symbol, dict, dict -> number)
        The kernel function for two nodes based on label dictionaries.

    ke : function (tuple, tuple, dict, dict, dict, dict -> number)
        The kernel function for edges (tuples) of two graphs based
        on their edges (dict 1,2) and edge labels (dict 3,4).

    Returns
    -------
    Vp : set
        The vertices of the weighted product graph.

    Ep : dict
        The edges of the weighted product graph.

    c : dict
        The cost dictionary for vertices and edges.

    """
    Gx.desired_format("all")
    Gy.desired_format("all")

    Lx = Gx.get_labels()
    Ly = Gy.get_labels()
    Lex = Gx.get_labels("edge")
    Ley = Gy.get_labels("edge")
    Ex = Gx.edge_dictionary
    Ey = Gy.edge_dictionary

    Kv = lambda x, y: kv(x, y, Lx, Ly)
    Ke = lambda x, y: ke(x, y, Ex, Ey, Lex, Ley)

    # initialise cost function
    c = dict()
    # Calculate valid vertices
    Vp = set()
    Ep = dict()

    # calculate product graph vertex set
    for i in range(0, Gx.n):
        for j in range(0, Gy.n):
            if(i != j):
                value = Kv(i, j)
                if(value > 0):
                    # add to vertex set
                    Vp.add((i, j))
                    # add cost
                    c[(i, j)] = value
                    # initialise an empty set for neighbors
                    if (i, j) not in Ep:
                        Ep[(i, j)] = set()
                else:
                    c[(i, j)] = 0
            else:
                c[(i, j)] = 0

    # calculate product graph valid edges
    for v in list(Vp):
        for w in list(Vp):
            if(v[0] != w[0] and v[1] != w[1]):
                value = Ke(v, w)
                if(value > 0):
                    # add edge
                    Ep[v].add(w)
                    # store value
                    c[(v, w)] = value
                else:
                    c[(v, w)] = 0
            else:
                c[(v, w)] = 0
    return Vp, Ep, c
