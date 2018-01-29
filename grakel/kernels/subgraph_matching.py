"""The sugraph mathing kernel as defined by :cite:`Kriege2012SubgraphMK`."""
import collections
import warnings

from sklearn.utils import Bunch

from grakel.kernels import kernel
from grakel.graph import graph

global kv_default, ke_default, lw_default

# Define default vertex, edge and lambda weight functions
kv_default = lambda x, y, Lx, Ly: 1 if Lx[x] == Ly[y] else 0


def ke_default(x, y, Ex, Ey, Lex, Ley):
    """Calculate the default edge kernel :eq:`sm_ke_def`."""
    cond_a, cond_b = x in Ex, y in Ey
    if (cond_a and cond_b and Lex[x] == Ley[y]) or \
            ((not cond_a) and (not cond_b)):
        return 1
    else:
        return 0


lw_default = lambda x: int(bool(x))


class subgraph_matching(kernel):
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

    Attributes
    ----------
    _kv : function (symbol, symbol, dict, dict -> number)
        The kernel function for two nodes based on label dictionaries.

    _ke : function (tuple, tuple, dict, dict, dict, dict -> number)
        The kernel function for edges (tuples) of two graphs based on their
        edges (dict 1,2) and edge labels (dict 3,4).

    _lw : function (set -> number)
        A lambda weight function for cliques.

    """

    _graph_format = "all"

    def __init__(self, **kargs):
        """Initialise a subtree_wl kernel."""
        self._valid_parameters |= {"kv", "ke", "lw"}

        super(subgraph_matching, self).__init__(**kargs)
        self._kv = kargs.get("kv", kv_default)
        self._ke = kargs.get("ke", ke_default)
        self._lw = kargs.get("lw", lw_default)

    def pairwise_operation(self, x, y):
        """Calculate the subgraph matching kernel.

        See :cite:`Kriege2012SubgraphMK`.

        Parameters
        ----------
        x, y : tuples
            *Vertex-set*, *edge-dictionary*, *node-label-dictionary*,
            *edge-labels-dictionary* tuple.

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
        Vp, Ep, c = self._weighted_product_graph(x, y)

        # Initialize values
        b = Bunch(value=.0)

        self._subgraph_matching_core(1, set(), Vp, Ep, c, b)

        return b.value

    def parse_input(self, X):
        """Parse and create features for graphlet_sampling kernel.

        Parameters
        ----------
        X : object
            For the input to pass the test, we must have:
            Each element must be an iterable with at most three features and at
            least one. The first that is obligatory is a valid graph structure
            (adjacency matrix or edge_dictionary) while the second is
            node_labels and the third edge_labels (that fitting the given graph
            format). If None the kernel matrix is calculated upon fit data.
            The test samples.

        Returns
        -------
        out : list
            The extracted adjacency matrices for any given input.

        """
        if not isinstance(X, collections.Iterable):
            raise ValueError('input must be an iterable\n')
        else:
            i = 0
            out = list()
            for x in iter(X):
                if len(x) == 0:
                    warnings.warn('Ignoring empty element on index: '+str(i))
                elif len(x) == 1 and type(x) is graph:
                    g = graph(x.get_adjacency_matrix(),
                              x.get_labels(purpose="adjacency"),
                              x.get_labels(purpose="adjacency",
                                           label_type="edge"),
                              self._graph_format)
                    n = g.nv()
                    E = g.get_edge_dictionary()
                    L = g.get_labels(purpose="dictionary")
                    Le = g.get_labels(purpose="dictionary",
                                      label_type="edge")
                    Er = set((a, b) for a in E.keys()
                             for b in E[a].keys() if a != b)

                    i += 1
                    out.append((n, Er, L, Le))
                elif len(x) == 3:
                    g = graph(x[0], x[1], x[2], "adjacency")
                    g.desired_format(self._graph_format, False)
                    n = g.nv()
                    E = g.get_edge_dictionary()
                    L = g.get_labels(purpose="dictionary")
                    Le = g.get_labels(purpose="dictionary",
                                      label_type="edge")
                    Er = set((a, b) for a in E.keys()
                             for b in E[a].keys() if a != b)
                    i += 1
                    out.append((n, Er, L, Le))
                else:
                    raise ValueError('each element of X must have at least' +
                                     ' one and at most 3 elements\n')

            if i == 0:
                raise ValueError('parsed input is empty')
            print("Input Parsed")
            return out

    def _weighted_product_graph(self, x, y):
        """Calculate the weighted product graph.

        For a definition see :cite:`Kriege2012SubgraphMK`
        (p.5, Definition 5).

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
        nx, Ex, Lx, Lex = x
        ny, Ey, Ly, Ley = y

        Kv = lambda x, y: self._kv(x, y, Lx, Ly)
        Ke = lambda x, y: self._ke(x, y, Ex, Ey, Lex, Ley)

        # initialise cost function
        c = dict()

        # Calculate valid vertices
        Vp = set()
        Ep = dict()

        # calculate product graph vertex set
        for i in range(nx):
            for j in range(ny):
                value = Kv(i, j)
                if(value > 0):
                    # add to vertex set
                    Vp.add((i, j))
                    # initialise an empty set for neighbors
                    if (i, j) not in Ep:
                        Ep[(i, j)] = set()
                c[(i, j)] = value
        # calculate product graph valid edges
        for v in list(Vp):
            for w in list(Vp):
                if v[0] == w[0] or v[1] == w[1]:
                    value = 0
                else:
                    value = Ke((v[0], w[0]),
                               (v[1], w[1]))
                    if(value > 0):
                        # add edge
                        Ep[v].add(w)
                    # store value
                c[(v, w)] = value

        return Vp, Ep, c

    def _subgraph_matching_core(self, w, C, P, Ep, c, b):
        """Calculate the core algorithm of the subgraph matching kernel.

        See :cite:`Kriege2012SubgraphMK` (p.5, Algorithm 1).

        Parameters
        ----------
        w : number
            The current weight.

        C : set
            The current clique.

        P : set
            Current vertices.

        c : dict
            Costs of vertices, edges.

        b : sklearn.utils.Bunch
            A mutable object with two fields on containing
            the *value* of the kernel algorithm and one
            *C* the clique.

        Returns
        -------
        None.

        """
        while len(P) > 0:
            v = P.pop()

            w = w*c[v]
            for u in C:
                w = w*c[(u, v)]

            Cp = C | {v}
            b.value += w*self._lw(Cp)

            # Create the new P
            PN = P & Ep[v]

            # apply subgraph matching for the new clique
            self._subgraph_matching_core(w, Cp, PN, Ep, c, b)


"""k = subgraph_matching()
print("fit")
k.fit([({(1, 2), (2, 3), (2, 1), (3, 2)},
       {1: 'N', 2: 'C', 3: 'O'},
       {(1, 2): ('N', 'C'), (2, 1): ('C', 'N'),
        (2, 3): ('C', 'O'), (3, 2): ('O', 'C')})])

print("transform")
print(k.transform([({(1, 2), (2, 3), (3, 4), (3, 5), (5, 6),
                     (2, 1), (3, 2), (4, 3), (5, 3), (6, 5)},
                    {1: 'O', 2: 'C', 3: 'N', 4: 'C', 5: 'C', 6: 'O'},
                    {(1, 2): ('O', 'C'), (2, 3): ('C', 'N'),
                     (3, 4): ('N', 'C'), (3, 5): ('N', 'C'),
                     (5, 6): ('C', 'O'), (2, 1): ('C', 'O'),
                     (3, 2): ('N', 'C'), (4, 3): ('C', 'N'),
                     (5, 3): ('C', 'N'), (6, 5): ('O', 'C')})]))"""
