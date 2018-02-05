"""The sugraph mathing kernel as defined by :cite:`Kriege2012SubgraphMK`."""
import collections
import warnings
import numpy as np

from numbers import Number

from grakel.kernels import kernel
from grakel.graph import Graph

global kv_default, ke_default, lw_default

# Define default vertex, edge and lambda weight functions
k_default = lambda a, b: 1 if a == b else 0


class subgraph_matching(kernel):
    r"""Calculate the subgraph matching kernel.

    See :cite:`Kriege2012SubgraphMK`.

    Parameters
    ----------
    X,Y : *valid-graph-format*
        The pair of graphs on which the kernel is applied.

    L{x,y} : dict
        Corresponding graph labels for vertices.

    k : int, default=5
        The upper bound for the maximum size of subgraphs.

    lw : str, valid_values={"uniform", "increasing"
    "decreasing", "strong_decreasing"}, default="uniform" | iterable, size=k+1,
    | callable, num_of_arguments=1, argument_type=int
        The lambda weights applied to the clique sizes.

    kv : function (`vertex_label, `vertex_label`, -> number),
    default=:math:`k_{v}^{default}(l(a), l(b))= \delta(l(a), l(b))`
        The kernel function between two vertex_labels.

    ke : function (`edge_label`, `edge_label` -> number),
    default=:math:`k_{e}^{default}(l(e), l(e'))= \delta(l(e), l(e'))`
        The kernel function between two edge_labels.

    Attributes
    ----------
    _kv : function (`vertex_label, `vertex_label`, -> number),
        The kernel function between two edge_labels.

    _ke : function (`edge_label`, `edge_label` -> number),
        The kernel function between two edge_labels.

    _k : int
        The kernel function between two edge_labels.

    _lambdas : np.array, shape=(1, k+1)
        All the lambdas corresponding to all the valid sizes of subgraphs.

    """

    _graph_format = "all"

    def __init__(self, **kargs):
        """Initialise a subtree_wl kernel."""
        self._valid_parameters |= {"kv", "ke", "lw"}

        super(subgraph_matching, self).__init__(**kargs)
        self._kv = kargs.get("kv", k_default)
        self._ke = kargs.get("ke", k_default)

        self._k = kargs.get("k", 5)
        if self._k < 1:
            raise ValueError('k must be greater-equal than 1')

        lw = kargs.get("lw", "uniform")

        self._k += 1
        if lw == "uniform":
            self._lambdas = np.full((1, self._k), 1.0)
        elif lw == "increasing":
            self._lambdas = np.reshape(np.arange(1.0, float(self._k) + 1.0),
                                       shape=(1, self._k))
        elif lw == "decreasing":
            self._lambdas = np.full((1, self._k), 1.0) / \
                            np.reshape(np.arange(1.0, float(self._k) + 1.0),
                                       shape=(1, self._k))
        elif lw == "strong_decreasing":
            self._lambdas = np.full((1, self._k), 1.0) / \
                            np.reshape(np.square(
                                np.arange(1.0, float(self._k) + 1.0)),
                                       shape=(1, self._k))
        elif isinstance(lw, collections.Iterable) is list \
                and len(lw) == self._k \
                and all(isinstance(x, Number.Real) for x in lw):
            np.reshape(np.array(lw), shape=(1, self._k))
        elif callable(lw):
            try:
                np.reshape(np.array([lw(i) for i in range(self._k)]),
                           shape=(1, self._k))
            except Exception as e:
                raise ValueError('Incorrect Callable: ' + str(e))
        else:
            raise ValueError('lw can either be str with values ' +
                             '"uniform", "increasing", "decreasing", ' +
                             '"strong_decreasing" or an iterable of k+1 ' +
                             'elements or a callable of one integer argument.')
        self._k -= 1

    def pairwise_operation(self, x, y):
        """Calculate the subgraph matching kernel.

        See :cite:`Kriege2012SubgraphMK`.

        Parameters
        ----------
        x, y : tuples
            *Vertex-set*, *edge-dictionary*, *node-label-dictionary*,
            *edge-labels-dictionary* tuple.

        Returns
        -------
        kernel : number
            The kernel value.

        """
        # Calculate product graph
        Vp, n, c = self._weighted_product_graph(x, y)

        # Initialize values
        tv = np.zeros(shape=(self._k + 1, 1))

        self._subgraph_matching_core(1, list(), Vp, c, tv, 0, n-1)
        return np.dot(self._lambdas, tv)

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
            for (idx, x) in enumerate(iter(X)):
                if type(x) is Graph:
                    g = Graph(x.get_adjacency_matrix(),
                              x.get_labels(purpose="adjacency"),
                              x.get_labels(purpose="adjacency",
                                           label_type="edge"),
                              self._graph_format)
                elif isinstance(x, collections.Iterable) and \
                        len(x) in [0, 3]:
                    if len(x) == 0:
                        warnings.warn('Ignoring empty element' +
                                      ' on index: '+str(idx))
                        continue
                    elif len(x) == 3:
                        g = Graph(x[0], x[1], x[2], "adjacency")
                        g.change_format(self._graph_format)
                else:
                    raise ValueError('each element of X must be either a ' +
                                     'graph object or a list with at least ' +
                                     'a graph like object and node, ' +
                                     'edge labels dict \n')
                n = g.nv()
                E = g.get_edge_dictionary()
                L = g.get_labels(purpose="dictionary")
                Le = g.get_labels(purpose="dictionary",
                                  label_type="edge")
                Er = set((a, b) for a in E.keys()
                         for b in E[a].keys() if a != b)

                i += 1
                out.append((n, Er, L, Le))

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
        x, y : tuples, size=4
            A tuple corresponding to the number of verices,
            the edge dictionarie starting from vertices of
            index zero, the labels for nodes, the labels for edges.

        Returns
        -------
        Vp : dict
            Enumeration of all the vertices of the weighted product graph.

        Ep : dict
            The edges of the weighted product graph.

        n : int
            The number of vertices.

        c : dict
            The cost dictionary for vertices and edges.

        """
        nx, Ex, Lx, Lex = x
        ny, Ey, Ly, Ley = y

        # initialise cost function
        c = dict()

        # Calculate valid vertices
        Vp = list()

        # calculate product graph vertex set
        nv = 0
        for i in range(nx):
            for j in range(ny):
                value = self._kv(Lx[i], Ly[j])
                if(value > .0):
                    # add to vertex set
                    Vp.append((i, j))
                    # initialise an empty set for neighbors
                    c[nv] = value
                    nv += 1

        # calculate product graph valid edges
        for (i, v) in enumerate(Vp):
            for (j, w) in enumerate(Vp):
                if i == j:
                    break
                if v[0] == w[0] or v[1] == w[1]:
                    value = .0
                else:
                    ea, eb = (v[0], w[0]), (v[1], w[1])
                    conda, condb = ea not in Ex, eb not in Ey
                    if conda and condb:
                        # d-edge
                        value = 1.
                    elif conda or condb:
                        value = .0
                    else:
                        # possible c-edge
                        value = self._ke(Lex[ea], Ley[eb])
                c[(j, i)] = c[(i, j)] = value

        return np.arange(nv), nv, c

    def _subgraph_matching_core(self, w, C, P, c, total_value, lbound, ubound):
        """Calculate the core algorithm of the subgraph matching kernel.

        See :cite:`Kriege2012SubgraphMK` (p.5, Algorithm 1).

        Parameters
        ----------
        w : number
            The current weight.

        C : list
            The current clique.

        P : np.array, shape=(nv)
            Current vertices.

        c : dict
            Costs of vertices, edges.

        total_value : np.array, shape=(self._k+1, 1)
            Holds values for every level of the cliques.

        lbound : int
            The lower bound on subgraph size.

        ubound : int
            The upper bound on subgraph size.

        Returns
        -------
        None.

        """
#       +-----------------------------+----------------+-----+
#     p |                             |   candidates   |     |
#       +-----------------------------+----------------+-----+
#                                  lBound           uBound
        for it in range(lbound, ubound + 1):
            v = P[it]

            w *= c[v]
            for u in C:
                w *= c[(u, v)]

            total_value[len(C)] += w

            if (len(C) + 1 < self._k):
                C.append(v)

                # prepare candidate set for recursive call
                newUBound = ubound
                while (c[P[newUBound]] == .0 and newUBound - 1 > it):
                    newUBound -= 1

                fm = it + 1
                # fm  is newUBound or the first element that
                # should stay in the candidate set
                while (fm <= newUBound and c[P[fm]] == .0):
                    fm += 1

                nm = fm + 1
                while (nm <= newUBound):
                    # nm is > fm and all elements in [fm, ..., nm]
                    # should stay in the candidate set
                    if c[P[nm]] == .0:
                        P[fm], P[nm] = P[nm], P[fm]
                        fm += 1
                    nm += 1

                # [fm, ..., newUBound] contains the new candidate set
                self._subgraph_matching_core(w, C, P, c,
                                             total_value, fm, newUBound)
                C.pop(-1)
        return None


if __name__ == "__main__":
    k = subgraph_matching()
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
                         (5, 3): ('C', 'N'), (6, 5): ('O', 'C')})]))
