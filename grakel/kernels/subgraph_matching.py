"""The sugraph mathing kernel as defined by :cite:`Kriege2012SubgraphMK`."""
import collections
import warnings

import numpy as np

from numbers import Real

from grakel.kernels import kernel
from grakel.graph import Graph
from grakel.kernels._c_functions import sm_kernel

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

        self._k = kargs.get("k", 3)
        if self._k < 1:
            raise ValueError('k must be greater-equal than 1')

        lw = kargs.get("lw", "uniform")

        self._k += 1

        if type(lw) is not str and isinstance(lw, collections.Iterable):
            lw = list(lw)

        if lw == "uniform":
            self._lambdas = np.full((1, self._k), 1.0)
        elif lw == "increasing":
            self._lambdas = np.arange(1.0,
                                      float(self._k) + 1.0).reshape(1, self._k)
        elif lw == "decreasing":
            self._lambdas = np.full((1, self._k), 1.0) / \
                            np.arange(1.0,
                                      float(self._k) + 1.0).reshape(1, self._k)
        elif lw == "strong_decreasing":
            self._lambdas = np.full((1, self._k), 1.0) / \
                            np.square(np.arange(1.0,
                                                float(self._k) + 1.0)
                                      ).reshape(1, self._k)
        elif len(lw) == self._k and all(isinstance(x, Real) for x in lw):
            np.array(lw).reshape((1, self._k))
        elif callable(lw):
            try:
                np.array([lw(i)
                          for i in range(self._k)]).reshape((1, self._k))
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
        tv = sm_kernel(x, y, self._kv, self._ke, self._k)
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
                is_iter = False
                if isinstance(x, collections.Iterable):
                    is_iter = True
                    x = list(x)

                if type(x) is Graph:
                    g = Graph(x.get_adjacency_matrix(),
                              x.get_labels(purpose="adjacency"),
                              x.get_labels(purpose="adjacency",
                                           label_type="edge"),
                              self._graph_format)
                elif is_iter and len(x) in [0, 3]:
                    x = list(x)
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
