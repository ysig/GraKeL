"""The sugraph mathing kernel as defined by :cite:`kriege2012subgraph`."""
# Author: Ioannis Siglidis <y.siglidis@gmail.com>
# License: BSD 3 clause
import collections
import warnings

import numpy as np

from numbers import Real

from grakel.kernels import Kernel
from grakel.graph import Graph
from grakel.kernels._c_functions import sm_kernel


# Define default vertex, edge and lambda weight functions
def _dirac(a, b):
    """Calculate the dirac function for labels."""
    return int(a == b)


class SubgraphMatching(Kernel):
    r"""Calculate the subgraph matching kernel.

    See :cite:`kriege2012subgraph`.

    Parameters
    ----------
    k : int, default=5
        The upper bound for the maximum size of subgraphs.

    lw : str, valid_values={"uniform", "increasing", "decreasing", "strong_decreasing"},
    default="uniform" | iterable, size=k+1,
    | callable, num_of_arguments=1, argument_type=int
        The lambda weights applied to the clique sizes.

    kv : function (`vertex_label, `vertex_label`, -> number), or None
    default=:math:`k_{v}^{default}(l(a), l(b))= \delta(l(a), l(b))`
        The kernel function between two vertex_labels.
        If no function is provided, this is equivalent with not taking into account node labels.

    ke : function (`edge_label`, `edge_label` -> number),
    default=:math:`k_{e}^{default}(l(e), l(e'))= \delta(l(e), l(e'))`
        The kernel function between two edge_labels.
        If no function is provided, this is equivalent with not taking into account edge labels.

    Attributes
    ----------
    lambdas_ : np.array, shape=(1, k+1)
        All the lambdas corresponding to all the valid sizes of subgraphs.

    """

    _graph_format = "all"

    def __init__(self, n_jobs=None, verbose=False,
                 normalize=False, k=5, kv=_dirac,
                 ke=_dirac, lw="uniform"):
        """Initialise a `subgraph_matching` kernel."""
        super(SubgraphMatching, self).__init__(
            n_jobs=n_jobs, verbose=verbose, normalize=normalize)

        self.k = k
        self.kv = kv
        self.ke = ke
        self.lw = lw
        self._initialized.update({"k": False, "kv": False, "ke": False, "lw": False})

    def initialize(self):
        """Initialize all transformer arguments, needing initialization."""
        super(SubgraphMatching, self).initialize()
        if not self._initialized["k"]:
            if type(self.k) is not int and self.k < 1:
                raise TypeError('k must be an integer greater-equal than 1')
            self._initialized["k"] = True

        if not self._initialized["kv"]:
            if not callable(self.kv) and self.kv is not None:
                raise TypeError('kv must be callable or None')
            self._initialized["kv"] = True

        if not self._initialized["ke"]:
            if not callable(self.ke) and self.ke is not None:
                raise TypeError('ke must be callable or None')
            self._initialized["ke"] = True

        if not self._initialized["lw"]:
            k = self.k + 1
            not_str_iter = type(self.lw) is not str and \
                isinstance(self.lw, collections.Iterable)
            if not_str_iter:
                lw = list(self.lw)

            if (not_str_iter and len(lw) == self.k and
                    all(isinstance(x, Real) for x in lw)):
                self.lambdas_ = np.array(lw).reshape((1, k))
            elif self.lw == "uniform":
                self.lambdas_ = np.full((1, k), 1.0)
            elif self.lw == "increasing":
                self.lambdas_ = np.arange(1.0,
                                          float(k) + 1.0).reshape(1, k)
            elif self.lw == "decreasing":
                self.lambdas_ = np.full((1, k), 1.0) / \
                                np.arange(1.0, float(k) + 1.0).reshape(1, k)
            elif self.lw == "strong_decreasing":
                self.lambdas_ = np.full((1, k), 1.0) / \
                                np.square(np.arange(1.0, float(k) + 1.0)
                                          ).reshape(1, k)
            elif callable(self.lw):
                try:
                    self.lambdas_ = \
                        np.array([self.lw(i) for i in range(k)]).reshape((1, k))
                except Exception as e:
                    raise TypeError('Incorrect Callable: ' + str(e))
            else:
                raise TypeError('lw can either be str with values '
                                '"uniform", "increasing", "decreasing", '
                                '"strong_decreasing" or an iterable of k+1 '
                                'elements or a callable of one integer '
                                'argument.')

            self._initialized["lw"] = True

    def pairwise_operation(self, x, y):
        """Calculate the `subgraph_matching` kernel.

        See :cite:`kriege2012subgraph`.

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
        tv = sm_kernel(x, y, self.kv, self.ke, self.k)
        return np.dot(self.lambdas_, tv)

    def parse_input(self, X):
        """Parse and create features for the `subgraph_matching` kernel.

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
        out : list
            The extracted adjacency matrices for any given input.

        """
        if not isinstance(X, collections.Iterable):
            raise TypeError('input must be an iterable\n')
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
                    raise TypeError('each element of X must be either a ' +
                                    'graph object or a list with at least ' +
                                    'a graph like object and node, ' +
                                    'edge labels dict \n')
                n = g.nv()
                E = g.get_edge_dictionary()
                L = g.get_labels(purpose="dictionary", return_none=(self.kv is None))
                Le = g.get_labels(purpose="dictionary", label_type="edge",
                                  return_none=(self.ke is None))
                Er = set((a, b) for a in E.keys()
                         for b in E[a].keys() if a != b)

                i += 1
                out.append((n, Er, L, Le))

            if i == 0:
                raise ValueError('parsed input is empty')
            return out


if __name__ == "__main__":
    k = SubgraphMatching()
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
