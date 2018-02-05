"""The lovasz theta kernel as defined in :cite:`Johansson2015LearningWS`."""

import itertools
import collections
import warnings

from grakel.kernels import kernel
from grakel.graph import Graph


class lovasz_theta(kernel):
    """Lovasz theta kernel as proposed in :cite:`Johansson2015LearningWS`.

    Parameters
    ----------
    X,Y : *valid-graph-format*
        The pair of graphs on which the kernel is applied.

    n_samples : int, default=50
        The number of samples.

    subsets_size_range : tuple, len=2, default=(2,8)
        (min, max) size of the vertex set of sampled subgraphs.

    metric : function (number, number -> number), default=:math:`f(x,y) = x*y`
        The applied metric between the lovasz_theta numbers of the two graphs.

    Attributes
    ----------
    _n_samples : int
        Number of samples drawn for the computation of lovasz theta.

    _ssr : tuple, len=2
        A tuple containing two integers designating the minimum and the maximum
        size of the vertex set of considered subgraphs.

    _metric : function (number, number -> number)
        The applied metric between the lovasz_theta numbers of the two graphs.

    _metric_type : str, fixed="lovasz"
        The type of metric calculated from graphs.

    """

    _metric_type = "lovasz"
    _graph_format = "all"

    def __init__(self, **kargs):
        """Initialise a lovasz_theta kernel."""
        # setup valid parameters and initialise from parent
        self._valid_parameters |= {"n_samples", "subsets_size_range", "metric"}
        super(lovasz_theta, self).__init__(**kargs)

        self._n_samples = kargs.get("n_samples", 50)
        if self._n_samples <= 0 or type(self._n_samples) is not int:
            raise ValueError('n_samples must an integer be bigger than zero')

        self._ssr = kargs.get("subsets_size_range", (2, 8))
        if (type(self._ssr) is not tuple or len(self._ssr) != 2 or
                any(type(i) is not int for i in self._ssr) or
                self._ssr[0] > self._ssr[1]):
            raise ValueError('subsets_size_range subset size range must ' +
                             'be a tuple of two integers in increasing order')

        self._metric = kargs.get("metric", lambda x, y: x*y)
        if not callable(self._metric):
            raise ValueError('metric between arguments must be a function')

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
            The lovasz metrics for the given input.

        """
        if not isinstance(X, collections.Iterable):
            raise ValueError('input must be an iterable\n')
        else:
            i = 0
            out = list()
            for (idx, x) in enumerate(iter(X)):
                if isinstance(x, collections.Iterable) and \
                        len(x) in [0, 1, 2, 3]:
                    if len(x) == 0:
                        warnings.warn('Ignoring empty element ' +
                                      'on index: '+str(idx))
                        continue
                    else:
                        x = Graph(x[0], {}, {}, self._graph_format)
                elif type(x) is not Graph:
                    raise ValueError('each element of X must be either a ' +
                                     'graph or an iterable with at least 1 ' +
                                     'and at most 3 elements\n')
                i += 1
                out.append(
                    x.calculate_subgraph_samples_metric_dictionary(
                        self._metric_type,
                        n_samples=self._n_samples,
                        subsets_size_range=self._ssr))

            if i == 0:
                raise ValueError('parsed input is empty')

            return out

    def pairwise_operation(self, x, y):
        """Lovasz theta kernel as proposed in :cite:`Johansson2015LearningWS`.

        Parameters
        ----------
        x, y : dict
            Subgraph samples metric dictionaries for all levels.

        Returns
        -------
        kernel : number
            The kernel value.

        """
        kernel = 0

        for level in x.keys():
            if level in y and bool(x[level]) and bool(y[level]):
                    Z = len(x[level])*len(y[level])
                    kernel += sum(self._metric(a, b) for (a, b) in
                                  itertools.product(x[level], y[level]))/Z

        return kernel
