from warnings import warn
from collections import Counter, Iterable
from grakel import kernel, Graph

new_parameters = set()


class vertex_histogram(kernel):
    """Vertex Histogram kernel as found in :cite:`Sugiyama2015NIPS`

    Parameters
    ----------
    None.

    Attributes
    ----------
    None.

    """
    # Define the graph format that this kernel needs (if needed)
    # _graph_format = "auto" (default: "auto")

    def __init__(self, **kargs):
        """Initialise an `odd_sth` kernel."""

        # Add new parameters
        self._valid_parameters |= new_parameters

        super(vertex_histogram, self).__init__(**kargs)

        # Get parameters and check the new ones
        # @for i=1 to num_new_parameters
        # self._param_i = kargs.get("param_i", default_param_i_value)
        # if self._param_i is not type(param_i_type) and not satisfy_constraint(self.param_i, param_i_constraint):
        #     raise ValueError('Param_i should be of @param_i_type type and must satisfy the @param_i_constraint')

    def parse_input(self, X):
        """Parse and check the given input for vertex kernel.

        Parameters
        ----------
        X : iterable
            For the input to pass the test, we must have:
            Each element must be an iterable with at most three features and at
            least one. The first that is obligatory is a valid graph structure
            (adjacency matrix or edge_dictionary) while the second is
            node_labels and the third edge_labels (that fitting the given graph
            format).

        Returns
        -------
        out : list
            List of frequency-histogram for each Graph.

        """
        if not isinstance(X, Iterable):
            raise ValueError('input must be an iterable\n')
        else:
            out = list()
            for (i, x) in enumerate(iter(X)):
                is_iter = isinstance(x, Iterable)
                if is_iter:
                    x = list(x)
                if is_iter and len(x) in [0, 2, 3]:
                    if len(x) == 0:
                        warn('Ignoring empty element on index: '+str(i))
                        continue
                    else:
                        # Our element is an iterable of at least 2 elements
                        labels = x[1]
                elif type(x) is Graph:
                    # get labels in any existing format
                    labels = x.get_labels(purpose="any")
                else:
                    raise ValueError('each element of X must be either a ' +
                                     'graph object or a list with at least ' +
                                     'a graph like object and node labels ' +
                                     'dict \n')

                # Append frequencies for the current Graph
                out.append(Counter(labels.values()))

            if len(out) == 0:
                raise ValueError('parsed input is empty')
            return out

    def pairwise_operation(self, x, y):
        """Calculate sum of frequency products.

        Parameters
        ----------
        x, y : Counter
            Label-Frequency Counters as occur from `parse_input`.

        Returns
        -------
        kernel : number
            The kernel value.

        """
        return sum(x[k]*y[k] for k in x.keys())
