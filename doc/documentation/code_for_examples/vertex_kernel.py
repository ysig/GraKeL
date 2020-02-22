from warnings import warn
from collections import Counter, Iterable
from grakel import Kernel, Graph


class VertexHistogram(Kernel):
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

    def __init__(self,
                 n_jobs=n_jobs,
                 verbose=False,
                 normalize=False,
                 # kernel_param_1=kernel_param_1_default,
                 # ...
                 # kernel_param_n=kernel_param_n_default,
                 ):
        """Initialise an `odd_sth` kernel."""

        # Add new parameters
        self._valid_parameters |= new_parameters

        super(VertexHistogram, self).__init__(n_jobs=n_jobs, verbose=verbose, normalize=normalize)

        # Get parameters and check the new ones
        # @for i=1 to num_new_parameters
        #   self.kernel_param_i = kernel_param_i

        # self.initialized_.update({
        #    param_needing_initialization_1 : False
        #             ...
        #    param_needing_initialization_m : False
        # })

    def initialize_(self):
        """Initialize all transformer arguments, needing initialization."""
        # If you want to implement a parallelization by your self here is your chance
        # If there is a pairwise operation on the Kernel object there is parallelization is implemented
        # Just run the initialise from father to initialise a joblib Parallel (if n_jobs is not None).
        super(VertexHistogram, self).initialize_()

        # for i=1 .. m
        #     if not self.initialized_["param_needing_initialization_i"]:
        #         # Apply checks (raise ValueError or TypeError accordingly)
        #         # calculate derived fields stored on self._derived_field_ia .. z
        #         self.initialized_["param_needing_initialization_i"] = True
        pass

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
            raise TypeError('input must be an iterable\n')
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
                    raise TypeError('each element of X must be either a ' +
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
