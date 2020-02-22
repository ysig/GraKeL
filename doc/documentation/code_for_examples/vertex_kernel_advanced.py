from collections import Counter, Iterable
from sklearn.utils.validation import check_is_fitted
from grakel import Kernel, Graph
from numpy import zeros, einsum


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
        if not self.initialized_["n_jobs"]:
            # n_jobs is not used in this kernel
            # numpy utilises low-level parallelization for calculating matrix operations
            if self.n_jobs is not None:
                warnings.warn('no implemented parallelization for VertexHistogram')
            self.initialized_["n_jobs"] = True

        # for i=1 .. m
        #     if not self.initialized_["param_needing_initialization_i"]:
        #         # Apply checks (raise ValueError or TypeError accordingly)
        #         # calculate derived fields stored on self._derived_field_ia .. z
        #         self.initialized_["param_needing_initialization_i"] = True


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
        out : np.array, shape=(len(X), n_labels)
            A numpy array for frequency (cols) histograms for all Graphs (rows).

        """
        if not isinstance(X, Iterable):
            raise TypeError('input must be an iterable\n')
        else:
            rows, cols, data = list(), list(), list()
            if self._method_calling in [1, 2]:
                labels = dict()
                self._labels = labels
            elif self._method_calling == 3:
                labels = dict(self._labels)
            ni = 0
            for (i, x) in enumerate(iter(X)):
                is_iter = isinstance(x, Iterable)
                if is_iter:
                    x = list(x)
                if is_iter and len(x) in [0, 2, 3]:
                    if len(x) == 0:
                        warnings.warn('Ignoring empty element' +
                                      ' on index: '+str(i))
                        continue
                    else:
                        # Our element is an iterable of at least 2 elements
                        L = x[1]
                elif type(x) is Graph:
                    # get labels in any existing format
                    L = x.get_labels(purpose="any")
                else:
                    raise TypeError('each element of X must be either a ' +
                                     'graph object or a list with at least ' +
                                     'a graph like object and node labels ' +
                                     'dict \n')

                # construct the data input for the numpy array
                for (label, frequency) in Counter(L.values()).items():
                    # for the row that corresponds to that graph
                    rows.append(ni)

                    # and to the value that this label is indexed
                    col_idx = labels.get(label, None)
                    if col_idx is None:
                        # if not indexed, add the new index (the next)
                        col_idx = len(labels)
                        labels[label] = col_idx

                    # designate the certain column information
                    cols.append(col_idx)

                    # as well as the frequency value to data
                    data.append(frequency)
                ni += 1

            # Initialise the feature matrix
            features = zeros(shape=(ni, len(labels)))
            features[rows, cols] = data

            if ni == 0:
                raise ValueError('parsed input is empty')
            return features

    def _calculate_kernel_matrix(self, Y=None):
        """Calculate the kernel matrix given a target_graph and a kernel.

        Each a matrix is calculated between all elements of Y on the rows and
        all elements of X on the columns.

        Parameters
        ----------
        Y : np.array, default=None
            The array between samples and features.

        Returns
        -------
        K : numpy array, shape = [n_targets, n_inputs]
            The kernel matrix: a calculation between all pairs of graphs
            between targets and inputs. If Y is None targets and inputs
            are the taken from self.X. Otherwise Y corresponds to targets
            and self.X to inputs.

        """
        if Y is None:
            K = self.X.dot(self.X.T)
        else:
            K = Y[:, :self.X.shape[1]].dot(self.X.T)
        return K

    def diagonal(self):
        """Calculate the kernel matrix diagonal of the fitted data.

        Parameters
        ----------
        None.

        Returns
        -------
        X_diag : np.array
            The diagonal of the kernel matrix, of the fitted. This consists
            of each element calculated with itself.


        """
        # Check is fit had been called
        check_is_fitted(self, ['X', '_Y'])
        try:
            check_is_fitted(self, ['_X_diag'])
        except NotFittedError:
            # Calculate diagonal of X
            self._X_diag = einsum('ij,ij->i', self.X, self.X)


        Y_diag = einsum('ij,ij->i', self._Y, self._Y)

        return self._X_diag, Y_diag
