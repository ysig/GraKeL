"""The Edge Histogram kernel as defined in :cite:`sugiyama2015halting`."""
# Author: Ioannis Siglidis <y.siglidis@gmail.com>
# License: BSD 3 clause
from warnings import warn

from collections import Counter
from collections import Iterable

from sklearn.exceptions import NotFittedError
from sklearn.utils.validation import check_is_fitted

from grakel.kernels import Kernel
from grakel.graph import Graph

from numpy import array
from numpy import zeros
from numpy import squeeze
from numpy import einsum
from scipy.sparse import csr_matrix

# Python 2/3 cross-compatibility import
from six import iteritems
from six import itervalues


class EdgeHistogram(Kernel):
    """Edge Histogram kernel as found in :cite:`sugiyama2015halting`.

    Parameters
    ----------
    sparse : bool, or 'auto', default='auto'
        Defines if the data will be stored in a sparse format.
        Sparse format is slower, but less memory consuming and in some cases the only solution.
        If 'auto', uses a sparse matrix when the number of zeros is more than the half of the matrix size.
        In all cases if the dense matrix doesn't fit system memory, I sparse approach will be tried.

    Attributes
    ----------
    None.

    """

    def __init__(self, n_jobs=None, normalize=False, verbose=False, sparse='auto'):
        """Initialize an edge kernel."""
        super(EdgeHistogram, self).__init__(n_jobs=n_jobs, normalize=normalize, verbose=verbose)
        self.sparse = sparse
        self._initialized.update({'sparse': True})

    def initialize(self):
        """Initialize all transformer arguments, needing initialization."""
        if not self._initialized["n_jobs"]:
            if self.n_jobs is not None:
                warn('no implemented parallelization for EdgeHistogram')
            self._initialized["n_jobs"] = True
        if not self._initialized["sparse"]:
            if self.sparse not in ['auto', False, True]:
                TypeError('sparse could be False, True or auto')
            self._initialized["sparse"] = True

    def parse_input(self, X):
        """Parse and check the given input for EH kernel.

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
            A np array for frequency (cols) histograms for all Graphs (rows).

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
                if is_iter and len(x) in [0, 3]:
                    if len(x) == 0:
                        warn('Ignoring empty element on index: '+str(i))
                        continue
                    else:
                        # Our element is an iterable of at least 2 elements
                        L = x[2]
                elif type(x) is Graph:
                    # get labels in any existing format
                    L = x.get_labels(purpose="any", label_type="edge")
                else:
                    raise TypeError('each element of X must be either a ' +
                                    'graph object or a list with at least ' +
                                    'a graph like object and node labels ' +
                                    'dict \n')

                # construct the data input for the numpy array
                for (label, frequency) in iteritems(Counter(itervalues(L))):
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
            if self._method_calling in [1, 2]:
                if self.sparse == 'auto':
                    self.sparse_ = (len(cols)/float(ni * len(labels)) <= 0.5)
                else:
                    self.sparse_ = bool(self.sparse)

            if self.sparse_:
                features = csr_matrix((data, (rows, cols)), shape=(ni, len(labels)), copy=False)
            else:
                # Initialise the feature matrix
                try:
                    features = zeros(shape=(ni, len(labels)))
                    features[rows, cols] = data
                except MemoryError:
                    warn('memory-error: switching to sparse')
                    self.sparse_, features = True, csr_matrix((data, (rows, cols)), shape=(ni, len(labels)), copy=False)

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

        if self.sparse_:
            return K.toarray()
        else:
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
        check_is_fitted(self, ['X', 'sparse_'])
        try:
            check_is_fitted(self, ['_X_diag'])
        except NotFittedError:
            # Calculate diagonal of X
            if self.sparse_:
                self._X_diag = squeeze(array(self.X.multiply(self.X).sum(axis=1)))
            else:
                self._X_diag = einsum('ij,ij->i', self.X, self.X)
        try:
            # If transform has happened return both diagonals
            check_is_fitted(self, ['_Y'])
            if self.sparse_:
                Y_diag = squeeze(array(self._Y.multiply(self._Y).sum(axis=1)))
            else:
                Y_diag = einsum('ij,ij->i', self._Y, self._Y)
            return self._X_diag, Y_diag
        except NotFittedError:
            # Else just return both X_diag
            return self._X_diag
