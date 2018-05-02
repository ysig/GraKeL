"""RW-ker. as in :cite:`Kashima2003MarginalizedKB`, :cite:`Grtner2003OnGK`."""
import collections
import warnings

import numpy as np

from itertools import product
from numpy import array

from numpy import ComplexWarning
from numpy.linalg import inv
from numpy.linalg import eig
from scipy.linalg import expm

from grakel.kernels import Kernel
from grakel.graph import Graph

# Python 2/3 cross-compatibility import
from six.moves import filterfalse
from builtins import range


class RandomWalk(Kernel):
    """The random walk kernel class.

    See :cite:`Kashima2003MarginalizedKB`, :cite:`Grtner2003OnGK`
    and :cite:`Vishwanathan2006FastCO`.

    Parameters
    ----------
    lambda : float
        A lambda factor concerning summation.

    method_type : str, valid_values={"baseline", "fast"}
        The method to use for calculating random walk kernel:
            + "baseline" *Complexity*: :math:`O(|V|^6)`
              (see :cite:`Kashima2003MarginalizedKB`, :cite:`Grtner2003OnGK`)
            + "fast" *Complexity*: :math:`O((|E|+|V|)|V||M|)`
              (see :cite:`Vishwanathan2006FastCO`)

    kernel_type : str, valid_values={"geometric", "exponential"}
        Defines how inner summation will be applied.

    p : int or None
        If initialised defines the number of steps.

    Attributes
    ----------
    lamda : float, default=0.1
        A lambda factor concerning summation.

    kernel_type : str, valid_values={"geometric", "exponential"},
    default="geometric"
        Defines how inner summation will be applied.

    method_type : str valid_values={"baseline", "fast"},
    default="fast"
        The method to use for calculating random walk kernel:
            + "baseline" *Complexity*: :math:`O(|V|^6)`
              (see :cite:`Kashima2003MarginalizedKB`, :cite:`Grtner2003OnGK`)
            + "fast" *Complexity*: :math:`O((|E|+|V|)|V||M|)`
              (see :cite:`Vishwanathan2006FastCO`)

    p : int or None
        If initialised defines the number of steps.

    """

    _graph_format = "adjacency"

    def __init__(self, n_jobs=None,
                 normalize=False, verbose=False,
                 lamda=0.1, method_type="fast",
                 kernel_type="geometric", p=None):
        """Initialise a random_walk kernel."""
        # setup valid parameters and initialise from parent
        super(RandomWalk, self).__init__(
            n_jobs=n_jobs, normalize=normalize, verbose=verbose)

        # Ignores ComplexWarning as it does not signify anything problematic
        warnings.filterwarnings('ignore', category=ComplexWarning)

        # Setup method type and define operation.
        self.method_type = method_type
        self.kernel_type = kernel_type
        self.p = p
        self.lamda = lamda
        self.initialized_.update({"method_type": False, "kernel_type": False,
                                  "p": False, "lamda": False})

    def initialize_(self):
        """Initialize all transformer arguments, needing initialization."""
        super(RandomWalk, self).initialize_()

        if not self.initialized_["method_type"]:
            # Setup method type and define operation.
            if self.method_type == "fast":
                def invert(n, w, v):
                    return (np.sum(v, axis=0)/n, w, np.sum(inv(v), axis=1)/n)

                def add_input(x):
                    return invert(x.shape[0], *eig(x))

                self._add_input = add_input
            elif self.method_type == "baseline":
                def add_input(x):
                    return x
                self._add_input = add_input
            else:
                raise ValueError('unsupported method_type')
            self.initialized_["method_type"] = True

        if not self.initialized_["kernel_type"]:
            if self.kernel_type not in ["geometric", "exponential"]:
                raise ValueError('unsupported kernel type: either "geometric" '
                                 'or "exponential"')

        if not self.initialized_["p"]:
            if self.p is not None:
                if type(self.p) is int and self.p > 0:
                    if self.kernel_type == "geometric":
                        self._mu = [1]
                        fact = 1
                        power = 1
                        for k in range(1, self.p + 1):
                            fact *= k
                            power *= self.lamda
                            self._mu.append(fact/power)
                    else:
                        self._mu = [1]
                        power = 1
                        for k in range(1, self.p + 1):
                            power *= self.lamda
                            self._mu.append(power)
                else:
                    raise TypeError('p must be a positive integer bigger than '
                                    'zero or nonetype')
                self.initialized_["kernel_type"] = True

        if not self.initialized_["lamda"]:
            if self.lamda <= 0:
                raise TypeError('lambda must be positive bigger than equal')
            elif self.lamda > 0.5 and self.p is None:
                warnings.warn('ranodm-walk series may fail to converge')
            self.initialized_["lamda"] = True

    def parse_input(self, X):
        """Parse and create features for random_walk kernel.

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
                is_iter = isinstance(x, collections.Iterable)
                if is_iter:
                    x = list(x)
                if is_iter and len(x) in [0, 1, 2, 3]:
                    if len(x) == 0:
                        warnings.warn('Ignoring empty element' +
                                      ' on index: '+str(idx))
                        continue
                    else:
                        A = Graph(x[0], {}, {},
                                  self._graph_format).get_adjacency_matrix()
                elif type(x) is Graph:
                    A = x.get_adjacency_matrix()
                else:
                    raise TypeError('each element of X must be either a ' +
                                    'graph or an iterable with at least 1 ' +
                                    'and at most 3 elements\n')
                i += 1
                out.append(self._add_input(A))

            if i == 0:
                raise ValueError('parsed input is empty')

            return out

    def pairwise_operation(self, X, Y):
        """Calculate the random walk kernel.

        Fast:
        Spectral demoposition algorithm as presented in
        :cite:`Vishwanathan2006FastCO` p.13, s.4.4, with
        complexity of :math:`O((|E|+|V|)|E||V|^2)` for graphs witout labels.

        Baseline:
        Algorithm presented in :cite:`Kashima2003MarginalizedKB`,
        :cite:`Grtner2003OnGK` with complexity of :math:`O(|V|^6)`

        Parameters
        ----------
        X, Y : Objects
            Objects as produced from parse_input.

        Returns
        -------
        kernel : number
            The kernel value.

        """
        if self.method_type == "baseline":
            # calculate the product graph
            XY = np.kron(X, Y)

            # algorithm presented in
            # [Kashima et al., 2003; Gartner et al., 2003]
            # complexity of O(|V|^6)

            # XY is a square matrix
            s = XY.shape[0]
            Id = np.identity(s)

            if self.kernel_type == "geometric":
                return np.linalg.multi_dot(
                    (np.ones(s),
                     inv(Id - self.lamda*XY).T, np.ones(shape=(s))))
            elif self.kernel_type == "exponential":
                return np.linalg.multi_dot((np.ones(s),
                                            expm(self.lamda*XY).T,
                                            np.ones(shape=(s))))

        elif self.method_type == "fast":
            # Spectral demoposition algorithm as presented in
            # [Vishwanathan et al., 2006] p.13, s.4.4, with
            # complexity of O((|E|+|V|)|E||V|^2) for graphs
            # witout labels

            # calculate kernel
            qi_Pi, wi, Pi_inv_pi = X
            qj_Pj, wj, Pj_inv_pj = Y

            # calculate left right flanking factors
            fl = np.kron(qi_Pi, qj_Pj)
            fr = np.kron(Pi_inv_pi, Pj_inv_pj)

            # calculate D based on the method
            Dij = np.kron(wi, wj)
            if self.p is not None:
                Q = np.diagonal(Dij)
                D = np.eye(Q.shape[0])
                S = self._mu[0] * Q
                for k in self._mu[1:]:
                    D *= Q
                    S += k*D

            else:
                if self.kernel_type == "geometric":
                    D = np.diagflat(1/(1-self.lamda*Dij))
                elif self.kernel_type == "exponential":
                    D = np.diagflat(np.exp(self.lamda*Dij))
            return np.linalg.multi_dot((fl, D, fr))


class RandomWalkLabeled(RandomWalk):
    """The labeled random walk kernel class.

    See :cite:`Kashima2003MarginalizedKB_lbld`, :cite:`Grtner2003OnGK_lbld`
    and :cite:`Vishwanathan2006FastCO_lbld`.

    Parameters
    ----------
    lambda : float
        A lambda factor concerning summation.

    method_type : str, valid_values={"baseline", "fast"}
        The method to use for calculating random walk kernel:
            + "baseline" *Complexity*: :math:`O(|V|^6)`
              (see :cite:`Kashima2003MarginalizedKB_lbld`,
              :cite:`Grtner2003OnGK_lbld`)
            + "fast" *Complexity*: :math:`O((|E|+|V|)|V||M|)`
              (see :cite:`Vishwanathan2006FastCO_lbld`)

    kernel_type : str, valid_values={"geometric", "exponential"}
        Defines how inner summation will be applied.

    p : int, optional
        If initialised defines the number of steps.

    Attributes
    ----------
    _lamda : float, default=0.1
        A lambda factor concerning summation.

    _kernel_type : str, valid_values={"geometric", "exponential"},
    default="geometric"
        Defines how inner summation will be applied.

    _method_type : str valid_values={"baseline", "fast"},
    default="fast"
        The method to use for calculating random walk kernel:
            + "baseline" *Complexity*: :math:`O(|V|^6)`
              (see :cite:`Kashima2003MarginalizedKB_lbld`,
              :cite:`Grtner2003OnGK_lbld`)
            + "fast" *Complexity*: :math:`O((|E|+|V|)|V||M|)`
              (see :cite:`Vishwanathan2006FastCO_lbld`)

    _p : int, default=1
        If not -1, the number of steps of the random walk kernel.

    """

    _graph_format = "adjacency"

    def __init__(self, n_jobs=None,
                 normalize=False, verbose=False,
                 lamda=0.1, method_type="fast",
                 kernel_type="geometric", p=None):
        """Initialise a labeled random_walk kernel."""
        # Initialise from parent
        super(RandomWalkLabeled, self).__init__(
            n_jobs=n_jobs, normalize=normalize, verbose=verbose,
            lamda=lamda, method_type=method_type, kernel_type=kernel_type,
            p=p)

    def parse_input(self, X):
        """Parse and create features for graphlet_sampling kernel.

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
                is_iter = isinstance(x, collections.Iterable)
                if is_iter:
                    x = list(x)
                if is_iter and len(x) in [1, 2, 3]:
                    if len(x) == 0:
                        warnings.warn('Ignoring empty element' +
                                      ' on index: '+str(idx))
                        continue
                    else:
                        x = Graph(x[0], x[1], {}, self._graph_format)
                elif type(x) is not Graph:
                    raise TypeError('each element of X must be either a ' +
                                    'graph or an iterable with at least 2 ' +
                                    'and at most 3 elements\n')
                i += 1
                x.desired_format("adjacency")
                out.append((self._add_input(x.get_adjacency_matrix()),
                            x.get_labels(purpose="adjacency")))

            if i == 0:
                raise ValueError('parsed input is empty')

            return out

    def pairwise_operation(self, X, Y):
        """Calculate the labeled random walk kernel.

        Fast:
        Spectral demoposition algorithm as presented in
        :cite:`Vishwanathan2006FastCO_lbld` p.13, s.4.4, with
        complexity of :math:`O((|E|+|V|)|E||V|^2)` for graphs witout labels.

        Baseline:
        Algorithm presented in :cite:`Kashima2003MarginalizedKB_lbld`,
        :cite:`Grtner2003OnGK_lbld` with complexity of :math:`O(|V|^6)`

        Parameters
        ----------
        X, Y : tuples
            Tuples of adjacency matrices and labels.

        Returns
        -------
        kernel : number
            The kernel value.

        """
        X, Lx = X
        Y, Ly = Y

        if self.method_type == "baseline":
            # calculate the labeled kronecker product graph
            XY = [[X[i, l] * Y[j, k] for l, k in filterfalse(
                    lambda x: Lx[x[0]] != Ly[x[1]],
                    product(range(X.shape[1]), range(Y.shape[1])))]
                  for i, j in filterfalse(
                    lambda x: Lx[x[0]] != Ly[x[1]],
                    product(range(X.shape[0]), range(Y.shape[0])))]

            # Check that the matrix is not empty.
            if any(len(l) == 0 for l in XY):
                return .0
            else:
                XY = np.array(XY)

            # algorithm presented in
            # [Kashima et al., 2003; Gartner et al., 2003]
            # complexity of O(|V|^6)

            # XY is a square matrix
            s = XY.shape[0]
            Id = np.identity(s)

            if self.kernel_type == "geometric":
                return np.linalg.multi_dot(
                    (np.ones(s), inv(Id - self.lamda*XY).T,
                     np.ones(shape=(s))))
            elif self.kernel_type == "exponential":
                return np.linalg.multi_dot((np.ones(s),
                                            expm(self.lamda*XY).T,
                                            np.ones(shape=(s))))

        elif self.method_type == "fast":
            # Spectral demoposition algorithm as presented in
            # [Vishwanathan et al., 2006] p.13, s.4.4, with
            # complexity of O((|E|+|V|)|E||V|^2) for graphs
            # witout labels

            # calculate kernel
            qi_Pi, wi, Pi_inv_pi = X
            qj_Pj, wj, Pj_inv_pj = Y

            # calculate left right flanking factors
            label_idx = [(i, j) for i in range(qi_Pi.shape[0])
                         for j in range(qj_Pj.shape[0])]
            fl = array([qi_Pi[i]*qj_Pj[j] for i, j in label_idx])
            fr = array([Pi_inv_pi[i]*Pj_inv_pj[j] for i, j in label_idx])

            # calculate D based on the method
            Dij = array([wi[i]*wj[j] for i, j in label_idx])
            if self.p is not None:
                Q = np.diagonal(Dij)
                D = np.eye(Q.shape[0])
                S = self._mu[0] * Q
                for k in self._mu[1:]:
                    D *= Q
                    S += k*D

            else:
                if self.kernel_type == "geometric":
                    D = np.diagflat(1/(1-self.lamda*Dij))
                elif self.kernel_type == "exponential":
                    D = np.diagflat(np.exp(self.lamda*Dij))
            return np.linalg.multi_dot((fl, D, fr))
