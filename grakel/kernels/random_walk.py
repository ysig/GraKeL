"""RW-ker. as in :cite:`Kashima2003MarginalizedKB`, :cite:`Grtner2003OnGK`."""
import collections
import warnings

import numpy as np

from numpy import ComplexWarning
from numpy.linalg import inv
from numpy.linalg import eig
from scipy.linalg import expm

from grakel.kernels import kernel
from grakel.graph import Graph


class random_walk(kernel):
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

    p : int, optional
        If initialised defines the number of steps.

    Attributes
    ----------
    _lambda : float, default=0.1
        A lambda factor concerning summation.

    _kernel_type : str, valid_values={"geometric", "exponential"},
    default="geometric"
        Defines how inner summation will be applied.

    _method_type : str valid_values={"baseline", "fast"},
    default="fast"
        The method to use for calculating random walk kernel:
            + "baseline" *Complexity*: :math:`O(|V|^6)`
              (see :cite:`Kashima2003MarginalizedKB`, :cite:`Grtner2003OnGK`)
            + "fast" *Complexity*: :math:`O((|E|+|V|)|V||M|)`
              (see :cite:`Vishwanathan2006FastCO`)

    _p : int, default=1
        If not -1, the number of steps of the random walk kernel.

    """

    _graph_format = "adjacency"

    _method_type = "fast"
    _kernel_type = "geometric"
    _p = -1
    _lambda = 0.1

    def __init__(self, **kargs):
        """Initialise a random_walk kernel."""
        # setup valid parameters and initialise from parent
        self._valid_parameters |= {"lambda", "method_type"}
        super(random_walk, self).__init__(**kargs)

        warnings.filterwarnings('ignore', category=ComplexWarning)

        # Setup method type and define operation.
        method_type = kargs.get("method_type", "fast")
        if method_type == "fast":
            invert = lambda n, w, v: (np.sum(v, axis=0)/n,
                                      w,
                                      np.sum(inv(v), axis=1)/n)
            self._add_input = lambda x: invert(x.shape[0], *eig(x))
        else:
            self._add_input = lambda x: x

        self._kernel_type = kargs.get("kernel_type", "geometric")
        if self._kernel_type not in ["geometric", "exponential"]:
            raise ValueError('unsupported kernel type: either "geometric" ' +
                             'or "exponential"')

        if "p" in kargs:
            if kargs["p"] > 0:
                self._p = kargs["p"]
                if self._kernel_type == "geometric":
                    self._mu = [1]
                    fact = 1
                    power = 1
                    for k in range(1, self._p + 1):
                        fact *= k
                        power *= self._lambda
                        self._mu.append(fact/power)
                else:
                    self._mu = [1]
                    power = 1
                    for k in range(1, self._p + 1):
                        power *= self._lambda
                        self._mu.append(power)

        self._lambda = kargs.get("lambda", 0.1)
        if self._lambda <= 0:
            raise ValueError('lambda must be positive bigger than equal')
        elif self._lambda > 0.5 and self._p == -1:
                warnings.warn('ranodm-walk series may fail to converge')

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
                    A = x.get_adjacency_matrix()
                elif isinstance(x, collections.Iterable) and \
                        len(x) in [0, 1, 2, 3]:
                    if len(x) == 0:
                        warnings.warn('Ignoring empty element' +
                                      ' on index: '+str(idx))
                        continue
                    else:
                        A = Graph(x[0], {}, {},
                                  self._graph_format).get_adjacency_matrix()
                else:
                    raise ValueError('each element of X must be either a ' +
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
        if(self._method_type == "baseline"):
            # calculate the product graph
            XY = np.kron(X, Y)

            # algorithm presented in
            # [Kashima et al., 2003; Gartner et al., 2003]
            # complexity of O(|V|^6)

            # XY is a square matrix
            s = XY.shape[0]
            Id = np.identity(s)

            if self._kernel_type == "geometric":
                return np.dot(np.dot(np.ones(s), inv(Id - self._lambda*XY)).T,
                              np.ones(shape=(s)))
            elif self._kernel_type == "exponential":
                return np.dot(np.dot(np.ones(s), expm(self._lambda*XY)).T,
                              np.ones(shape=(s)))

        else:
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
            if self._p > 0:
                Q = np.diagonal(Dij)
                D = np.eye(Q.shape[0])
                S = self._mu[0] * Q
                for k in self._mu[1:]:
                    D *= Q
                    S += k*D

            else:
                if self._kernel_type == "geometric":
                    D = np.diagflat(1/(1-self._lambda*Dij))
                elif self._kernel_type == "exponential":
                    D = np.diagflat(np.exp(self._lambda*Dij))
            return np.dot(fl, np.dot(D, fr))
