"""RW-ker. as in :cite:`Kashima2003MarginalizedKB`, :cite:`Grtner2003OnGK`."""

import warnings

import numpy as np

from numpy.linalg import inv
from scipy.linalg import solve_sylvester

from grakel.graph import graph


def random_walk(X, Y, lamda=0.1, method_type="sylvester"):
    """Calculate the random walk kernel.

    See :cite:`Kashima2003MarginalizedKB`, :cite:`Grtner2003OnGK`
    and :cite:`Vishwanathan2006FastCO`.

    Parameters
    ----------
    X,Y : *valid-graph-format*
        The pair of graphs on which the kernel is applied.

    lamda : float
        A lambda factor concerning summation.

    method_type : str, valid_values={"simple", "sylvester"}
        The method to use for calculating random walk kernel:
            + "simple" *Complexity*: :math:`O(|V|^6)`
              (see :cite:`Kashima2003MarginalizedKB`, :cite:`Grtner2003OnGK`)
            + "sylvester" *Complexity*: :math:`O(|V|^3)`
              (see :cite:`Vishwanathan2006FastCO`)

    Returns
    -------
    kernel : number
        The kernel value.

    """
    g_x = graph(X)
    g_y = graph(Y)

    return random_walk_pair(g_x, g_y, lamda, method_type)


def random_walk_pair(Gx, Gy, lamda=0.1, method_type="sylvester"):
    """Calculate the random walk kernel.

    See :cite:`Kashima2003MarginalizedKB`, :cite:`Grtner2003OnGK`
    and in :cite:`Vishwanathan2006FastCO`.

    Parameters
    ----------
    G{x,y} : graph
        The pair of graphs on which the kernel is applied.

    lamda : float
        A lambda factor concerning summation.

    method_type : str, valid_values={"simple", "sylvester"}
        The method to use for calculating random walk kernel:
            + "simple" *Complexity*: :math:`O(|V|^6)`
              (see :cite:`Kashima2003MarginalizedKB`, :cite:`Grtner2003OnGK`)
            + "sylvester" *Complexity*: :math:`O(|V|^3)`
              (see :cite:`Vishwanathan2006FastCO`)

    Returns
    -------
    kernel : number
        The kernel value.

    """
    Gx.desired_format("adjacency")
    Gy.desired_format("adjacency")
    X = Gx.adjacency_matrix
    Y = Gy.adjacency_matrix

    if lamda > 0.5:
        warnings.warn('by increasing lambda the \
            summation series may not converge')

    if(method_type == "simple"):
        # calculate the product graph
        XY = np.kron(X, Y)

        # algorithm presented in [Kashima et al., 2003; Gartner et al., 2003]
        # complexity of O(|V|^6)

        # XY is a square matrix
        s = XY.shape[0]
        Id = np.identity(s)
        k = np.dot(np.dot(np.ones(s), inv(Id - lamda*XY)).T,
                   np.ones(shape=(s)))
    elif(method_type == "sylvester"):
        # algorithm presented in [Vishwanathan et al., 2006]
        # complexity of O(|V|^3)

        X_dimension = X.shape[0]
        Y_dimension = Y.shape[0]

        # For efficiency reasons multiply lambda
        # with the smallest e_{x,y} in dimension
        e_x = np.ones(shape=(X_dimension, 1))
        e_y = np.ones(shape=(1, Y_dimension))

        # Prepare parameters for sylvester equation
        A = Y

        try:
            B = np.divide(inv(X.T), -lamda)
        except np.linalg.LinAlgError as err:
            if 'Singular matrix' in err.message:
                raise ValueError('Adjacency matrix of the first graph \
                                is not invertible')
            else:
                raise

        C = -np.dot(e_x, np.dot(e_y, B))

        try:
            R = solve_sylvester(A, B, C)
        except np.linalg.LinAlgError as err:
            raise ValueError('Solution was not found for the Sylvester \
                    Equation. Check Input!')

        # calculate kernel
        k = - np.sum(np.sum(R, axis=1), axis=0)
    else:
        pass
        # raise exception?
        # such method does not exist

    return k
