"""This file contains the standard and simple dirac kernel."""

from grakel.graph import graph
from grakel.tools import inv_dict


def dirac(X, Y, Lx, Ly):
    r"""Labelled graph dirac kernel.

    Parameters
    ----------
    X,Y : *valid-graph-format*
        The pair of graphs on which the kernel is applied.

    L{x,y} : dict
        Corresponding graph labels for vertices.

    Returns
    -------
    kernel : number
        The kernel value.
        .. math::
            k(X,Y) = \sum_{v \in V_{1}}\sum_{u \in V_{2}}\delta(l(v),l(u))

    """
    Gx = graph(X, Lx)
    Gy = graph(Y, Ly)
    return dirac_pair(Gx, Gy)


def dirac_pair(Gx, Gy):
    r"""Labelled graph dirac kernel.

    Parameters
    ----------
    G_{x,y} : graph
        The pair of graphs on which the kernel is applied

    Returns
    -------
    kernel : number
        The kernel value.
        .. math::
            k(X,Y) = \sum_{v \in V_{1}}\sum_{u \in V_{2}}\delta(l(v),l(u))

    """
    Gx.desired_format("dictionary")
    Gy.desired_format("dictionary")

    # Calculate kernel
    linv_x = inv_dict(Gx.get_labels(purpose="dictionary"))
    linv_y = inv_dict(Gy.get_labels(purpose="dictionary"))

    kernel = 0
    for lx in linv_x:
        if lx in linv_y:
            kernel += len(linv_x[lx])*len(linv_y[lx])

    return kernel
