import numpy as np
import cython

cimport functions
cimport numpy as np

from libc.string cimport const_char
from libc.stdlib cimport malloc, free

def APHash(word):
    """C++ wrapped implementation of Arash Partov Hashing."""
    bs = word.encode('UTF-8')
    cdef int length = len(word);
    cdef const_char* string = bs;
    return functions.ArashPartov(string, length)


def sm_kernel(x, y, kv, ke, k):
    """Calculate the weighted product graph and calculate the sm_kernel.

    For a definition of the weighted product graph see
    :cite:`Kriege2012SubgraphMK` (p.5, Definition 5).

    Parameters
    ----------
    x, y : tuples, size=4
        A tuple corresponding to the number of verices,
        the edge dictionarie starting from vertices of
        index zero, the labels for nodes, the labels for edges.

    kv : function
        A kernel for vertex labels.

    ke : function
        A kernel for edge labels.

    k : int
        The upper bound for the maximum size of subgraphs.

    Returns
    -------
    tv_np : np.array
        An array holding values for all clique sizes from 0 to k.

    """
    nx, Ex, Lx, Lex = x
    ny, Ey, Ly, Ley = y

    # Costs for vertices
    cv_l = list()

    # Calculate valid vertices
    Vp = list()

    # calculate product graph vertex set
    nv = 0
    for i in range(nx):
        for j in range(ny):
            value = kv(Lx[i], Ly[j])
            if(value != .0):
                # add to vertex set
                Vp.append((i, j))
                # initialise an empty set for neighbors
                cv_l.append(value)
                nv += 1

    # Initialise c arrays.
    cdef int *enum
    cdef double *cv
    cdef double **ce
    enum = <int *>malloc(nv*cython.sizeof(int))
    ce = <double **>malloc(nv*cython.sizeof(cython.p_double))
    cv = <double *>malloc(nv*cython.sizeof(cython.double))
    for (i, v) in enumerate(cv_l):
        enum[i] = i
        cv[i] = v
        ce[i] = <double *>malloc(nv*cython.sizeof(cython.double))

    with cython.boundscheck(False):
        # calculate product graph valid edges
        for (i, v) in enumerate(Vp):
            for (j, w) in enumerate(Vp):
                if i == j:
                    ce[j][i] = .0
                    break
                if v[0] == w[0] or v[1] == w[1]:
                    value = .0
                else:
                    ea, eb = (v[0], w[0]), (v[1], w[1])
                    conda, condb = ea not in Ex, eb not in Ey
                    if conda and condb:
                        # d-edge
                        value = -1.
                    elif conda or condb:
                        value = .0
                    else:
                        # possible c-edge
                        try:
                          lea = Lex[ea]
                          leb = Ley[eb]
                        except KeyError as key_error:
                          raise KeyError(key_error +
                                         '\nEdge labels must exist for '
                                         'all edges.')
                        value = ke(lea, leb)

                ce[j][i] = ce[i][j] = value


    # Initialize values
    cdef np.ndarray[double, ndim=1] tv_np = np.zeros(shape=(k + 1))
    cdef double *tv = &tv_np[0]

    try:
        # Run the core function
        functions.sm_core_init(1, enum, nv, k, cv, ce, tv)

        tv_np.reshape((k+1, 1))
        return tv_np
    finally:
        # Deallocate memory
        free(enum)
        free(cv)
        for i in range(nv):
            free(ce[i])
        free(ce)
