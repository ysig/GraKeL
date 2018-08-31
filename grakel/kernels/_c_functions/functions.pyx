"""A file that wraps c++ functions used in kernels submodule of grakel"""
# Author: Ioannis Siglidis <y.siglidis@gmail.com>
# License: BSD 3 clause"
import numpy as np
import cython

from functools import reduce as freduce
from itertools import combinations
from collections import defaultdict

from numpy import floor, sqrt

cimport numpy as np

from libc.string cimport const_char
from libc.stdlib cimport malloc, free

from grakel.kernels._c_functions.header cimport ArashPartov, sm_core_init

def APHash(word):
    """C++ wrapped implementation of Arash Partov Hashing."""
    bs = word.encode('UTF-8')
    cdef int length = len(word);
    cdef const_char* string = bs;
    return ArashPartov(string, length)


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

    if kv is not None:
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
    else:
        Vp = [(i, j) for i in range(nx) for j in range(ny)]
        nv = nx*ny
        cv_l = iter(1. for _ in range(nv))

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
        if ke is not None:
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
        else:
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
                            value = 1.
                    ce[j][i] = ce[i][j] = value
            

    # Initialize values
    cdef np.ndarray[double, ndim=1] tv_np = np.zeros(shape=(k + 1))
    cdef double *tv = &tv_np[0]

    try:
        # Run the core function
        sm_core_init(1, enum, nv, k, cv, ce, tv)

        tv_np.reshape((k+1, 1))
        return tv_np
    finally:
        # Deallocate memory
        free(enum)
        free(cv)
        for i in range(nv):
            free(ce[i])
        free(ce)

def k_to_ij_triangular(k, dim):
    i = int(dim - 1 - floor(sqrt(-8*k + 4*(dim+1)*dim-7)/2.0 - 0.5))
    j = int(k + i - (dim+1)*(dim)//2 + (dim-i+1)*(dim-i)//2)
    return (i, j)

def k_to_ij_rectangular(k, dim):
    i = k % dim
    j = k // dim
    return (i, j)


# ConSubg from:
# Karakashian, Shant Kirakos et al. “An Algorithm for Generating All Connected Subgraphs with k Vertices of a Graph.” (2013).
def ConSubg(G, k, symmetric):
    # G: dict of sets
    l = set()
    if symmetric:
        sG = G
        for u in G.keys():
            l |= CombinationsWithV(u, k, sG)    
            sGP = dict()
            for v in sG.keys():
                if u != v:
                    sGP[v] = sG[v] - {u}
            sG = sGP
    else:
        for u in G.keys():
            l |= CombinationsWithV(u, k, G)

    return l

def CombinationsWithV(u, k, G_init):
    l = list()
    tree = defaultdict(set)
    treeL = {0: u}
    MarkN = dict()
    def CombinationTree(u, k, G):
        root = u
        l = [set() for i in range(k)]
        l[0].add(u)

        MarkV = dict()
        def BuildTree(nt, depth, k):
            # globals l, MarkN, MarkV, tree
            l[depth] = set(l[depth-1])
            for v in G[treeL[nt]]:
                if v != nt and v not in l[depth]:
                    ntp = len(treeL)
                    treeL[ntp] = v
                    tree[nt].add(ntp)
                    l[depth].add(v)
                    if not MarkV.get(v, False):
                        MarkN[ntp], MarkV[v] = True, True
                    else:
                        MarkN[ntp] = False
                    if depth + 1 <= k-1:
                        BuildTree(ntp, depth + 1, k)

        BuildTree(0, 1, k)

    def unionProduct(S1, S2):
        # globals tree, MarkN
        # print("To compare", S1, S2)
        if not len(S1):
            return set()
        elif not len(S2):
            return {S1}
        else:
            return {s1 | s2 for s1 in S1 for s2 in S2 for s1p, s2p in [({treeL[i] for i in s1}, {treeL[i] for i in s2})] if not len(s1p & {treeL[i] for i in s2}) and (any(MarkN[j] for j in s2) or all(not len({treeL[j] for j in tree[i]} & s2p) for i in s1))}
      
    # Memoization
    CFM = dict()

    def CombinationsFromTree(root, k):
        # Globals tree
        t = root
        lnodesets = set()
        if k == 1:
            return {frozenset({t})}
        for i in range(1, min(len(tree[t]), k - 1) + 1):
            for NodeComb in combinations(tree[t], i):
                for string in compositions(k - 1, i):
                    fail = False
                    S = list()
                    for pos in range(i):
                        stRoot = NodeComb[pos]
                        size = string[pos]
                        m = CFM.get((stRoot, size), None)
                        if m is None:
                            m = CFM[stRoot, size] = CombinationsFromTree(stRoot, size)
                        
                        S.append(m)
                        if not len(S[-1]):
                           fail = True
                           break
                    if fail:
                       continue
                    for combProduct in freduce(unionProduct, S):
                        lnodesets.add(frozenset(combProduct | {t}))
        return lnodesets

    CombinationTree(u, k, G_init)
    return {frozenset({treeL[f] for f in fs}) for fs in CombinationsFromTree(0, k)}

def compositions(n, k):
  if n < 0 or k < 0:
    return
  elif k == 0:
    if n == 0:
      yield []
    return
  elif k == 1:
    yield [n]
    return
  else:
    for i in range(1, n):
      for comp in compositions(n-i, k-1):
        yield [i] + comp
