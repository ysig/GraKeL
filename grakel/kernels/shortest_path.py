""" This file contains the shortest path kernel
    as defined by [Borgwardt & Kriegel, 2005]
"""

from ..graph import graph
from ..tools import inv_dict

def shortest_path(X, Y, Lx, Ly, algorithm_type="dijkstra"):
    """ This is a function that implements shortest path kernel
        as proposed by [Borgwardt & Kriegel, 2005]

        X, Y: to valid graph types i.e. edge dictionary or adjacency_mat
        algorithm_type: "dijkstra" or "floyd_warshall" or "auto" (for the best concerning current format)
    """
    g_x = graph(X,Lx)
    g_y = graph(Y,Ly)

    return shortest_path_inner(g_x, g_y, algorithm_type)

def shortest_path_inner(g_x, g_y, algorithm_type="dijkstra"):
    """ A function that calculates the inner summation part of the
        shortest path kernel as proposed by [Borgwardt & Kriegel, 2005]

        g_{x,y}: graph type objects
        algorithm_type: "dijkstra" or "floyd_warshall" or "auto" (for the best concerning current format)

        Complexity: O(n^2*min(1,m)^2), where the number of vertices and m the mean maximum
        branching factor of the surjective labeling function of the Graph with minimum product
        over the two (not yet implemented!).
    """
    # calculate shortest path matrix
    S_x, Lx = g_x.build_shortest_path_matrix(algorithm_type)
    S_y, Ly = g_y.build_shortest_path_matrix(algorithm_type)
    kernel = 0
    Ly_inv = inv_dict(Ly)
    dim = S_x.shape[0]
    for i in range(0,dim):
        for j in range(0, dim):
            if(i!=j):
                if((Lx[i] in Ly_inv) and (Lx[j] in Ly_inv)):
                    for k in Ly_inv[Lx[i]]:
                        for m in Ly_inv[Lx[j]]:
                            if(k!=m and S_x[i,j] == S_y[k,m]):
                                kernel +=1

    return kernel
