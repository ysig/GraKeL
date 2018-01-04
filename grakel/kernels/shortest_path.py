""" This file contains the shortest path kernel as defined in :cite:`Borgwardt2005ShortestpathKO`.
"""

from ..graph import graph
from ..tools import inv_dict

def shortest_path(X, Y, Lx, Ly, algorithm_type="dijkstra"):
    """ This is a function that implements shortest path kernel as proposed in :cite:`Borgwardt2005ShortestpathKO`.

    Parameters
    ----------
    X,Y : *valid-graph-format*
        The pair of graphs on which the kernel is applied.
        
    algorithm_type : str, default={"dijkstra", "floyd_warshall", "auto"}
        Apply the dijkstra or floyd_warshall algorithm for calculating shortest path
        ot for the best concerning current graph format ("auto").

    Returns
    -------
    kernel : number
        The kernel value.
    """
    g_x = graph(X,Lx)
    g_y = graph(Y,Ly)

    return shortest_path_inner(g_x, g_y, algorithm_type)

def shortest_path_inner(g_x, g_y, algorithm_type="dijkstra"):
    """ A function that calculates the inner summation part of the shortest path kernel as proposed in :cite:`Borgwardt2005ShortestpathKO`.

    Parameters
    ----------
    g_{x,y} : graph
        The pair graphs on which the kernel is applied.
        
    algorithm_type : str, default={"dijkstra", "floyd_warshall", "auto"}
        Apply the dijkstra or floyd_warshall algorithm for calculating shortest path
        ot for the best concerning current graph format ("auto").
    
    Returns
    -------
    kernel : number.
        The kernel value.
        
    Complexity
    ----------
    :math:`O(n^2*min(1,m)^2)`, where the number of vertices is :math:`n` and
    :math:`m` the mean maximum branching factor of the surjective labeling function of the Graph.
    """
    # calculate shortest path matrix
    S_x, Lx = g_x.build_shortest_path_matrix(algorithm_type)
    S_y, Ly = g_y.build_shortest_path_matrix(algorithm_type)
    
    # Initialise
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
