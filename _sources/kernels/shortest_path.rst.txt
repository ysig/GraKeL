.. _shortest_path:

Shortest Path Kernel
====================
The shortest-path kernel decomposes graphs into shortest paths and compares pairs of shortest paths
according to their lengths and the labels of their endpoints.
The first step of the shortest-path kernel is to transform the input graphs into shortest-paths graphs.  
Given an input graph :math:`G=(V,E)`, we create a new graph :math:`S=(V,E_s)` (i.e. its shortest-path graph). 
The shortest-path graph :math:`S` contains the same set of vertices as the graph from which it originates.
The edge set of the former is a superset of that of the latter, since in the shortest-path graph :math:`S`,
there exists an edge between all vertices which are connected by a walk in the original graph :math:`G`.
To complete the transformation, we assign labels to all the edges of the shortest-path graph :math:`S`.  
The label of each edge is set equal to the shortest distance between its endpoints in the original graph :math:`G`.

Given the above procedure for transforming a graph into a shortest-path graph, the shortest-path kernel is defined as follows.

Definition: Shortest-Path Kernel
--------------------------------
Let :math:`G_i`, :math:`G_j` be two graphs, and :math:`S_i`, :math:`S_j` their corresponding shortest-path graphs.
	
The shortest-path kernel is then defined on :math:`S_i=(V_i,E_i)` and :math:`S_j=(V_j,E_j)` as

.. math::
    :nowrap:

    \begin{equation}
        k(S_i,S_j) = \sum_{e_i \in E_i} \sum_{e_j \in E_j} k_{walk}^{(1)}(e_i, e_j)
    \end{equation}

where :math:`k_{walk}^{(1)}(e_i, e_j)` is a positive semidefinite kernel on edge walks of length :math:`1`.

In labeled graphs, the :math:`k_{walk}^{(1)}(e_i, e_j)` kernel is designed to compare both the lengths
of the shortest paths corresponding to edges :math:`e_i` and :math:`e_j`, and the labels of their endpoint vertices.

Let :math:`e_i = \{v_i, u_i\}` and :math:`e_j = \{v_j, u_j\}`.
Then, :math:`k_{walk}^{(1)}(e_i, e_j)` is usually defined as:

.. math::
    :nowrap:
    
    \begin{equation}
    \begin{split}
        k_{walk}^{(1)}(e_i, e_j) &= k_v(\ell(v_i),\ell(v_j)) \ k_e(\ell(e_i),\ell(e_j)) \ k_v(\ell(u_i),\ell(u_j)) \\
        &+ k_v(\ell(v_i),\ell(u_j)) \ k_e(\ell(e_i),\ell(e_j)) \ k_v(\ell(u_i),\ell(v_j))
    \end{split}
    \end{equation}

where :math:`k_v` is a kernel comparing vertex labels, and :math:`k_e` a kernel comparing shortest path lengths.
Vertex labels are usually compared via a dirac kernel, while shortest path lengths may also be compared via
a dirac kernel or, more rarely, via a brownian bridge kernel :cite:`borgwardt2005shortest`.

In terms of runtime complexity, the shortest-path kernel is very expensive since its computation takes :math:`\mathcal{O}(n^4)` time.


Two versions of this kernel can be found implemented below. The first takes as input graphs with discrete node labels and applies a speed-up technique for faster kernel calculation.

.. currentmodule:: grakel

.. autosummary::

   ShortestPath
   ShortestPathAttr

Bibliography
------------
.. bibliography:: graph_kernels.bib
   :filter: docname in docnames
