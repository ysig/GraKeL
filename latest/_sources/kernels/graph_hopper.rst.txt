.. _graph_hopper:

Graph Hopper Kernel
===================

Given two graphs, the GraphHopper kernel compares shortest paths between pairs of vertices from the two graphs :cite:`feragen2013scalable`.
The kernel takes into account both path lengths and the vertices encountered while ``hopping'' along shortest paths.
The kernel is equivalent to a weighted sum of node kernels.


Let :math:`G=(V,E)` be a graph.
The graph contains either discrete node labels or continuous node attributes.
Let :math:`\ell : \mathcal{V} \rightarrow \mathcal{L}` be a labeling function that assigns either discrete labels or continuous attributes to vertices.
The kernel compares node labels/attributes using a kernel :math:`k_n` (\eg delta kernel in the case of node labels, and linear or gaussian kernel in the case of node attributes).
Given two vertices :math:`v,u \in V`, a path :math:`\pi` from :math:`v` to :math:`u` in :math:`G` is defined as a sequence of vertices

.. math::

    \pi = [v_1, v_2, v_3, \ldots, v_l]

where :math:`v_1 = v`, :math:`v_l = u` and :math:`(v_i, v_{i+1}) \in E` for all :math:`i=1,\ldots,l-1`.
Let :math:`\pi(i) = v_i` denote the :math:`i^{th}` vertex encountered when ``hopping'' along the path.
Denote by :math:`l(\pi)` the weighted length of :math:`\pi` and by :math:`|\pi|` its discrete length, defined as the number of vertices in :math:`\pi`.
The shortest path :math:`\pi_{ij}` from :math:`v_i` to :math:`v_j` is defined in terms of weighted length.
The diameter :math:`\delta(G)` of :math:`G` is the maximal number of nodes in a shortest path in :math:`G`, with respect to weighted path length.

The GraphHopper kernel is defined as a sum of path kernels :math:`k_p` over the families :math:`P, P'` of shortest
paths in :math:`G,G'`

.. math::

    k(G,G') = \sum_{\pi \in P} \sum_{\pi' \in P'} k_p(\pi, \pi')

The path kernel :math:`k_p(\pi, \pi')` is a sum of node kernels :math:`k_n` on vertices simultaneously encountered while simultaneously hopping along paths :math:`\pi` and :math:`\pi'` of equal discrete length, that is

.. math::

    k_p(\pi, \pi') = \begin{cases}
        \sum_{j=1}^{|\pi|} k_n(\pi(j), \pi'(j)), & \text{if $|\pi| = |\pi'|$},\\
        0, & \text{otherwise.} 
        \end{cases}

The :math:`k(G,G')` kernel can be decomposed into a weighted sum of node kernels

.. math::

    k(G,G') = \sum_{v \in V} \sum_{v' \in V'} w(v,v') k_n(v, v')

where :math:`w(v,v')` counts the number of times :math:`v` and :math:`v'` appear at the same hop, or coordinate, :math:`i` of shortest paths :math:`\pi,\pi'` of equal discrete length :math:`|\pi| = |\pi'|`.
We can decompose the weight :math:`w(v,v')` as

.. math::

    w(v,v') = \sum_{j=1}^\delta \sum_{i=1}^\delta | \{ (\pi,\pi') : \pi(i)=v, \pi'(i)=v', |\pi|=|\pi'|=j \} | = \sum_{j=1}^\delta \sum_{i=1}^\delta [\mathbf{M_v}]_{ij} [\mathbf{M_{v'}}]_{ij}

where :math:`\mathbf{M_v}` is a :math:`\delta \times \delta` matrix whose entry :math:`[\mathbf{M_v}]_{ij}` counts how many times :math:`v` appears at the :math:`i^{th}` coordinate of a shortest path in :math:`G` of discrete length :math:`j`, and :math:`\delta = \max(\delta(G), \delta(G'))`.
The components of these matrices can be computed efficiently using recursive message-passing algorithms. 
The total complexity of computing :math:`k(G,G')` is :math:`\mathcal{O}(n^2(m + \log n + d + \delta^2))` where :math:`n` is the number of vertices, :math:`m` is the number of edges and :math:`d` is the dimensionality of the node attributes (:math:`d=1` in the case of discrete node labels).

The implementation of the neighborhood hash kernel can be found below

.. currentmodule:: grakel

.. autosummary::

   GraphHopper

Bibliography
------------
.. bibliography:: graph_kernels.bib
   :filter: docname in docnames
