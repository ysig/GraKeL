.. _core_framework:

Core Kernel Framework
=====================

The core framework is a tool for increasing the expressive power of graph kernels :cite:`nikolentzos2018degeneracy`.
The framework is not restricted to graph kernels, but can be applied to any graph comparison algorithm.
It capitalizes on the :math:`k`-core decomposition which is capable of uncovering topological and hierarchical properties of graphs.
Specifically, the :math:`k`-core decomposition is a powerful tool for network analysis and it is commonly used as a measure of importance and well connectedness for vertices in a broad spectrum of applications.
The notion of :math:`k`-core was first introduced by Seidman to study the cohesion of social networks :cite:`seidman1983network`.
In recent years, the :math:`k`-core decomposition has been established as a standard tool in many application domains such as in network visualization :cite:`alvarez2006large`, in protein function prediction :cite:`wuchty2005peeling` and in graph clustering :cite:`giatsidis2014corecluster`.

Core Decomposition
------------------

Let :math:`G = (V,E)` be an undirected and unweighted graph.
Let :math:`n` and :math:`m` denote the number of vertices and number of edges, respectively.
Given a subset of vertices :math:`S \subseteq V`, let :math:`E(S)` be the set of edges that have both end-points in :math:`S`.
Then, :math:`G'=(S,E(S))` is the subgraph induced by :math:`S`.
We use :math:`G' \subseteq G` to denote that :math:`G'` is a subgraph of :math:`G`.
The degree of a vertex :math:`v \in S`, :math:`d_{G'}(v)`, is equal to the number of vertices that are adjacent to :math:`v` in :math:`G'`.
Let :math:`G` be a graph and :math:`G'` a subgraph of :math:`G` induced by a set of vertices :math:`S`.
Then, :math:`G'` is defined to be a :math:`k`-core of :math:`G`, denoted by :math:`C_k`, if it is a maximal subgraph of :math:`G` in which all vertices have degree at least :math:`k`.
Hence, if :math:`G'` is a :math:`k`-core of :math:`G`, then :math:`\forall v \in S`, :math:`d_{G'}(v) \geq k`.
Each :math:`k`-core is a unique subgraph of :math:`G`, and it is not necessarily connected.
The core number :math:`c(v)` of a vertex :math:`v` is equal to the highest-order core that :math:`v` belongs to.
In other words, :math:`v` has core number :math:`c(v) = k`, if it belongs to a :math:`k`-core but not to a :math:`(k+1)`-core.
The degeneracy :math:`\delta^*(G)` of a graph :math:`G` is defined as the maximum :math:`k` for which graph :math:`G` contains a non-empty :math:`k`-core subgraph, :math:`\delta^*(G) = \max_{v \in V}c(v)`.
Furthermore, assuming that :math:`\mathcal{C} = \{  C_0, C_1, \ldots, C_{\delta^*(G)} \}` is the set of all :math:`k`-cores, then :math:`\mathcal{C}` forms a nested chain

.. math::

    C_{\delta^*(G)} \subseteq \ldots \subseteq C_1 \subseteq C_0 = G

Therefore, the :math:`k`-core decomposition is a very useful tool for discovering the hierarchical structure of graphs.
The :math:`k`-core decomposition of a graph can be computed in :math:`\mathcal{O}(n+m)` time \cite{matula1983smallest,batagelj2011fast}. 
The underlying idea is that we can obtain the :math:`i`-core of a graph if we recursively remove all vertices with degree less than :math:`i` and their incident edges from the graph until no other vertex can be removed.


Core Kernels
------------

The :math:`k`-core decomposition builds a hierarchy of nested subgraphs, each having stronger connectedness properties compared to the previous ones.
The core framework measures the similarity between the corresponding according to the hierarchy subgraphs and aggregates the results.
Let :math:`G=(V,E)` and :math:`G'=(V',E')` be two graphs.
Let also :math:`k` be any kernel for graphs.
Then, the core variant of the base kernel :math:`k` is defined as

.. math::

    k_c(G, G') = k(C_0,C'_0) + k(C_1,C'_1) + \ldots + k(C_{\delta^*_{min}},C'_{\delta^*_{min}}) 

where :math:`\delta^*_{min}` is the minimum of the degeneracies of the two graphs, and :math:`C_0,C_1,\ldots,C_{\delta^*_{min}}` and :math:`C'_0,C'_1,\ldots,C'_{\delta^*_{min}}` are the :math:`0`-core, :math:`1`-core,:math:`\ldots`, :math:`\delta^*_{min}`-core subgraphs of :math:`G` and :math:`G'`, respectively.
By decomposing graphs into subgraphs of increasing importance, the algorithm is capable of more accurately capturing their underlying structure.

The computational complexity of the core framework depends on the complexity of the base kernel and the degeneracy of the graphs under comparison.
Given a pair of graphs :math:`G, G'` and an algorithm :math:`A` for comparing the two graphs, let :math:`\mathcal{O}_A` be the time complexity of algorithm :math:`A`.
Let also :math:`\delta^*_{min} = \min \big( \delta^*(G),\delta^*(G') \big)` be the minimum of the degeneracies of the two graphs.
Then, the complexity of computing the core variant of algorithm :math:`A` is :math:`\mathcal{O}_{c}=\delta^*_{min}\mathcal{O}_A`.
It is well-known that the degeneracy of a graph is upper bounded by the maximum of the degrees of its vertices and by the largest eigenvalue of its adjacency matrix :math:`\lambda_1`.
Since in most real-world graphs it holds that :math:`\lambda_1 \ll n`, it also holds that :math:`\delta^*_{max} \ll n`, and hence, the time complexity added by the core framework is not very high.

The implementation of the core framework can be found below

.. currentmodule:: grakel

.. autosummary::

   CoreFramework

Bibliography
------------
.. bibliography:: graph_kernels.bib
   :filter: docname in docnames
