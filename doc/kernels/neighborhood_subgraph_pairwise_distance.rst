.. _nspdk:

Neighborhood Subgraph Pairwise Distance Kernel
==============================================
The neighborhood subgraph pairwise distance kernel extracts pairs of rooted subgraphs from each graph whose roots are located at a certain distance from each other, and which contain vertices up to a certain distance from the root :cite:`costa2010fast`.
It then compares graphs based on these pairs of rooted subgraphs.
To avoid isomorphism checking, graph invariants are employed to encode each rooted subgraph.

Let :math:`G=(V,E)` be a graph.
The distance between two vertices :math:`u,v \in V`, denoted :math:`D(u,v)`, is the length of the shortest path between them.
The neighborhood of radius :math:`r` of a vertex :math:`v` is the set of vertices at a distance less than or equal to :math:`r` from :math:`v`, that is :math:`\{ u \in V : D(u,v) \leq r\}`.
Given a subset of vertices :math:`S \subseteq V`, let :math:`E(S)` be the set of edges that have both end-points in :math:`S`.
Then, the subgraph with vertex set :math:`S` and edge set :math:`E(S)` is known as the subgraph induced by :math:`S`.
The neighborhood subgraph of radius :math:`r` of vertex :math:`v` is the subgraph induced by the neighborhood of radius :math:`r` of :math:`v` and is denoted by :math:`N_r^v`.
Let also :math:`R_{r,d}(A_v,B_u,G)` be a relation between two rooted graphs :math:`A_v`, :math:`B_u` and a graph :math:`G=(V,E)` that is true if and only if both :math:`A_v` and :math:`B_u` are in :math:`\{N_r^v : v \in V \}`, where we require :math:`A_v, B_u` to be isomorphic to some :math:`N_r^v` to verify the set inclusion, and that :math:`D(u,v) = d`.
We denote with :math:`R^{-1}(G)` the inverse relation that yields all the pairs of rooted graphs :math:`A_v`, :math:`B_u` satisfying the above constraints.
Hence, :math:`R^{-1}(G)` selects all pairs of neighborhood graphs of radius :math:`r` whose roots are at distance :math:`d` in a given graph :math:`G`.
The neighborhood subgraph pairwise distance kernel utilizes the following kernel

.. math::

    k_{r,d}(G, G') = \sum_{A_v, B_v \in R_{r,d}^{-1}(G)} \quad \sum_{A'_{v'}, B'_{v'} \in R_{r,d}^{-1}(G')} \delta(A_v, A'_{v'}) \delta(B_v, B'_{v'})

where :math:`\delta` is :math:`1` if its input subgraphs are isomorphic, and :math:`0` otherwise.
The above kernel counts the number of identical pairs of neighboring graphs of radius :math:`r` at distance :math:`d` between two graphs.
Then, the neighborhood subgraph pairwise distance kernel is defined as

.. math::

    k(G, G') = \sum_{r=0}^{r^*} \sum_{d=0}^{d^*} \hat{k}_{r,d}(G, G')

where :math:`\hat{k}_{r,d}` is a normalized version of :math:`k_{r,d}`, that is

.. math::

    \hat{k}_{r,d}(G,G') = \frac{k_{r,d}(G,G')}{\sqrt{k_{r,d}(G,G) k_{r,d}(G',G')}}

The above version ensures that relations of all orders are equally weighted regardless of the size of the induced part sets.

The neighborhood subgraph pairwise distance kernel includes an exact matching kernel over two graphs (\ie the :math:`\delta` kernel) which is equivalent to solving the graph isomorphism problem.
Solving the graph isomorphism problem is not feasible.
Therefore, the kernel produces an approximate solution to it instead.
Given a subgraph :math:`G_S` induced by the set of vertices :math:`S`, the kernel computes a graph invariant encoding for the subgraph via a label function :math:`\mathcal{L}^g : \mathcal{G} \rightarrow \Sigma^*`, where :math:`\mathcal{G}` is the set of rooted graphs and :math:`\Sigma^*` is the set of strings over a finite alphabet :math:`\Sigma`.
The function :math:`\mathcal{L}^g` makes use of two other label functions: (:math:`1`) a function :math:`\mathcal{L}^n` for vertices, and (:math:`2`) a function :math:`\mathcal{L}^e` for edges.
The :math:`\mathcal{L}^n` function assigns to vertex :math:`v` the concatenation of the lexicographically sorted list of distance-distance from root-label triplets :math:`\langle D(v,u), D(v,h), \mathcal{L}(u) \rangle` for all :math:`u \in S`, where :math:`h` is the root of the subgraph and :math:`\mathcal{L}` is a function that maps vertices/edges to their label symbol.
Hence, the above function relabels each vertex with a string that encodes the initial label of the vertex, the vertex distance from all other labeled vertices, and the distance from the root vertex.
The :math:`\mathcal{L}^e(u,v)` function assigns to edge :math:`(u,v)` the label :math:`\langle \mathcal{L}^n(u)`, :math:`\mathcal{L}^n(v)`, :math:`\mathcal{L}((u,v)) \rangle`.
The :math:`\mathcal{L}^e(u,v)` function thus annotates each edge based on the new labels of its endpoints, and its initial label, if any.
Finally, the function :math:`\mathcal{L}^g(G_S)` assigns to the rooted graph induced by :math:`S` the concatenation of the lexicographically sorted list of :math:`\mathcal{L}^e(u,v)` for all :math:`(u,v) \in E(S)`.
The kernel then employs a hashing function from strings to natural numbers :math:`H : \Sigma^* \rightarrow \mathbb{N}` to obtain a unique identifier for each subgraph.
Hence, instead of testing pairs of subgraphs for isomorphism, the kernel just checks if the subgraphs share the same identifier.

The computational complexity of the neighborhood subgraph pairwise distance kernel is :math:`\mathcal{O}(|V| |S| |E(S)| \log |E(S)|)` and is dominated by the repeated computation of the graph invariant for each vertex of the graph.
Since this is a constant time procedure, for small values of :math:`d^*` and :math:`r^*`, the complexity of the kernel is in practice linear in the size of the graph.

The implementation of the neighborhood subgraph pairwise distance kernel can be found below

.. currentmodule:: grakel

.. autosummary::

   NeighborhoodSubgraphPairwiseDistance

Bibliography
------------
.. bibliography:: graph_kernels.bib
   :filter: docname in docnames
