.. _subgraph_matching:

Subgraph Matching Kernel
========================

The subgraph matching kernel counts the number of matchings between subgraphs of bounded size in two graphs :cite:`kriege2012subgraph`.
The kernel is very general since it can be applied to graphs that contain node labels, edge labels, node attributes or edge attributes.

Let :math:`\mathcal{G}` be a set of graphs.
We assume that the graphs that are contained in the set are labeled or attributed.
Specifically, let :math:`\ell : \mathcal{V} \cup \mathcal{E} \rightarrow \mathcal{L}` be a labeling function that assigns either discrete labels or continuous attributes to vertices and edges.
A graph isomorphism between two labeled/attributed graphs :math:`G=(V,E)` and :math:`G'=(V',E')` is a bijection :math:`\phi : V \rightarrow V'` that preserves adjacencies, \ie :math:`\forall v,u \in V : (v,u) \in E \Leftrightarrow (\phi(v), \phi(u)) \in E'`, and labels, \ie if :math:`\psi \in V \times V \rightarrow V' \times V'` is the mapping of vertex pairs implicated by the bijection :math:`\phi` such that :math:`\psi((v,u)) = (\phi(v), \phi(u))`, then, the conditions :math:`\forall v \in V : \ell(v) \equiv \ell(\phi(v))` and :math:`\forall e \in E : \ell(e) \equiv \ell(\psi(e))` must hold, where :math:`\equiv` denotes that two labels are considered equivalent.

Given two graphs :math:`G=(V,E)` and :math:`G'=(V',E')`, let :math:`\mathcal{B}(G,G')` denote the set of all bijections between sets :math:`S \subseteq V` and :math:`S' \subseteq V'`, and let :math:`\lambda : \mathcal{B}(G,G') \rightarrow \mathbb{R}^+` be a weight function.
The subgraph matching kernel is defined as

.. math::

    k(G, G') = \sum_{\phi \in \mathcal{B}(G,G')} \lambda(\phi) \prod_{v \in S} \kappa_V(v, \phi(v)) \prod_{e \in S \times S} \kappa_E(e, \psi(e))

where :math:`S = dom(\phi)` and :math:`\kappa_V, \kappa_E` are kernel functions defined on vertices and edges, respectively.

The instance of the subgraph matching kernel that is obtained if we set the :math:`\kappa_V, \kappa_E` functions as follows

.. math::

    \begin{split}
        \kappa_V(v,v') &= \begin{cases}
        1, & \text{if } \ell(v) \equiv \ell(v'),\\
        0, & \text{otherwise and} 
        \end{cases}\\
        \kappa_E(e,e') &= \begin{cases}
        1, & \text{if } e \in E \wedge e' \in E' \wedge \ell(e) \equiv \ell(e') \text{ or } e \not \in E \wedge e' \not \in E',\\
        0, & \text{otherwise.}
        \end{cases}
    \end{split}

is known as the common subgraph isomorphism kernel.
This kernel counts the number of isomorphic subgraphs contained in two graphs.

To count the number of isomorphisms between subgraphs, the kernel capitalizes on a classical result of Levi :cite:`levi1973note` which makes a connection between common subgraphs of two graphs and cliques in their product graph.
More specifically, each maximum clique in the product graph is associated with a maximum common subgraph of the factor graphs.
This allows someone to compute the common subgraph isomorphism kernel by enumerating the cliques of the product graph.

The general subgraph matching kernel extends the theory of Levi and builds a weighted product graph to allow a more flexible scoring of bijections.
Given two graphs :math:`G=(V,E)`, :math:`G'=(V',E')`, and vertex and edge kernels :math:`\kappa_V` and :math:`\kappa_E`, the weighted product graph :math:`G_P=(V_P, E_P)` of :math:`G` and :math:`G'` is defined as

.. math::

    \begin{split}
        V_P &= \{ (v,v') \in V \times V' : \kappa_V(v,v') > 0 \} \\
        E_P &= \{ ((v,v'),(u,u')) \in V_P \times V_P : v \neq u \wedge v' \neq u' \wedge \kappa_E((v,v'),(u,u')) > 0 \} \\
        c(u) &= \kappa_V(v,v') \quad \forall u=(v,v') \in V_P \\
        c(e) &= \kappa_E((v,u),(v',u')) \quad \forall e \in E_P, \\
        \text{where } &e=((v,v'),(u,u')) 
    \end{split}

After creating the weighted product graph, the kernel enumerates its cliques.
The kernel starts from an empty clique and extends it stepwise by all vertices preserving the clique property.
Let :math:`w` be the weight of a clique :math:`C`.
Whenever the clique :math:`C` is extended by a new vertex :math:`v`, the weight of the clique is updated as follows: first it is multiplied by the weight of the vertex :math:`w' = w \cdot c(v)`, and then, it is multiplied by all the edges connecting :math:`v` to a vertex in :math:`C`, that is :math:`w' = \sum_{u \in C} w \cdot c((v,u))`.
The algorithm effectively avoids duplicates by removing a vertex from the candidate set after all cliques containing it have been exhaustively explored.

The runtime of the subgraph matching kernel depends on the number of cliques in the product graph.
The worst-case runtime complexity of the kernel when considering subgraphs of size up to :math:`k` is :math:`\mathcal{O}(kn^{k+1})`, where :math:`n=|V|+|V'|` is the sum of the number of vertices of the two graphs.

The implementation of the subgraph matching kernel can be found below

.. currentmodule:: grakel

.. autosummary::

   SubgraphMatching

Bibliography
------------
.. bibliography:: graph_kernels.bib
   :filter: docname in docnames
