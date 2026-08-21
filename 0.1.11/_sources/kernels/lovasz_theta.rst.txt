.. _lovasz_theta:

Lovasz Theta Kernel
===================

The Lovász number :math:`\vartheta(G)` of a graph :math:`G=(V,E)` is a real number that is an upper bound on the Shannon capacity of the graph.
It was introduced by László Lovász in :math:`1979` :cite:`lovasz1979shannon`.
The Lovász number is intimately connected with the notion of orthonormal representations of graphs.
An orthonormal representation of a graph :math:`G` consists of a set of unit vectors :math:`U_G = \{ \mathbf{u}_i \in \mathbb{R}^d : || \mathbf{u}_i || = 1 \}_{i \in V}` where each vertex :math:`i` is assigned a unit vector :math:`\mathbf{u}_i` such that :math:`(i,j) \not \in E \implies \mathbf{u}_i^\top \mathbf{u}_j = 0`.
Specifically, the Lovász number of a graph :math:`G` is defined as

.. math::

    \vartheta(G) = \min_{\mathbf{c}, U_G} \max_{i \in V} \frac{1}{(\mathbf{c}^\top \mathbf{u}_i)^2}

where :math:`\mathbf{c} \in \mathbb{R}^d` is a unit vector and :math:`U_G` is an orthonormal representation of :math:`G`. 
Geometrically, :math:`\vartheta(G)` is defined by the smallest cone enclosing a valid orthonormal representation :math:`U_G`.
The Lovász number :math:`\vartheta(G)` of a graph :math:`G` can be computed to arbitrary precision in polynomial time by solving a semidefinite program.

The Lovász :math:`\vartheta` kernel utilizes the orthonormal representations associated with the Lovász number to compare graphs :cite:`johansson2014global`.
The kernel is applicable only to unlabeled graphs.
Given a collection of graphs, it first generates orthonormal representations for the vertices of each graph by computing the Lovász :math:`\vartheta` number.
Hence, :math:`U_G` is a set that contains the orthonormal representations of :math:`G`.
Let :math:`S \subseteq V` be a subset of the vertex set of :math:`G`.
Then, the Lovász value of the set of vertices :math:`S` is defined as

.. math::

    \vartheta_S(G) = \min_{\mathbf{c}} \max_{i \in S} \frac{1}{(\mathbf{c}^\top \mathbf{u}_i)^2}

where :math:`\mathbf{c} \in \mathbb{R}^d` is a unit vector and :math:`\mathbf{u}_i` is the representation of vertex :math:`i` obtained by computing the Lovász number :math:`\vartheta(G)` of :math:`G`.
The Lovász value of a set of vertices :math:`S` represents the angle of the smallest cone enclosing the set of orthonormal representations of these vertices (\ie subset of :math:`U_G` defined as :math:`\{ \mathbf{u}_i : \mathbf{u}_i \in U_G, i \in S \}`).

The Lovász :math:`vartheta` kernel between two graphs :math:`G, G'` is then defined as follows:

.. math::

    k_{Lo}(G, G') = \sum_{S \subseteq V} \sum_{S' \subseteq V'} \delta(|S|, |S'|) \frac{1}{Z_{|S|}} k(\vartheta_S(G), \vartheta_{S'}(G'))

where :math:`Z_{|S|} = \binom{|V|}{|S|} \binom{|V'|}{|S|}`, :math:`\delta(|S|, |S'|)` is a delta kernel (equal to :math:`1` if :math:`|S|=|S'|`, and :math:`0` otherwise), and :math:`k` is a positive semi-definite kernel between Lovász values (\eg linear kernel, gaussian kernel).

The Lovász :math:`\vartheta` kernel consists of two main steps: (:math:`1`) computing the Lovász number :math:`\vartheta` of each graph and obtaining the associated orthonormal representations, and (:math:`2`) computing the Lovász value for all subgraphs (\ie subsets of vertices :math:`S \subseteq V`) of each graph.
Exact computation of the Lovász :math:`\vartheta` kernel is in most real settings infeasible since it requires computing the minimum enclosing cones of :math:`2^n` sets of vertices.

When dealing with large graphs, it is thus necessary to resort to sampling.
Given a graph :math:`G`, instead of evaluating the Lovász value on all :math:`2^n` sets of vertices, the algorithm evaluates it in on a smaller number of subgraphs :math:`\mathfrak{S} \in 2^V`.
Then, the Lovász :math:`\vartheta` kernel is defined as follows

.. math::

    \hat{k}_{Lo}(G, G') = \sum_{S \subseteq \mathfrak{S}} \sum_{S' \subseteq \mathfrak{S}'} \delta(|S|, |S'|) \frac{1}{\hat{Z}_{|S|}} k(\vartheta_S(G), \vartheta_{S'}(G'))

where :math:`\hat{Z}_{|S|} = |\mathfrak{S}_{|S|}| |\mathfrak{S}'_{|S|}|` and :math:`\mathfrak{S}_{|S|}` denotes the subset of :math:`\mathfrak{S}` consisting of all sets of cardinality :math:`|S|`, that is :math:`\mathfrak{S}_{|S|} = \{ B \in \mathfrak{S} : |B| = |S| \}`.

The time complexity of computing :math:`\hat{k}_{Lo}(G, G')` is :math:`\mathcal{O}(n^2 m \epsilon^{-1} + s^2 T(k) + sn)` where :math:`T(k)` is the complexity of computing the base kernel :math:`k`, :math:`n = |V|`, :math:`m = |E|` and :math:`s = \max(|\mathfrak{S}|, |\mathfrak{S}'|)`.
The first term represents the cost of solving the semi-definite program that computes the Lovász number :math:`\vartheta`.
The second term corresponds to the worst-case complexity of computing the sum of the Lovász values.
And finally, the third term is the cost of computing the Lovász values of the sampled subsets of vertices.

The implementation of the Lovász :math:`\vartheta` kernel can be found below

.. currentmodule:: grakel

.. autosummary::

   LovaszTheta

Bibliography
------------
.. bibliography:: graph_kernels.bib
   :filter: docname in docnames
