.. _multiscale_laplacian:

Multiscale Laplacian Kernel
===========================

The multiscale Laplacian graph kernel can handle unlabeled graphs, graphs with discrete node labels and graphs with continuous node attributes :cite:`kondor2016multiscale`.
It takes into account structure in graphs at a range of different scales by building a hierarchy of nested subgraphs.
These subgraphs are compared to each other using another graph kernel, called the feature space laplacian graph kernel.
This kernel is capable of lifting a base kernel defined on the vertices of two graphs to a kernel between the graphs themselves.
Since exact computation of the multiscale laplacian graph kernel is a very expensive operation, the kernel uses a randomized projection procedure  similar to the popular Nystr{\"o}m approximation for kernel matrices :cite:`williams2001using`.

Let :math:`G=(V,E)` be an undirected graph such that :math:`n = |V|`.
The Laplacian of :math:`G` is a :math:`n \times n` matrix defined as

.. math::

    \mathbf{L} = \mathbf{D} - \mathbf{A} 

where :math:`\mathbf{A}` is the adjacency matrix of :math:`G` and :math:`\mathbf{D}` is a diagonal matrix such that :math:`\mathbf{D}_{ii} = \sum_j \mathbf{A}_{ij}`.

Given two graphs :math:`G_1` and :math:`G_2` of :math:`n` vertices, we can define the kernel between them to be a kernel between the corresponding normal distributions :math:`p_1 = \mathcal{N}(\mathbf{0}, \mathbf{L_1}^{-1})` and :math:`p_2 = \mathcal{N}(\mathbf{0}, \mathbf{L_2}^{-1})` where :math:`\mathbf{0}` is the :math:`n`-dimensional all-zeros vector.
More specifically, given two graphs :math:`G_1` and :math:`G_2` of :math:`n` vertices with Laplacians :math:`\mathbf{L_1}` and :math:`\mathbf{L_2}` respectively, the Laplacian graph kernel with parameter :math:`\gamma` between the two graphs is

.. math::

    k_{LG}(G_1, G_2) = \frac{| (\frac{1}{2} \mathbf{S}_1^{-1} + \frac{1}{2} \mathbf{S}_2^{-1} )^{-1} |^{1/2}}{|\mathbf{S}_1|^{1/4} |\mathbf{S}_2|^{1/4}} 

where :math:`\mathbf{S}_1 = \mathbf{L}_1^{-1} + \gamma \mathbf{I}`, :math:`\mathbf{S}_2 = \mathbf{L}_2^{-1} + \gamma \mathbf{I}` and 
:math:`\mathbf{I}` is the :math:`n \times n` identity matrix.
The Laplacian graph kernel captures similarity between the overall shapes of the two graphs.
However, it assumes that both graphs have the same size, and it is not invariant to permutations of the vertices.

To achieve permutation invariance, the multiscale Laplacian graph kernel represents each vertex as an :math:`m`-dimensional vector whose components correspond to local and permutation invariant vertex features.
Such features may include for instance the degree of the vertex or the number of triangles in which it participates.
Then, it performs a linear transformation and represents each graph as a distribution of the considered features instead of a distribution of its vertices.
Let :math:`\mathbf{U}_1, \mathbf{U}_2 \in \mathbb{R}^{m \times n}` be the feature mapping matrices of the two graphs, that is the matrices whose columns contain the vector representations of the vertices of the two graphs. 
Then, the feature space Laplacian graph kernel is defined as

.. math::

    k_{FLG}(G_1, G_2) = \frac{| (\frac{1}{2} \mathbf{S}_1^{-1} + \frac{1}{2} \mathbf{S}_2^{-1} )^{-1} |^{1/2}}{|\mathbf{S}_1|^{1/4} |\mathbf{S}_2|^{1/4}} 

where :math:`\mathbf{S}_1 = \mathbf{U}_1 \mathbf{L}_1^{-1} \mathbf{U}_1^\top + \gamma \mathbf{I}`, :math:`\mathbf{S}_2 = \mathbf{U}_2 \mathbf{L}_2^{-1} \mathbf{U}_2^\top + \gamma \mathbf{I}` and :math:`\mathbf{I}` is the :math:`m \times m` identity matrix.
Since the vertex features are local and invariant to vertex reordering, the feature space Laplacian graph kernel is permutation invariant.
Furthermore, since the distributions now live in the space of features rather than the space of vertices, the feature space Laplacian graph kernel can be applied to graphs of different sizes.

Let :math:`\phi(v)` be the representation of vertex :math:`v` constructed from local vertex features as described above.
The base kernel :math:`\kappa` between two vertices :math:`v_1` and :math:`v_2` corresponds to the dot product of their feature vectors

.. math::

    \kappa(v_1, v_2) = \phi(v_1)^\top \phi(v_2) 

Let :math:`G_1` and :math:`G_2` be two graphs with vertex sets :math:`V_1 = \{ v_1, \ldots, v_{n_1}\}` and :math:`V_2 = \{ u_1, \ldots, u_{n_2} \}` respectively, and let :math:`\bar{V} = \{ \bar{v}_1, \ldots, \bar{v}_{n_1+n_2} \}` be the union of the two vertex sets.
Let also :math:`\mathbf{K} \in \mathbb{R}^{(n_1+n_2) \times (n_1+n_2)}` be the kernel matrix defined as

.. math::

    \mathbf{K}_{ij} = \kappa(\bar{v}_i, \bar{v}_j) = \phi(\bar{v}_i)^\top \phi(\bar{v}_j)

Let :math:`\mathbf{u}_1, \ldots, \mathbf{u}_p` be a maximal orthonormal set of the non-zero eigenvalue eigenvectors of :math:`\mathbf{K}`
with corresponding eigenvalues :math:`\lambda_1, \ldots, \lambda_p`.
Then the vectors

.. math::

    \xi_i = \frac{1}{\sqrt{\lambda_i}} \sum_{l=1}^{n_1+n_2} [\mathbf{u}_i]_l \phi(\bar{v}_l)

where :math:`[\mathbf{u}_i]_l` is the :math:`l^{th}` component of vector :math:`\mathbf{u}_i` form an orthonormal basis for the subspace :math:`\{ \phi(\bar{v}_1), \ldots, \phi(\bar{v}_{n_1+n_2}) \}`.
Moreover, let :math:`\mathbf{Q} = [ \lambda_1^{1/2} \mathbf{u}_1, \ldots,\lambda_p^{1/2} \mathbf{u}_p ] \in \mathbb{R}^{p \times p}` and :math:`\mathbf{Q}_1, \mathbf{Q}_2` denote the first :math:`n_1` and last :math:`n_2 ` rows of matrix :math:`\mathbf{Q}` respectively.
Then, the generalized feature space Laplacian graph kernel induced from the base kernel :math:`\kappa` is defined as

.. math::

    k_{FLG}^\kappa(G_1, G_2) = \frac{| (\frac{1}{2} \mathbf{S}_1^{-1} + \frac{1}{2} \mathbf{S}_2^{-1} )^{-1} |^{1/2}}{|\mathbf{S}_1|^{1/4} |\mathbf{S}_2|^{1/4}} 

where :math:`\mathbf{S}_1 = \mathbf{Q}_1 \mathbf{L}_1^{-1} \mathbf{Q}_1^\top + \gamma \mathbf{I}` and :math:`\mathbf{S}_2 = \mathbf{Q}_2 \mathbf{L}_2^{-1} \mathbf{Q}_2^\top + \gamma \mathbf{I}` where :math:`\mathbf{I}` is the :math:`p \times p` identity matrix.

The multiscale Laplacian graph kernel builds a hierarchy of nested subgraphs, where each subgraph is centered around a vertex and computes the generalized feature space Laplacian graph kernel between every pair of these subgraphs.
Let :math:`G` be a graph with vertex set :math:`V`, and :math:`\kappa` a positive semi-definite kernel on :math:`V`.
Assume that for each :math:`v \in V`, we have a nested sequence of :math:`L` neighborhoods

.. math::

    v \in N_1(v) \subseteq N_2(v) \subseteq \ldots \subseteq N_L(v)

and for each :math:`N_l(v)`, let :math:`G_l(v)` be the corresponding induced subgraph of :math:`G`.
The multiscale Laplacian subgraph kernels are defined as :math:`\mathfrak{K}_1, \ldots, \mathfrak{K}_L : V \times V \rightarrow \mathbb{R}` as follows

1. :math:`\mathfrak{K}_1` is just the generalized feature space Laplacian graph kernel :math:`k_{FLG}^\kappa` induced from the base kernel :math:`\kappa` between the lowest level subgraphs (\ie the vertices)
    
.. math::

        \mathfrak{K}_1(v,u) = k_{FLG}^\kappa(v, u)
    
2. For :math:`l=2,3,\ldots,L`, :math:`\mathfrak{K}_l` is the the generalized feature space Laplacian graph kernel induced from :math:`\mathfrak{K}_{l-1}` between :math:`G_l(v)` and :math:`G_l(u)`
    
.. math::

        \mathfrak{K}_l(v,u) = k_{FLG}^{\mathfrak{K}_{l-1}}(G_l(v), G_l(u))
    

Then, the multiscale Laplacian graph Kernel between two graphs :math:`G_1, G_2` is defined as follows

.. math::

    k_{MLG}(G_1, G_2) = k_{FLG}^{\mathfrak{K}_L}(G_1, G_2)

The multiscale Laplacian graph kernel computes :math:`\mathfrak{K}_1` for all pairs of vertices, then computes :math:`\mathfrak{K}_2` for all pairs of vertices, and so on.
Hence, it requires :math:`\mathcal{O}(Ln^2)` kernel evaluations.
At the top levels of the hierarchy each subgraph centered around a vertex :math:`G_l(v)` may have as many as :math:`n` vertices.
Therefore, the cost of a single evaluation of the generalized feature space Laplacian graph kernel may take :math:`\mathcal{O}(n^3)` time.
This means that in the worst case, the overall cost of computing :math:`k_{MLG}` is :math:`\mathcal{O}(Ln^5)`.
Given a dataset of :math:`N` graphs, computing the kernel matrix requires repeating this for all pairs of graphs, which takes :math:`\mathcal{O}(LN^2n^5)` time and is clearly problematic for real-world settings.

The solution to this issue is to compute for each level :math:`l=1,2,\ldots,L+1` a single joint basis for all subgraphs at the given level across all graphs.
Let :math:`G_1, G_2, \ldots, G_N` be a collection of graphs, :math:`V_1, V_2, \ldots, V_N` their vertex sets, and assume that :math:`V_1, V_2, \ldots, V_N \subseteq \mathcal{V}` for some general vertex space :math:`\mathcal{V}`.
The joint vertex feature space of the whole graph collection is :math:`W = span \big\{ \bigcup_{i=1}^N \bigcup_{v \in V_i} \{ \phi(v) \} \big\}`.
Let :math:`c = \sum_{i=1}^N |V_i|` be the total number of vertices and :math:`\bar{V} = (\bar{v}_1, \ldots, \bar{v}_c)` be the concatenation of the vertex sets of all graphs.
Let :math:`\mathbf{K}` be the corresponding joint kernel matrix and :math:`\mathbf{u}_1, \ldots, \mathbf{u}_p` be a maximal orthonormal set of non-zero eigenvalue eigenvectors of :math:`\mathbf{K}` with corresponding eigenvalues :math:`\lambda_1,\ldots,\lambda_p` and :math:`p=dim(W)`.
Then the vectors

.. math::

    \xi_i = \frac{1}{\sqrt{\lambda_i}} \sum_{l=1}^c [\mathbf{u}_i]_l \phi(\bar{v}_l) \qquad i=1,\ldots,p

form an orthonormal basis for :math:`W`.
Moreover, let :math:`\mathbf{Q} = [ \lambda_1^{1/2} \mathbf{u}_1, \ldots, \lambda_p^{1/2} \mathbf{u}_p ] \in \mathbb{R}^{p \times p}` and :math:`\mathbf{Q}_1` denote the first :math:`n_1` rows of matrix :math:`\mathbf{Q}`, :math:`\mathbf{Q}_2` denote the next :math:`n_2 ` rows of matrix :math:`\mathbf{Q}` and so on.
For any pair of graphs :math:`G_i, G_j` of the collection, the generalized feature space Laplacian graph kernel induced from :math:`\kappa` can be expressed as

.. math::

    k_{FLG}^\kappa(G_i, G_j) = \frac{| (\frac{1}{2} \bar{\mathbf{S}}_i^{-1} + \frac{1}{2} \bar{\mathbf{S}}_j^{-1} )^{-1} |^{1/2}}{|\bar{\mathbf{S}}_i|^{1/4} |\bar{\mathbf{S}}_j|^{1/4}} 

where :math:`\bar{\mathbf{S}}_i = \mathbf{Q}_i \mathbf{L}_i^{-1} \mathbf{Q}_i^\top + \gamma \mathbf{I}`, :math:`\bar{\mathbf{S}}_j = \mathbf{Q}_j \mathbf{L}_j^{-1} \mathbf{Q}_j^\top + \gamma \mathbf{I}` and :math:`\mathbf{I}` is the :math:`p \times p` identity matrix.

The implementation of the multiscale Laplacian kernel can be found below

.. currentmodule:: grakel

.. autosummary::

   MultiscaleLaplacian


Low Rank Approximation
----------------------

Computing the kernel matrix between all vertices of all graphs (:math:`c` vertices in total) and storing it is a very costly procedure.
Computing its eigendecomposition is even worse in terms of the required runtime.
Morever, :math:`p` is also very large.
Managing the :math:`\bar{\mathbf{S}}_1, \ldots, \bar{\mathbf{S}}_N` matrices (each of which is of size :math:`p \times p`) becomes infeasible.
Hence, the multiscale Laplacian graph kernel replaces :math:`W` with a smaller, approximate joint features space.
Let :math:`\tilde{V} = (\tilde{v}_1, \ldots, \tilde{v}_{\tilde{c}})` be :math:`\tilde{c} \ll c` vertices sampled from the joint vertex set.
Then, the corresponding subsampled vertex feature space is :math:`\tilde{W} = span \{ \phi(v) : v \in \tilde{V} \}`.
Let :math:`\tilde{p} = dim(\tilde{W})`.
Similarly to before, the kernel constructs an orthonormal basis :math:`\{ \xi_1, \ldots, \xi_{\tilde{p}} \}` for :math:`\tilde{W}` by forming the (now much smaller) kernel matrix :math:`\mathbf{K}_{ij} = \kappa(\tilde{v}_i, \tilde{v}_j)`, computing its eigenvalues and eigenvectors, and setting :math:`\xi_i = \frac{1}{\sqrt{\lambda_i}} \sum_{l=1}^{\tilde{c}} [\mathbf{u}_i]_l \phi(\tilde{v}_l)`. 
The resulting approximate generalized feature space Laplacian graph kernel is

.. math::

    k_{FLG}^\kappa(G_1, G_2) = \frac{| (\frac{1}{2} \tilde{\mathbf{S}}_1^{-1} + \frac{1}{2} \tilde{\mathbf{S}}_2^{-1} )^{-1} |^{1/2}}{|\tilde{\mathbf{S}}_1|^{1/4} |\tilde{\mathbf{S}}_2|^{1/4}} 

where :math:`\tilde{\mathbf{S}}_1 = \tilde{\mathbf{Q}}_1 \mathbf{L}_1^{-1} \tilde{\mathbf{Q}}_1^\top + \gamma \mathbf{I}`, :math:`\tilde{\mathbf{S}}_2 = \tilde{\mathbf{Q}}_2 \mathbf{L}_2^{-1} \tilde{\mathbf{Q}}_2^\top + \gamma \mathbf{I}` are the projections of :math:`\bar{\mathbf{S}}_1` and :math:`\bar{\mathbf{S}}_2` to :math:`\tilde{W}` and :math:`\mathbf{I}` is the :math:`\tilde{p} \times \tilde{p}` identity matrix. Finally, the kernel introduces a further layer of approximation by restricting :math:`\tilde{W}` to be the space spanned by the first :math:`\hat{p} < \tilde{p}` basis vectors (ordered by descending eigenvalue), effectively doing kernel PCA on :math:`\{ \phi(\tilde{v}) \}_{\tilde{v} \in \tilde{V}}`.
The combination of these two factors makes computing the entire stack of kernels feasible, reducing the complexity of computing the kernel matrix for a dataset of :math:`N` graphs to :math:`\mathcal{O}(NL \tilde{c}^2 \hat{p}^3 + NL \tilde{c}^3 + N^2 \hat{p}^3)`.

The approximate multiscale Laplacian graph kernel can be found below

.. currentmodule:: grakel

.. autosummary::

   MultiscaleLaplacianFast

Bibliography
------------
.. bibliography:: graph_kernels.bib
   :filter: docname in docnames
