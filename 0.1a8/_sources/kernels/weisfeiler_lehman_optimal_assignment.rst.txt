.. _weisfeiler_lehman:

Weisfeiler-Lehman Optimal Assignment
====================================

The Weisfeiler-Lehman optimal assignment kernel capitalizes on the theory of valid assignment kernels to improve the performance of the Weisfeiler-Lehman subtree kernel :cite:`kriege2016valid`.
Before we delve into the details of the kernel, it is necessary to introduce the theory of valid optimal assignment kernels.

Definition: Valid Assignment Kernels
------------------------------------

Let :math:`\mathcal{X}` be a set, and :math:`[\mathcal{X}]^n` denote the set of all :math:`n`-element subsets of :math:`\mathcal{X}`.
Let also :math:`X,X' \in [\mathcal{X}]^n` for :math:`n \in \mathbb{N}`, and :math:`\mathfrak{B}(X,X')` denote the set of all bijections between :math:`X` and :math:`X'`.
The optimal assignment kernel on :math:`[\mathcal{X}]^n` is defined as

.. math::
    :nowrap:

    \begin{equation}
      K_\mathfrak{B}^k(X,X') = \max_{B \in \mathfrak{B}(X,X')} \sum_{(x,x') \in B} k(x,x')
    \end{equation}

where :math:`k` is a kernel between the elements of :math:`X` and :math:`X'`.
Kriege et al. showed that the above function :math:`K_\mathfrak{B}(\mathcal{X},\mathcal{X}')` is a valid kernel only if the base kernel :math:`k` is strong.
Strong kernels are equivalent to kernels obtained from a hierarchy defined on set :math:`\mathcal{X}`.
More specifically, let :math:`T` be a rooted tree such that the leaves of :math:`T` are the elements of :math:`\mathcal{X}`.
Let :math:`V(T)` be the set of vertices of :math:`T`.
Each inner vertex :math:`v` in :math:`T` corresponds to a subset of :math:`\mathcal{X}` comprising all leaves of the subtree rooted at :math:`v`.
Let :math:`w : V(T) \rightarrow \mathbb{R}_{\geq 0}` be a weight function such that :math:`w(v) \geq w(p(v))` for all :math:`v` in :math:`T` where :math:`p(v)` is the parent of vertex :math:`v`.
Then, the tuple :math:`(T,w)` defines a hierarchy.
Let :math:`LCA(u,v)` be the lowest common ancestor of vertices :math:`u` and :math:`v`, that is, the unique vertex with maximum depth that is an ancestor of both :math:`u` and :math:`v`.
Let :math:`H = (T,w)` be a hierarchy on :math:`\mathcal{X}`, then the function defined as :math:`k(x,y) = w(LCA(x,y))` for all :math:`x,y` in :math:`\mathcal{X}` is the kernel on :math:`\mathcal{X}` induced by :math:`H`.
Interestingly, strong kernels are equivalent to kernels obtained from a hierarchical partition of the domain of the kernel. 
Hence, by constructing a hierarchy on :math:`\mathcal{X}`, we can derive a strong kernel :math:`k` and ensure that the emerging assignment function is a valid kernel.

Based on the property that every strong kernel is induced by a hierarchy, we can derive explicit feature maps for strong kernels.
Let :math:`\omega : V(T) \rightarrow \mathbb{R}_{\geq 0}` be an additive weight function defined as :math:`\omega(v) = w(v) - w(p(v))` and :math:`\omega(r) = w(r)` for the root :math:`r`.
Note that the property of a hierarchy assures that the values of the :math:`\omega` function are nonnegative.
For :math:`v \in V(T)`, let :math:`P(v) \subseteq V(T)` denote the vertices on the path from :math:`v` to the root :math:`r`.
The strong kernel :math:`k` induced by the hierarchy :math:`H` can be defined using the mapping :math:`\phi : \mathcal{X} \rightarrow \mathbb{R}^n`, where :math:`n = |V(T)|` and the components indexed by :math:`v \in V(T)` are

.. math::
    :nowrap:

    \begin{equation}
      \phi(v) = \left\{
      \begin{array}{lr}
        \sqrt{\omega(u)} & \text{if }u \in P(v),\\
        \quad 0 & \text{otherwise}
      \end{array}
    \right.
    \end{equation}

The Figure below shows an example of a strong kernel, an associated hierarchy and the derived feature vectors.

.. figure:: ../_figures/optimal_assignment_example.png

    The matrix of a strong kernel on objects :math:`a,b` and :math:`c` (a) induced by the hierarchy (b) and the derived feature vectors (c). A vertex :math:`v` in (b) is annotated by its weights :math:`w(v);\omega(v)`.

Let :math:`H = (T,w)` be a hierarchy on :math:`\mathcal{X}`.
As mentioned above, the hierarchy :math:`H` induces a strong kernel :math:`k`.
Since :math:`k` is strong, the function :math:`K_\mathfrak{B}^k` is a valid kernel.
The kernel :math:`K_\mathfrak{B}^k` can be computed in linear time in the number of vertices :math:`n` of the tree :math:`T` using the histogram intersection kernel as follows

.. math::
    :nowrap:

    \begin{equation}
        K_\mathfrak{B}^k(X, X') = \sum_{i=1}^n \min\big(H_{X}(i),H_{X'}(i)\big)
    \end{equation}

which is known to be a valid kernel on :math:`\mathbb{R}^n` :cite:`barla2003histogram`.
Hence, the complexity of the proposed kernel depends on the size of the tree :math:`T`.
The following Figure illustrates the relation between the optimal assignment kernel employing a strong base kernel and the histogram intersection kernel.

.. figure:: ../_figures/optimal_assignment_histograms.png

    An assignment instance (a) for :math:`X,Y \in [\mathcal{X}]^5` and the derived histograms (b). The set :math:`X` contains three distinct vertices labeled :math:`a` and the set :math:`Y` two distinct vertices labeled :math:`b` and :math:`c`. Taking the multiplicities into account the histograms are obtained from the hierarchy of the base kernel :math:`k` depicted in the previous Figure.  The optimal assignment yields a value of :math:`K_\mathfrak{B}^k(X, Y)` :math:`= \sum_{i=1}^n \min\big(H_{X}(i),H_{Y}(i)\big)` :math:`= \min\{5,5\}` :math:`+ \min\{8,6\}` :math:`+ \min\{3,1\}` :math:`+ \min\{2,4\}` :math:`+ \min\{1s,2\}=15`.


Definition: Weisfeiler-Lehman Optimal Assignment
------------------------------------------------

We next present the Weisfeiler-Lehman optimal assignment kernel.
Let :math:`G=(V,E)` and :math:`G'=(V',E')` be two graphs.
The Weisfeiler-Lehman optimal assignment kernel is defined as

.. math::
    :nowrap:

    \begin{equation}
        k(G,G') = K_\mathfrak{B}^k(V,V')
    \end{equation}

where :math:`k` is the following base kernel

.. math::
    :nowrap:
    
    \begin{equation}
        k(v,v') = \sum_{i=0}^h \delta(\tau_i(v), \tau_i(v'))
    \end{equation}

where :math:`\tau_i(v)` is the label of node :math:`v` at the end of the :math:`i^{th}` iteration of the Weisfeiler-Lehman relabeling procedure.

The base kernel value reflects to what extent two vertices :math:`v` and :math:`v'` have a similar neighborhood.
It can be shown that the colour refinement process of the Weisfeiler-Lehman algorithm defines a hierarchy on the set of all vertices of the input graphs.
Specifically, the sequence :math:`(\tau_i)_{0\leq i \leq h}` gives rise to a family of nested subsets, which can naturally be represented by a hierarchy :math:`(T,w)`.
When assuming :math:`\omega(v) = 1` for all vertices :math:`v \in V(T)`, the hierarchy induces the kernel defined above.
Such a hierarchy for a graph on six vertices is illustrated in the following Figure.

.. figure:: ../_figures/wl_optimal_assignment.png

    A graph :math:`G` with uniform initial labels :math:`\tau_0` and refined labels :math:`\tau_i` for :math:`i \in \{1,\ldots,3\}` (a), and the associated hierarchy (b).


The above kernel is implemented below

.. currentmodule:: grakel

.. autosummary::

   WeisfeilerLehmanOptimalAssignment

Bibliography
------------
.. bibliography:: graph_kernels.bib
   :filter: docname in docnames
