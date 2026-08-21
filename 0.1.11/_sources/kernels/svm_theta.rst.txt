.. _svm_theta:

SVM Theta Kernel
================

The SVM-:math:`\vartheta` kernel is very related to the Lovász :math:`\vartheta` kernel :cite:`johansson2014global`.
The Lovász :math:`\vartheta` kernel suffers from high computational complexity, and the SVM-:math:`\vartheta` kernel was developed as a more efficient alternative. 
Similar to the Lovász :math:`\vartheta` kernel, this kernel also assumes unlabeled graphs.

Given a graph :math:`G=(V,E)` such that :math:`|V| = n`, the Lovász number of :math:`G` can be defined as

.. math::

    \vartheta(G) = \min_{\mathbf{K} \in L} \omega(\mathbf{K})

where :math:`\omega(\mathbf{K})` is the one-class SVM given by

.. math:: \omega(\mathbf{K}) = \max_{\alpha_i > 0} 2\sum_{i=1}^{n} \alpha_i - \sum_{i=1}^{n} \sum_{j=1}^{n} \alpha_i \alpha_j \mathbf{K}_{ij}
   :label: oneclass_svm


and :math:`L` is a set of positive semidefinite matrices defined as

.. math::

    L = \{ \mathbf{K} \in S_{n}^+ : \mathbf{K}_{ii} = 1, \mathbf{K}_{ij}=0 \: \forall (i,j) \not \in E \}

where :math:`S_{n}^+` is the set of all :math:`n \times n` positive semidefinite matrices.

The SVM-:math:`\vartheta` kernel first computes the matrix :math:`\mathbf{K}_{LS}` which is equal to

.. math::

    \mathbf{K}_{LS} = \frac{\mathbf{A}}{\rho} + \mathbf{I}

where :math:`\mathbf{A}` is the adjacency matrix of :math:`G`, :math:`\mathbf{I}` is the :math:`n \times n` identity matrix, and :math:`\rho \geq -\lambda_n` with :math:`\lambda_n` the minimum eigenvalue of :math:`\mathbf{A}`.
The matrix :math:`\mathbf{K}_{LS}` is positive semidefinite by construction and it has been shown in :cite:`jethava2013lovasz` that

.. math::

    \omega(\mathbf{K}_{LS}) = \sum_{i=1}^n \alpha_i

where :math:`\alpha_i` are the maximizers of Equation :eq:`oneclass_svm`. 
Furthermore, it was shown that on certain families of graphs (e.g. Erdős–Rényi random graphs), :math:`\omega(\mathbf{K}_{LS})` is with high probability a constant factor approximation to :math:`\vartheta(G)`.

Then, the SVM-:math:`\vartheta` kernel is defined as follows

.. math::

    k_{SVM}(G, G') = \sum_{S \subseteq V} \sum_{S' \subseteq V'} \delta(|S|, |S'|) \frac{1}{Z_{|S|}} k \Big(\sum_{i \in S} \alpha_i, \sum_{j \in S'} \alpha_j \Big)

where :math:`Z_{|S|} = \binom{|V|}{|S|} \binom{|V'|}{|S|}`, :math:`\delta(|S|, |S'|)` is a delta kernel (equal to :math:`1` if :math:`|S|=|S'|`, and :math:`0` otherwise), and :math:`k` is a positive semi-definite kernel between real values (\eg linear kernel, gaussian kernel).

The SVM-:math:`\vartheta` kernel consists of three main steps: (:math:`1`) constructing matrix :math:`\mathbf{K}_{LS}` of :math:`G` which takes :math:`\mathcal{O}(n^3)` time (:math:`2`) solving the one-class SVM problem in :math:`\mathcal{O}(n^2)` time to obtain the :math:`\alpha_i` values, and (:math:`3`) computing the sum of the :math:`\alpha_i` values for all subgraphs (\ie subsets of vertices :math:`S \subseteq V`) of each graph.
Computing the above quantity for all :math:`2^n` sets of vertices is not feasible in real-world scenarios.

To address the above issue, the SVM-:math:`\vartheta` kernel employs sampling schemes.
Given a graph :math:`G`, the kernel samples a specific number of subgraphs :math:`\mathfrak{S} \in 2^V`.
Then, the SVM-:math:`\vartheta` kernel is defined as follows

.. math::

    \hat{k}_{SVM}(G, G') = \sum_{S \subseteq \mathfrak{S}} \sum_{S' \subseteq \mathfrak{S}'} \delta(|S|, |S'|) \frac{1}{\hat{Z}_{|S|}} \Big(\sum_{i \in S} \alpha_i, \sum_{j \in S'} \alpha_j \Big)

where :math:`\hat{Z}_{|S|} = |\mathfrak{S}_{|S|}| |\mathfrak{S}'_{|S|}|` and :math:`\mathfrak{S}_{|S|}` denotes the subset of :math:`\mathfrak{S}` consisting of all sets of cardinality :math:`|S|`, that is :math:`\mathfrak{S}_{|S|} = \{ B \in \mathfrak{S} : |B| = |S| \}`.

The time complexity of computing :math:`\hat{k}_{SVM}(G, G')` is :math:`\mathcal{O}(n^3 + s^2 T(k) + sn)` where :math:`T(k)` is the complexity of computing the base kernel :math:`k` and :math:`s = \max(|\mathfrak{S}|, |\mathfrak{S}'|)`.
The first term represents the cost of computing :math:`\mathbf{K}_{LS}` (dominated by the eigenvalue decomposition).
The second term corresponds to the worst-case complexity of comparing the sums of the :math:`\alpha_i` values.
And finally, the third term is the cost of computing the sum of the :math:`\alpha_i` values for the sampled subsets of vertices.

The implementation of the SVM-:math:`\vartheta` kernel can be found below

.. currentmodule:: grakel

.. autosummary::

   SvmTheta

Bibliography
------------
.. bibliography:: graph_kernels.bib
   :filter: docname in docnames
