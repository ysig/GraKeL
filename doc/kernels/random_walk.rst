.. _random_walk:

.. raw:: latex

    \newtheorem{definition}{Definition}

Random Walk Kernel
==================
The most well-studied family of graph kernels is probably the *random walk kernels* which quantify the similarity between a pair of graphs based on the number of common walks in the two graphs 
:cite:`kashima2003marginalized`, :cite:`gartner2003graph`, :cite:`mahe2004extensions`, :cite:`borgwardt2005protein`, :cite:`vishwanathan2010graph`, :cite:`sugiyama2015halting`.

Kernels belonging to this family have concentrated mainly on counting matching walks in the two input graphs. There are several variations of random walk kernels. The :math:`k`-step random walk kernel compares random walks up to length :math:`k` in the two graphs. The most widely-used kernel from this family is the geometric random walk kernel :cite:`gartner2003graph` which compares walks up to infinity assigning a weight :math:`\lambda^k` (:math:`\lambda < 1`) to walks of length :math:`k` in order to ensure convergence of the corresponding geometric series. We next give the formal definition of the geometric random walk kernel.
Given two node-labeled graphs :math:`G_i=(V_i,E_i)` and :math:`G_j=(V_j,E_j)`, their direct product
:math:`G_\times=(V_\times,E_\times)` is a graph with vertex set:

.. math::
    :nowrap:   

    \begin{equation}
	    V_{\times} = \{(v_i,v_j) : v_i \in V_i \wedge v_j \in V_j \wedge \ell(v_i) = \ell(v_j) \} 
    \end{equation}

and edge set:

.. math::
    :nowrap:

    \begin{equation}
	    E_{\times} = \{\{(v_i,v_j),(u_i,u_j)\} : \{v_i,u_i\} \in E_i \wedge \{v_j,u_j\} \in E_j\}
    \end{equation}

Performing a random walk on :math:`G_{\times}` is equivalent to performing a simultaneous random walk
on :math:`G_i` and :math:`G_j`.
The geometric random walk kernel counts common walks (of potentially infinite length) in two graphs
and is defined as follows.

Definition: Geometric Random Walk Kernel
----------------------------------------

Let :math:`G_i` and :math:`G_j` be two graphs, let :math:`A_\times` denote the adjacency matrix of their
product graph :math:`G_\times`, and let :math:`V_\times` denote the vertex set of the product
graph :math:`G_\times`.

Then, the geometric random walk kernel is defined as

.. math::
   :nowrap:

   \begin{equation}
       	K_{\times}^{\infty}(G_i,G_j) = \sum_{p,q=1}^{|V_{\times}|} \Big[ \sum_{l=0}^{\infty} \lambda^l A_{\times}^l \Big]_{pq} = e^T(I - \lambda A_{\times})^{-1} e
   \end{equation}

where :math:`I` is the identity matrix, :math:`e` is the all-ones vector, and :math:`\lambda`
is a positive, real-valued weight. The geometric random walk kernel converges only if
:math:`\lambda < \frac{1}{\lambda_\times}` where :math:`\lambda_\times` is the largest eigenvalue of
:math:`A_{\times}`.

Direct computation of the geometric random walk kernel requires :math:`\mathcal{O}(n^6)` time.
The computational complexity of the method severely limits its applicability to real-world applications.
To account for this, Vishwanathan et al. proposed in :cite:`vishwanathan2010graph` four efficient 
methods to compute random walk graph kernels which generally reduce the computational complexity from 
:math:`\mathcal{O}(n^6)` to :math:`\mathcal{O}(n^3)`.
Mahé et al. proposed in :cite:`mahe2004extensions` some other extensions of random walk kernels.
Specifically, they proposed a label enrichment approach which increases specificity and in most
cases also reduces computational complexity.
They also employed a second order Markov random walk to deal with the problem of "tottering".
Sugiyama and Borgwardt focused in :cite:`sugiyama2015halting` on a different problem of random walk
kernels, a phenomenon referred to as "halting".

Next follow two implementations of this kernel (one for unlabeled graphs and one for graphs with discrete node labels)

.. currentmodule:: grakel

.. autosummary::
   RandomWalk
   RandomWalkLabeled


Bibliography
------------
.. bibliography:: graph_kernels.bib
   :filter: docname in docnames
