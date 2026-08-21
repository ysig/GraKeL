.. _weisfeiler_lehman:

Weisfeiler Lehman Framework
===========================

This Weisfeiler Lehman framework operates on top of existing graph kernels and is inspired by the
Weisfeiler-Lehman test of graph isomorphism :cite:`weisfeiler1968reduction`.
The key idea of the Weisfeiler-Lehman algorithm is to replace the label of each vertex with a multiset
label consisting of the original label of the vertex and the sorted set of labels of its neighbors.
The resultant multiset is then compressed into a new, short label.
This relabeling procedure is then repeated for :math:`h` iterations.
Note that this procedure is performed simultaneously on all input graphs.
Therefore, two vertices from different graphs will get identical new labels
if and only if they have identical multiset labels.

More formally, given a graph :math:`G=(V,E)` endowed with a labeling function :math:`\ell=\ell_0`,
the Weisfeiler-Lehman graph of :math:`G` at height :math:`i` is a graph :math:`G_i=(V,E)` endowed
with a labeling function :math:`\ell_i` which has emerged after :math:`i` iterations of the
relabeling procedure described above.

The Weisfeiler-Lehman sequence up to height :math:`h` of :math:`G` consists of the Weisfeiler-Lehman
graphs of :math:`G` at heights from :math:`0` to :math:`h`, :math:`\{ G_0,G_1,\ldots,G_h\}`.


Definition: Weisfeiler-Lehman Framework
---------------------------------------

Let :math:`k` be any kernel for graphs, that we will call the base kernel.
Then the Weisfeiler-Lehman kernel with :math:`h` iterations with the base
kernel :math:`k` between two graphs :math:`G` and :math:`G'` is defined as

.. math::
    :nowrap:

	\begin{equation}
		k_{WL}(G,G') = k(G_0,G_0') + k(G_1,G_1') + \ldots + k(G_h,G_h')
	\end{equation}

where :math:`h` is the number of Weisfeiler-Lehman iterations, and 
:math:`\{ G_0,G_1,\ldots,G_h\}` and :math:`\{ G_0',G_1',\ldots,G_h'\}`
are the Weisfeiler-Lehman sequences of :math:`G` and :math:`G'` respectively.

From the above definition, it is clear that any graph kernel that takes into
account discrete node labels can take advantage of the Weisfeiler-Lehman framework
and compare graphs based on the whole Weisfeiler-Lehman sequence.

The general implementation of this framework can be found here:

.. currentmodule:: grakel

.. autosummary::

   WeisfeilerLehman

It should support all :code:`Kernel` objects inside its parameter :code:`base_kernel` (formulated in the correct way).


Weisfeiler-Lehman Subtree Kernel
--------------------------------
The Weisfeiler-Lehman subtree kernel is a very popular algorithm, and is considered the state-of-the-art in graph classification.
Let :math:`G`, :math:`G'` be two graphs. Define :math:`\Sigma_i \subseteq \Sigma` as the set of letters that occur as node labels
at least once in :math:`G` or :math:`G'` at the end of the :math:`i^{th}` iteration of the Weisfeiler-Lehman algorithm.
Let :math:`\Sigma_0` be the set of original node labels of :math:`G` and :math:`G'`.
Assume all :math:`\Sigma_i` are pairwise disjoint.
Without loss of generality, assume that every :math:`\Sigma_i = \{ \sigma_{i1},\ldots,\sigma_{i|\Sigma_i|} \}` is ordered.
Define a map :math:`c_i : \{ G,G' \} \times \Sigma_i \rightarrow \mathbb{N}` such that :math:`c_i(G, \sigma_{ij})`
is the number of occurrences of the letter :math:`\sigma_{ij}` in the graph :math:`G`.

The Weisfeiler-Lehman subtree kernel on two graphs :math:`G` and :math:`G'` with :math:`h` iterations is defined as

.. math::
    :nowrap:

    \begin{equation}
	    k(G,G') = \langle \phi(G),\phi(G') \rangle 
    \end{equation}

where

.. math::
    :nowrap:

    \begin{equation}
	    \phi(G) = (c_0(G,\sigma_{01}),\ldots,c_0(G,\sigma_{0|\Sigma_0|}),\ldots,c_h(G,\sigma_{h1}),\ldots,c_h(G,\sigma_{h|\Sigma_h|}))
    \end{equation}

and

.. math::
    :nowrap:

    \begin{equation}
	    \phi(G') = (c_0(G',\sigma_{01}),\ldots,c_0(G',\sigma_{0|\Sigma_0|}),\ldots,c_h(G',\sigma_{h1}),\ldots,c_h(G',\sigma_{h|\Sigma_h|}))
    \end{equation}

It can be shown that the above definition is equivalent to comparing the number of shared subtrees between the two input graphs :cite:`shervashidze2011weisfeiler`.
It is interesting to note that the Weisfeiler-Lehman subtree kernel exhibits an attractive computational complexity since it can be computed in :math:`\mathcal{O}(hm)` time.

.. note::
    
    To create an instance of the above kernel use the :ref:`vertex_histogram` as the :code:`base_kernel`.

Bibliography
------------
.. bibliography:: graph_kernels.bib
   :filter: docname in docnames
