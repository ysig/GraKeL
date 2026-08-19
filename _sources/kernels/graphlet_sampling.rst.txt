.. _graphlet_sampling:

Graphlet Sampling Kernel
========================

The graphlet sampling kernel decomposes graphs into graphlets (i.e. small subgraphs with :math:`k` nodes
where :math:`k \in \{ 3,4,5, \ldots \}`) :cite:`prvzulj2007biological` and counts matching graphlets
in the input graphs. Let :math:`\mathcal{G} = \{ graphlet_1,graphlet_2, \ldots, graphlet_r\}` be the set
of size-:math:`k` graphlets.
Let also :math:`f_G \in \mathbb{N}^r` be a vector such that its :math:`i`-th entry is equal to the
frequency of occurrence of :math:`graphlet_i` in :math:`G`, :math:`f_{G,i} = \#(graphlet_i \sqsubseteq G)`.
Then, the graphlet kernel is defined as follows.

Graphlet of size :math:`k` Kernel
---------------------------------
Let :math:`G_i`, :math:`G_j` be two graphs of size :math:`n \geq k`, and :math:`f_{G_i}, f_{G_j}` vectors
that count the occurrence of each graphlet of size :math:`k` in the two graphs. 
Then the graphlet kernel is defined as

.. math::
    :nowrap:

	\begin{equation}
  		k(G_i,G_j) = f_{G_i}^\top \ f_{G_j}
  	\end{equation}

As is evident from the above definition, the graphlet kernel is computed by explicit feature maps.
First, the representation of each graph in the feature space is computed.
And then, the kernel value is computed as the dot product of the two feature vectors.
The main problem of graphlet kernel is that an exaustive enumeration of graphlets is very expensive.
Since there are :math:`\binom{n}{k}` size-:math:`k` subgraphs in a graph, computing the feature vector
for a graph of size :math:`n` requires :math:`\mathcal{O}(n^k)` time.
To account for that, Shervashidze et al. resorted to sampling :cite:`shervashidze2009efficient`.
Following Weissman et al. :cite:`weissman2003inequalities`, they showed that by sampling a fixed number
of graphlets the empirical distribution of graphlets will be sufficiently close to their actual distribution
in the graph.

Below follows the implemented graphlet sampling kernel. By using the parameter *sampling* the user
can explore various possibilities from sampling all graphlets, to sampling probabilistically based
either on the number of samples or the satisfaction of a certain probabilistic on error.

.. currentmodule:: grakel

.. autosummary::

   GraphletSampling

Bibliography
------------
.. bibliography:: graph_kernels.bib
   :filter: docname in docnames
