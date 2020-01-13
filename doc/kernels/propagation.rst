.. _propagation:

Propagation Kernel
======================
Propagation kernels where introduced as a general framework in :cite:`neumann2015propagation`. They are based in the idea of propagating label information between nodes of the graph, based on the graph structure.
A graph is considered to have **attributes** on nodes, where in the case of labels they correspond to One-Hot-Vectors of the full dictionary of labels.
The totality of nodes for each graph, can be seen as a probability distribution :math:`P` of size :math:`n \times d` where :math:`n` corresponds to the number of nodes and :math:`d` to the size of attributes.
After the idea of diffusion is applied in order to construct the algorithmic framework of propagation kernels.
Given a transition matrix :math:`T` that is row normalized, an iterative propagation scheme will be build on the basis of the following simple substitution rule:

.. math::

    P_{t+1} \leftarrow T P_{t}

Other than a user given transition matrix, :math:`T = D^{−1}A`, was considered as default for each graph, where :math:`D = diag(\sum_{j} A_{ij})` and :math:`A` corresponds to the adjacency matrix of this graph.
The general algorithm for propagation kernels is as follows:

.. figure:: ../_static/marion1.png

    The general algorithmic scheme for propagation kernels.


The kernel computation :math:`\langle \Phi, \Phi \rangle_{ij}`, at iteration :math:`t` between two graphs :math:`i, j` is equivalent with the:

.. math::

    K(G^{(i)}_{t}, G^{(j)}_{t}) = \sum_{u \in G^{(i)}_{t}} \sum_{v \in G^{(j)}_{t}} k(u, v)

where the node kernel :math:`k(u, v)` is resolved through binning.
In order to bin nodes a method should be found that was both efficient and expressive.
A simple hashing function was considered a bad idea as it would separate values that where much more common than others. A sense of *locality* was needed when binning in order to group similar diffusion patterns in the same bin, similar to what is shown in the following:

.. figure:: ../_static/marion2.png

    A binning scheme between a two step label propagation.

For that the technique of locally sensitive hashing [**LSH**] was utilized, which was applied to all the input graphs as shown in the following:

.. figure:: ../_static/marion3.png

    The locally sensitive function implemented inside the kernel.

Finally the following algorithm was implemented, in our case where we consider all graphs to be *fully-labeled*:

.. figure:: ../_static/marion4.png

    The propagation kernel algorithm, which was implemented inside the package.

In case we have an attributed graph the :math:`P_{0} \leftarrow \delta_{l(V)}` is replaced by :math:`P_{0} \leftarrow attr(V)` considering all the node attributes have the same dimension.

Both for attributed and for labeled graphs, implementations of the propagation can be found below:

.. currentmodule:: grakel

.. autosummary::

   Propagation
   PropagationAttr

Bibliography
------------
.. bibliography:: graph_kernels.bib
   :filter: docname in docnames