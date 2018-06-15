.. _edge_histogram:

Edge Histogram Kernel
=====================

The edge histogram kernel is a basic linear kernel on edge label histograms.
The kernel assumes edge-labeled graphs.
Let :math:`\mathcal{G}` be a collection of graphs, and assume that each of their edges comes from an abstract edge space :math:`\mathcal{E}`.
Given a set of node labels :math:`\mathcal{L}`, :math:`\ell : \mathcal{E} \rightarrow \mathcal{L}` is a function that assigns labels to the edges of the graphs.
Assume that there are :math:`d` labels in total, that is :math:`d = |\mathcal{L}|`.
Then, the edge label histogram of a graph :math:`G=(V,E)` is a vector :math:`\mathbf{f} = (f_1, f_2, \ldots, f_d)`, such that :math:`f_i = |\{ (v,u) \in E : \ell(v,u) = i \}|` for each :math:`i \in \mathcal{L}`.
Let :math:`\mathbf{f}, \mathbf{f}'` be the edge label histograms of two graphs :math:`G, G'`, respectively.
The edge histogram kernel is then defined as the linear kernel between :math:`\mathbf{f}` and :math:`\mathbf{f}'`, that is

.. math::

    k(G, G') = \langle \mathbf{f}, \mathbf{f}' \rangle

The complexity of the edge histogram kernel is linear in the number of edges of the graphs.

An implementation of that kernel can be found below

.. currentmodule:: grakel

.. autosummary::

   EdgeHistogram
