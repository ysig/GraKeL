.. _vertex_histogram:

Vertex Histogram Kernel
=======================

The vertex histogram kernel is a basic linear kernel on vertex label histograms.
The kernel assumes node-labeled graphs.
Let :math:`\mathcal{G}` be a collection of graphs, and assume that each of their vertices comes from an abstract vertex space :math:`\mathcal{V}`.
Given a set of node labels :math:`\mathcal{L}`, :math:`\ell : \mathcal{V} \rightarrow \mathcal{L}` is a function that assigns labels to the vertices of the graphs.
Assume that there are :math:`d` labels in total, that is :math:`d = |\mathcal{L}|`.
Then, the vertex label histogram of a graph :math:`G=(V,E)` is a vector :math:`\mathbf{f} = (f_1, f_2, \ldots, f_d)`, such that :math:`f_i = |\{ v \in V : \ell(v) = i \}|` for each :math:`i \in \mathcal{L}`.
Let :math:`\mathbf{f}, \mathbf{f}'` be the vertex label histograms of two graphs :math:`G, G'`, respectively.
The vertex histogram kernel is then defined as the linear kernel between :math:`\mathbf{f}` and :math:`\mathbf{f}'`, that is

.. math::

    k(G, G') = \langle \mathbf{f}, \mathbf{f}' \rangle

The complexity of the vertex histogram kernel is linear in the number of vertices of the graphs.

An implementation of that kernel can be found below

.. currentmodule:: grakel

.. autosummary::

   VertexHistogram
