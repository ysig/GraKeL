.. _api_ref:

=============
API Reference
=============

This is the class and function reference of **grakel**.
In order for the user to understand how to use the **grakel** project and how it works we suggest him
to read :ref:`user_manual` section.

:mod:`grakel.graph`: Graph class with its utility functions
===========================================================

Base Class
----------
.. currentmodule:: grakel

.. autosummary::
   :toctree: generated/
   :template: class.rst

   Graph


Utility Functions
-----------------
.. currentmodule:: grakel

.. autosummary::
   :toctree: generated/
   :template: function.rst

   graph.is_adjacency
   graph.is_edge_dictionary
   graph.laplacian
   graph.floyd_warshall

**User guide:** See the :ref:`graph` section for further details.

:mod:`grakel.graph_kernels`: A kernel decorator
===============================================
.. currentmodule:: grakel

Graph Kernel (decorator)
------------------------
.. autosummary::
   :toctree: generated/
   :template: class.rst

   grakel.GraphKernel

**User guide:** See the :ref:`graph_kernel` section for further details.

:mod:`grakel.kernels`: A collection of graph kernels
====================================================

Kernels
-------

.. currentmodule:: grakel

.. autosummary::
   :toctree: generated/
   :template: kernel.rst

   kernel
   random_walk
   random_walk_labeled
   pyramid_match
   neighborhood_hash
   shortest_path
   shortest_path_attr
   subgraph_matching
   weisfeiler_lehman
   hadamard_code
   neighborhood_subgraph_pairwise_distance
   lovasz_theta
   svm_theta
   propagation
   odd_sth
   multiscale_laplacian
   multiscale_laplacian_fast
   hadamard_code
   vertex_histogram
   edge_histogram

**User guide:** See the :ref:`kernels` section for further details.

:mod:`grakel.datasets`: Datasets
=================================

Fetch
-----

.. currentmodule:: grakel.datasets

.. autosummary::
   :toctree: generated/
   :template: function.rst
   
   fetch_dataset
   get_dataset_info

**User guide:** See the :ref:`datasets` section for further details.

.. _gd:	https://ls11-www.cs.tu-dortmund.de/staff/morris/graphkerneldatasets