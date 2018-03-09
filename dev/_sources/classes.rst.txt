.. _api_ref:

=============
API Reference
=============

This is the class and function reference of grakel.
.. Please refer to the :ref:`full user guide <user_guide>` for further details, as the class and
.. function raw specifications may not be enough to give full guidelines on their
.. uses.
.. For reference on concepts repeated across the API, see :ref:`glossary`.

:mod:`grakel.graph`: Graph class with its utility functions
===========================================================

Base Class
----------
.. currentmodule:: grakel

.. autosummary::
   :toctree: generated/
   :template: class.rst

   grakel.Graph


Utility Functions
-----------------
.. currentmodule:: grakel.graph

.. autosummary::
   :toctree: generated/
   :template: function.rst

   grakel.graph.is_adjacency
   grakel.graph.is_edge_dictionary
   grakel.graph.laplacian
   grakel.graph.floyd_warshall

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

   grakel.kernel

.. currentmodule:: grakel.kernels

.. autosummary::
   :toctree: generated/
   :template: kernel.rst

   grakel.kernels.kernel
   grakel.kernels.random_walk
   grakel.kernels.pyramid_match
   grakel.kernels.neighborhood_hash
   grakel.kernels.shortest_path
   grakel.kernels.shortest_path_attr
   grakel.kernels.subgraph_matching
   grakel.kernels.weisfeiler_lehman
   grakel.kernels.subtree_wl
   grakel.kernels.hadamard_code
   grakel.kernels.neighborhood_subgraph_pairwise_distance
   grakel.kernels.lovasz_theta
   grakel.kernels.svm_theta
   grakel.kernels.propagation
   grakel.kernels.odd_sth
   grakel.kernels.multiscale_laplacian
   grakel.kernels.multiscale_laplacian_fast
   grakel.kernels.hadamard_code
   grakel.kernels.vertex_histogram
   grakel.kernels.edge_histogram

**User guide:** See the :ref:`kernels` section for further details.

:mod:`grakel.datasets`: Datasets
=================================

Fetch
-----

.. currentmodule:: grakel.datasets

.. autosummary::
   :toctree: generated/
   :template: function.rst
   
   grakel.datasets.fetch_dataset
   grakel.datasets.get_dataset_info

**User guide:** See the :ref:`datasets` section for further details.

.. _gd:	https://ls11-www.cs.tu-dortmund.de/staff/morris/graphkerneldatasets