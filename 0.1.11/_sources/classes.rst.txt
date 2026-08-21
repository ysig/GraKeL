.. _api_ref:

=============
API Reference
=============

This is the class and function reference of *GraKeL*. In order for the user to understand how to use the package, we suggest he reads :ref:`documentation` section.

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

   Kernel
   RandomWalk
   RandomWalkLabeled
   PyramidMatch
   NeighborhoodHash
   ShortestPath
   ShortestPathAttr
   GraphletSampling
   SubgraphMatching
   WeisfeilerLehman
   HadamardCode
   NeighborhoodSubgraphPairwiseDistance
   LovaszTheta
   SvmTheta
   Propagation
   PropagationAttr
   OddSth
   MultiscaleLaplacian
   MultiscaleLaplacianFast
   HadamardCode
   VertexHistogram
   EdgeHistogram
   GraphHopper
   CoreFramework
   WeisfeilerLehmanOptimalAssignment

**User guide:** See the :ref:`kernels` section for further details.

:mod:`grakel.datasets`: Datasets
=================================

Fetch
-----

.. currentmodule:: grakel.datasets

.. autosummary::
   :toctree: generated/
   :template: function_bib.rst

   fetch_dataset

.. autosummary::
   :toctree: generated/
   :template: function.rst

   get_dataset_info


**User guide:** See the :ref:`datasets` section for further details.


:mod:`grakel`: Utils
=================================

.. currentmodule:: grakel

Use a kernel matrix as a transformer
------------------------------------

.. autosummary::
   :toctree: generated/
   :template: class.rst

   KMTransformer

Cross Validation
----------------

.. autosummary::
   :toctree: generated/
   :template: function.rst

   cross_validate_Kfold_SVM

Load from other file formats
----------------------------

.. autosummary::
   :toctree: generated/
   :template: function.rst

   graph_from_networkx
   graph_from_pandas
   graph_from_csv

**User guide:** Usefull functions for applying to existing datasets, of other formats.

.. _gd:	https://ls11-www.cs.tu-dortmund.de/staff/morris/graphkerneldatasets
