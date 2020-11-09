.. grakel documentation master file, created by
   sphinx-quickstart on Mon Jan 18 14:44:12 2016.

========
Overview
========

*GraKeL* is a Python package which provides implementations of several graph kernels, a family of powerful methods which allow kernel-based learning approaches such as SVMs to work directly on graphs.

Getting Started

  .. toctree::
    :maxdepth: 2

    documentation

==========
Benchmarks
==========

To demonstrate the efficiency of the algorithms implemented in *GraKeL*, we present a comparison of the running times of the implementations of some graph kernels from *GraKeL* and from other packages. We also compare the running times of the different kernels to each other.

  .. toctree::
    :maxdepth: 2

    benchmarks

=================
Package Reference
=================

A collection of all classes and functions important for the use and understanding of the *GraKeL* package.

GraKeL provides

  .. toctree::
    :maxdepth: 1

    api
    classes
    auto_examples/index
    tutorials


==========
What's New
==========

- Version **0.1a8**

  + Added a new kernel: [Weisfeiler-Lehman-Optimal-Assignment](https://ysig.github.io/GraKeL/0.1a8/kernels/weisfeiler_lehman_optimal_assignment.html).
  + Removed MultiScaleLaplacian (as being really slow and useless) and renamed MultiScaleLaplacianFast to MultiScaleLaplacian.
  + Fixed minor issues (joblib deprecation, skbunch etc) from `0.1a7`.

- Version **0.1a7**

  + Detailed installation instructions for c++ extensions in windows.
  + Changed `base_kernel` alias in frameworks with `base_graph_kernel` to disambiguate with vectorial kernels.
  + Speed-up for floyd_warshall calculation in graph.py.
  + Large update throughout all the documentation.

- Version **0.1a6**

  + More scikit-learn compatibility:

    1. Initialise kernels by name and alias on GraphKernel (as GraphKernel(kernel="shortest_path").
    2. Fit and instantion by default parameters.
    3. Random number generator standardized `check_random_state`. `random_seed` are now `random_state` arguments.
    4. Doctests.

  + Miscelanous: 

    1. Detailed unsupported kernel output.
    2. More detailed licensing information considering **cvxopt** and **BLISS**
    3. Small bugfix inside the (Count Sensitive) Neighborhood Hash Kernel.
    4. Added sparse-compatibility for VertexHistogram and for EdgeHistogram.

- Version **0.1a5**

  + Various bugfixes in kernel implementations.
  + Added a bunch of :code:`utils` functions for external operations: transforming existing *graph formats* (csv, pandas, networkx) to the grakel native, *k-fold cross validation* with an SVM and *kernel matrix transformer* for manipulating precomputed kernel matrices in an :code:`Transformer` fashion.
  + **Conda** compatibility: visit `<https://anaconda.org/ysig/grakel-dev>`_.

==================
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
