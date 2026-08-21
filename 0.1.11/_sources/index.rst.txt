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

- Version **0.1.11**

  + Python support is now **3.9 to 3.12**; wheels are built for all of them on
    Linux, macOS and Windows. Python 2 is gone from CI, the docs and the conda
    recipe.
  + Dependency minimums raised to the first versions supporting 3.9: numpy
    1.19, scikit-learn 0.24, Cython 0.29.21 and scipy 1.12. **scipy** is now a
    declared dependency rather than one picked up transitively.
  + Fixes for modern numpy and scipy, which the kernels had stopped working
    against: big-endian sparse dtypes in :code:`VertexHistogram`, the removed
    :code:`eigvals` keyword in :code:`SvmTheta` and :code:`LovaszTheta`, and
    scalar extraction under numpy 2 in :code:`SubgraphMatching`.
  + :code:`WeisfeilerLehman` no longer fails on a node with no outgoing edges
    (`#119 <https://github.com/ysig/GraKeL/issues/119>`_). Note that it
    relabels from outgoing edges, so on directed graphs incoming edges do not
    take part.
  + :code:`NeighborhoodSubgraphPairwiseDistance` can be used as a base kernel
    of a framework again (`#120
    <https://github.com/ysig/GraKeL/issues/120>`_): its :code:`diagonal()`
    returned a scalar where one value per graph was expected.
  + Printing a :code:`Graph` works again (`#102
    <https://github.com/ysig/GraKeL/issues/102>`_).
  + The geometric :code:`RandomWalk` kernel keeps its series convergent: the
    decay factor is checked against the spectral radius of the graphs it is
    given and lowered, with a warning, when it would diverge. Previously this
    produced negative self similarities and NaNs once normalized.
  + :code:`GraphHopper` says when it has been handed discrete labels instead of
    continuous node attributes (`#106
    <https://github.com/ysig/GraKeL/issues/106>`_), rather than failing with a
    shape mismatch.
  + Tooling: the project builds with **uv**, CI runs on github actions with the
    docs built and published from there, and the release is a single manual
    workflow. CircleCI is gone.

- Version **0.1.10**

  + Fixes for a batch of reported bugs: :code:`EdgeHistogram` errors (`#97
    <https://github.com/ysig/GraKeL/issues/97>`_), one against many comparison
    with :code:`WeisfeilerLehman.transform` (`#95
    <https://github.com/ysig/GraKeL/issues/95>`_),
    :code:`RandomWalkLabeled` returning an all ones matrix (`#96
    <https://github.com/ysig/GraKeL/issues/96>`_),
    :code:`NeighborhoodSubgraphPairwiseDistance` diagonals below one (`#94
    <https://github.com/ysig/GraKeL/issues/94>`_), missing arguments on
    :code:`fit_transform` (`#75
    <https://github.com/ysig/GraKeL/issues/75>`_) and the initialisation of
    :code:`self.mu_` in the random walk kernel (`#71
    <https://github.com/ysig/GraKeL/issues/71>`_).
  + :code:`Graph.copy()` and a :code:`__repr__` for graphs (`#90
    <https://github.com/ysig/GraKeL/issues/90>`_).
  + Newer python and **cvxopt** versions; cp312 wheels were not yet buildable
    at the time of this release.

- Version **0.1.9**

  + Wheels are built and published from CI for every supported python, with
    win32 dropped following scipy, and musllinux skipped.
  + :code:`LovaszTheta` warns on windows, where the underlying solver is
    unreliable.

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

================
Acknowledgements
================

We would like to thank `@SneachChea <https://github.com/SneachChea>`_ for a round of tooling, CI, and docs modernization.

We would like to thank `@eddiebergman <https://github.com/eddiebergman>`_ for modernizing our CI and extending our python support.

==================
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
