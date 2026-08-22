.. grakel documentation master file, created by
   sphinx-quickstart on Mon Jan 18 14:44:12 2016.

========
Overview
========

*GraKeL* is a Python package which provides implementations of several graph kernels, a family of powerful methods which allow kernel-based learning approaches such as SVMs to work directly on graphs.

.. toctree::
   :maxdepth: 2

   documentation

==========
What's New
==========

- Version **0.1.12**

  + :code:`networkx_from_graph` converts a :code:`grakel.Graph` back into a
    NetworkX graph, completing the round trip with
    :code:`graph_from_networkx`. NetworkX is now an optional dependency,
    installable with :code:`pip install grakel[networkx]`, and support is
    pinned to NetworkX >= 3.0.
  + Fixed a label-conversion bug in :code:`Graph` that silently dropped edge
    labels when converting an adjacency-format graph with edge-only labels to
    the dictionary format.

- Version **0.1.11**

  .. note::

     From this release onwards **GraKeL is Python 3 only**. The last release
     that carried any Python 2 support was 0.1.10; if you are still on Python
     2, pin ``grakel<=0.1.10``. Supported versions are 3.9 to 3.12.

  + The Python 2 compatibility layer is gone: :code:`six` and :code:`future`
    are no longer dependencies, the :code:`__future__` imports and the
    :code:`six.moves` shims have been replaced by their standard library
    equivalents, and the Python 2 C-API wrapper for BLISS has been removed
    along with it. No behaviour changes beyond the fix below.
  + Fixed a latent error in :code:`priority_dict.smallest`, which used the
    Python 2 :code:`raise(IndexError, ...)` form. On Python 3 that raises
    :code:`TypeError: exceptions must derive from BaseException` instead of
    the intended :code:`IndexError`. The path is reached through the Dijkstra
    and shortest-path kernels on an empty queue.
  + :code:`setup.py` no longer imports :code:`distutils`, which was dropped
    from the standard library in 3.12.
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

  + Added a new kernel:
    `Weisfeiler-Lehman-Optimal-Assignment
    <https://ysig.github.io/GraKeL/0.1a8/kernels/weisfeiler_lehman_optimal_assignment.html>`_.
  + Removed :code:`MultiScaleLaplacian` as being really slow, and renamed
    :code:`MultiScaleLaplacianFast` to :code:`MultiScaleLaplacian`.
  + Fixed minor issues carried over from 0.1a7 (joblib deprecation, skbunch and
    the like).

- Version **0.1a7**

  + Detailed installation instructions for the c++ extensions on windows.
  + Renamed the :code:`base_kernel` argument of the frameworks to
    :code:`base_graph_kernel`, to disambiguate it from the vectorial kernels.
  + Faster floyd-warshall calculation in :code:`graph.py`.
  + Large update throughout the documentation.

- Version **0.1a6**

  + More scikit-learn compatibility: kernels can be initialised by name and
    alias on :code:`GraphKernel`, fitting and instantiation work from the
    default parameters, the random number generator is standardised on
    :code:`check_random_state` so :code:`random_seed` arguments are now
    :code:`random_state`, and the docstrings carry doctests.
  + More detailed output when a kernel is unsupported.
  + More detailed licensing information covering **cvxopt** and **BLISS**.
  + Bugfix in the (count sensitive) neighborhood hash kernel.
  + Sparse input support for :code:`VertexHistogram` and :code:`EdgeHistogram`.

- Version **0.1a5**

  + Various bugfixes in the kernel implementations.
  + Added :code:`utils` functions for external operations: converting existing
    *graph formats* (csv, pandas, networkx) to the grakel native one, *k-fold
    cross validation* with an SVM, and a *kernel matrix transformer* for
    manipulating precomputed kernel matrices in a :code:`Transformer` fashion.
  + **Conda** compatibility: visit `anaconda.org/ysig/grakel-dev
    <https://anaconda.org/ysig/grakel-dev>`_.

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

.. toctree::
   :maxdepth: 2

   api
   classes
   auto_examples/index
   tutorials


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
