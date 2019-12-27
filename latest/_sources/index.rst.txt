.. grakel documentation master file, created by
   sphinx-quickstart on Mon Jan 18 14:44:12 2016.

========
Overview
========

GraKeL is a Python package extension, for the study and use of an upcoming
area in data-mining and machine learning, known as graph kernels.

Getting Started

  .. toctree::
    :maxdepth: 2

    user_manual

For seeing this version of the **GraKeL** project on the relevant **GitHub** repository you can have a look on the `README <https://github.com/ysig/GraKeL/blob/develop/README.md>`_.

=============
Documentation
=============

A collection of all classes and functions important
for the use and understanding of the **GraKeL** package.

GrakeL provides

  .. toctree::
    :maxdepth: 1

    api
    classes
    auto_examples/index


==========
What's New
==========

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
