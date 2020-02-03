.. _installation:

=================
Installing GraKeL
=================
The GraKeL library requires the following packages to be installed:

* Python (>=2.7, >=3.5)
* NumPy (>=1.8.2)
* SciPy (>=0.13.3)
* Cython (>=0.27.3)
* cvxopt (>=1.2.0) [optional]
* future (>=0.16.0) (for python 2.7)

*GraKeL* is available via `PyPI`_ . You can install the latest release of *GraKeL* using the following command:

.. code-block:: bash

   $ pip install grakel

To also install the cvxopt package, which is a requirement of the Lovász-:math:`\vartheta` kernel, you can use the following command:

.. code-block:: bash

   $ pip install grakel[lovasz]

.. *GraKeL* is also available via `anaconda`_.

Building GraKeL
---------------

In order to build your own version of *GraKeL*, you need a C++ compiler since the package contains some C++ extensions. To build and install a local version of `GraKeL`, you need to execute :code:`pip install .` or :code:`python setup.py install` on the root folder. Furthermore, in case you want to build the extensions locally, execute :code:`python setup.py build_ext`.

In order for the C++ extensions to compile our extensions, a system-specific build environment should be configured. What you generally need is a C++ compiler and some python header files.

Unix Environment
^^^^^^^^^^^^^^^^^

In the case of Unix environments, you need to have installed:

- A C++ compiler like `g++`
- The package that contains the `Python.h` file such as `python-dev`

Windows Environment
^^^^^^^^^^^^^^^^^^^

In the case of a Windows environment, you need to install parts of the Windows Virtual Studio SDK (for C++) (for more details, please have a look here_).

.. note::

   If you have trouble building `GraKeL`, please raise an issue_ so that we can enrich our installation instructions, as well as addressing the problem.

Why so Many Packages?
---------------------

Graph kernels deal with the problem of graph comparison, a very challenging problem which has been studied for decades. Due to the complex nature of the problem, different types of approaches have been developed so far. Some approaches employ combinatorial algorithms, others formulate the graph comparison algorithm as a continuous optimization problem, while there are also other approaches that apply heuristics. The field of graph kernels is also characterized by such a large diversity of methods. For instance, the *graphlet kernel* solves the graph isomorphism problem to determine the identity of each graphet, while the *Lovász*-:math:`\vartheta` kernel solves a semidefinite programming problem to compute the Lovász number of each graph and the associated orthonormal representations. To solve such problems, *GraKeL* relies on well-established external libraries that provide optimized software that has been developed to address these problems. For example, *GraKeL* uses [bliss]_ to test graph isomorphism and the cvxopt_ library to optimize semidefinite programs.

.. _cvxopt: https://cvxopt.org/

.. [bliss] To test graph isomorphism, *GraKeL* extended `PyBliss`_, a Python wrapper for bliss. This allowed *GraKeL* to remain compatible with Python 2/3 and its installation on Windows. Among all the candidate packages, PyBliss was chosen thanks to the information shared by `Tamás Nepusz`_ (developer of the `iGraph`_ library), who pointed out that this package was the most efficient (both in terms of time and memory) for deciding isomorphism between small graphs in experiments conducted using the iGraph library. Other candidate packages include `pynauty`_ (a Python extension of `nauty`_) and `networkx`_ (contains an implementation of the `VF2`_ algorithm).

.. _PyBliss: http://www.tcs.hut.fi/Software/bliss/
.. _Tamás Nepusz: http://hal.elte.hu/~nepusz/
.. _iGraph: http://igraph.org/
.. _pynauty: https://web.cs.dal.ca/~peter/software/pynauty/html/
.. _nauty: http://users.cecs.anu.edu.au/~bdm/nauty/
.. _networkx: https://networkx.github.io/
.. _VF2: https://networkx.github.io/documentation/networkx-1.10/reference/algorithms.isomorphism.vf2.html
.. _PyPI: https://pypi.org/project/grakel-dev/
.. _anaconda: https://anaconda.org/ysig/grakel-dev
.. _issue: https://github.com/ysig/GraKeL/issues
.. _here: https://docs.microsoft.com/en-us/visualstudio/python/working-with-c-cpp-python-in-visual-studio?view=vs-2019#prerequisites
