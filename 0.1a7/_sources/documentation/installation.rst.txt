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

To also install the cvxopt package, which is a requirement of the Lov{\'a}sz-:math:`\vartheta` kernel, you can use the following command:

.. code-block:: bash

   $ pip install grakel[lovasz]

*GraKeL* is also available via `anaconda`_.

===============
Building GraKeL
===============

In order to build your own version of *GraKeL*, you need a C++, as the package relies on some C++ extensions.
Build and installing a local `grakel` can be done by executing :code:`pip install .` or :code:`python setup.py install` on the parent folder.
Additionally if you want to build your extensions locally, you can run :code:`python setup.py build_ext`.

In order for the C++ extensions to compile a system-specific building environment should be setup.


Linux Environment
-----------------


Windows Environment
-------------------


---------------------
Why so many packages?
---------------------
Graph kernels deal with the problem of graph comparison, a very challenging problem which has been studied for decades. Due to the complex nature of the problem, different types of approaches have been developed so far. Some approaches employ combinatorial algorithms, others formulate the graph comparison algorithm as a continuous optimization problem, while there are also other approaches that apply heuristics. The field of graph kernels is also characterized by such a large diversity of methods. For instance, the *graphlet kernel* solves the graph isomorphism problem to determine the identity of each graphet, while the *Lov{\'a}sz*-:math:`\vartheta` kernel solves a semidefinite programming problem to compute the Lov{\'a}sz number of each graph and the associated orthonormal representations. To solve such problems, *GraKeL* relies on well-established external libraries that provide optimized software that has been developed to address these problems. For example, *GraKeL* uses [bliss]_ to test graph isomorphism and the cvxopt_ library to optimize semidefinite programs.

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

