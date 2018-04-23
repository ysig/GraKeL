.. _installation:

=================
Installing grakel
=================
The grakel library requires:

* Python [>=2.7 or >= 3.5]
* NumPy [>= 1.8.2]
* SciPy [>= 0.13.3]
* Cython [>= 0.27.3]
* cvxopt [>= 1.2.0] [optional: lovasz]

For installing NumPy, SciPy, cvxopt and Cython the procedure
is the well known: :code:`pip install extension>=extension_version`.

---------------------
Why so many packages?
---------------------
The field of computationally efficient `Graph Kernels`_ can be considered
more as a collection of techniques calculating `PSD`_ similarity matrices, between
graph objects, than a field (with the compact sense of the term). The big diversity
of methods and ideas implied to extract a similarity score between to graphs, leads
to the use of very precise and studied algorithms as graph isomorphism 
(for determining *isomorphic-graphlets* on the graphlet-sampling kernel [bliss]_)
or convex optimization (through the use of *semidefinite-programming* 
for calculating *lovasz-theta* embeddings, in a *lovasz-theta* kernel)
which on the other hand will **not** appear in any other kernel.
As a result import of external libraries that have studied and optimized the solutions
of such problems in detail, gives a complexity - implementation "standard" to refer to, 
while follows the idea of standing on the shoulders of the scientific community, which
can be viable through the use of `free software`_.

.. _Graph Kernels: https://en.wikipedia.org/wiki/Graph_kernel
.. _PSD: https://en.wikipedia.org/wiki/Positive-definite_matrix
.. _free software: https://en.wikipedia.org/wiki/Free_software

.. [bliss] GraKeL extended `pybliss`_ in order to be compatible with Python2/3 and Windows installation
    while achieving specifically the task of deciding graph isomorphism. This package was choosen thanks to
    the information of `Tamás Nepusz`_ (developer of `iGraph`_), who pointed out the internal use of this
    package by `iGraph`_, for calculating decisional isomorphism in small graphs which was the package with
    the highest efficiency (both in time and memmory) inside experiments for this certain tasks. Other candidates
    where the `pynauty`_ package (using `nauty`_ practical isomorphism library) and `networkx`_ (implementing
    a `VF2`_ algorithm for graphs). Credits to the creators of this complex of useful tools!

.. _pybliss: http://www.tcs.hut.fi/Software/bliss/
.. _Tamás Nepusz: http://hal.elte.hu/~nepusz/
.. _iGraph: http://igraph.org/
.. _pynauty: https://web.cs.dal.ca/~peter/software/pynauty/html/
.. _nauty: http://users.cecs.anu.edu.au/~bdm/nauty/
.. _networkx: https://networkx.github.io/
.. _VF2: https://networkx.github.io/documentation/networkx-1.10/reference/algorithms.isomorphism.vf2.html
