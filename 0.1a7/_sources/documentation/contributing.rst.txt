.. _contributing:

============
Contributing
============

How you can contribute
----------------------
The project currently is on it's first circle of development and as so its form should not be considered static.
The most simple way a user can contribute is by proposing a new kernel that should be implemented or by implementing
it his own and send as the kernel by e-mail on :code:`ioannis.siglidis ~ inria.fr` or - if he knows how to use **git**
- to *fork* the `master`_ or the `develop`_ page and make a pull-request on the `develop`_ page. Integrating a kernel
takes a certain amount of time, because of unit-test, documentation and decoration is needed, as well as theoretical validation.

Areas you can contribute
------------------------
A list of areas where contribution seems important are the following:

* **Optimizing kernel computation**: A lot of kernels are implemented in this library, which will increase in the future
  if users get interested contributing to this project. Python is a language that has balanced high-level and good-looking
  code of an imperative and object-oriented manner with computational efficiency. Although this is true, usually it cannot
  achieve the same computational achievements of other low-level languages, such as C and C++ or either VM optimizations
  such as those done in Java, being restricted to the script/interpreter execution scheme, which in some cases is outperformed
  by the pre-referenced programming languages. So except than computationally optimizing kernels through algorithm analysis
  which is generally done by the users or through the scientific community in time, it can be important to either re-implement
  existing kernel functions using wrapped C++ packages, that we either write ourselves or import from existing libraries
  like numpy, scipy, sklearn etc.

* **Optimization of the Graph class**: As discussed at :ref:`longer_introduction`, :code:`Graph` solves as a uniform representation
  that supports both vertex-oriented and edge-oriented representations of graphs, while resolving conversions, importing of
  various types of inputs and finally a small group of important operations concerning graph objects, while being at the same
  time memory and computationally efficient. To achieve that goal the code of the :code:`Graph` should be optimized in python
  and if it is possible and seems more efficiently to be totally converted to a *Cython* extension (for the use of C++ as a back-end).

* **Redesigning the** :code:`Kernel` **class**: The :code:`kernel` class has been confined to the project constraints (namely sklearn compatibility)
  to the most abstract way possible. Other features and methods can be added to the kernel class to support frameworks that are currently
  not supported, such as *deep-graph-kernels* or *valid-assignment-kernels* if this is considered important by the community.

* **Unit-Testing**: As far as the kernel module is concerned we have not found a standard way of checking that the computations are valid
  either than doing a kernel validation by-hand or using code-of-others to see that we get similar results in computation, while controlling
  the classification accuracy concerning kernel-matrix-computation in benchmark datasets. Either than checking that a kernel is positive-semi-definite
  through the use of :code:`fit_transform`, we would like to be able to apply a more thorough unit-testing, that will expose the kernels calculation
  validity in relation to the theoretical foundation.

* **Utilizing Parallel/Concurrent computation**: It is important that the :code:`GraphKernel` generic-wrapper supports a parallel computation scheme,
  through the objects of the :code:`Kernel` class, for executing procedures faster if the user wants so. This feature appears as currently implemented,
  but actually its implementation is wrong. It is a major problem to tackle that on the other hand, should not be so difficult.

.. _master: https://github.com/ysig/GraKeL
.. _develop: https://github.com/ysig/GraKeL/tree/develop


Who to Blame for the GraKeL Project
-----------------------------------
The *GraKeL* project started in 2018 as part of a one year project funded by `Labex DigiCosme`_. The main contributor to *GraKeL*'s development is Ioannis Siglidis. Ioannis is also responsible for its maintenance. Giannis Nikolentzos is also an active contributor. The project was carried out under the supervision of Professor `Michalis Vazirgiannis`_ at the LIX laboratory of École Polytechnique. The following people have also contributed to the project: Christos Giatsidis, Stratis Limnios and Konstantinos Skianis.

License
-------
GraKeL is distributed under the __BSD 3-clause__ license. The library makes use of the C++ source code of BLISS_ (a tool for computing automorphism groups and canonical labelings of graphs) which is __LGPL__ licensed. Futhermore, the cvxopt_ package (a software package for convex optimization) which is an optional dependency of GraKeL is __GPL__ licensed.

.. _Labex DigiCosme: https://digicosme.lri.fr/tiki-index.php
.. _Michalis Vazirgiannis: http://www.lix.polytechnique.fr/~mvazirg/
.. _BLISS: http://www.tcs.hut.fi/Software/bliss
.. _cvxopt: https://cvxopt.org/
