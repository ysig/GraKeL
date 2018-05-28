.. _contributing:

============
Contributing
============

Ways you can contribute
-----------------------
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

* **Redesigning the** :code:`kernel` **class**: The :code:`kernel` class has been confined to the project constraints (namely sklearn compatibility)
  to the most abstract way possible. Other features and methods can be added to the kernel class to support frameworks that are currently
  not supported, such as *deep-graph-kernels* or *valid-assignment-kernels* if this is considered important by the community.

* **Unit-Testing**: As far as the kernel module is concerned we have not found a standard way of checking that the computations are valid
  either than doing a kernel validation by-hand or using code-of-others to see that we get similar results in computation, while controlling
  the classification accuracy concerning kernel-matrix-computation in benchmark datasets. Either than checking that a kernel is positive-semi-definite
  through the use of :code:`fit_transform`, we would like to be able to apply a more thorough unit-testing, that will expose the kernels calculation
  validity in relation to the theoretical foundation.

* **Utilizing Parallel/Concurrent computation**: It is important that the :code:`GraphKernel` decorator through the objects of the :code:`kernel`
  class supports a parallel computation scheme for executing procedures faster if the user wants so. This feature appears as currently implemented,
  but actually its implementation is wrong. It is a major problem to tackle that on the other hand, should not be so difficult.

.. _master: https://github.com/ysig/GraKeL
.. _develop: https://github.com/ysig/GraKeL/tree/develop


Who to blame for the GraKeL project
-----------------------------------
The GraKeL project is currently developed by the editor [Ioannis Siglidis] inside a one year project founded from the `Labex DigiCosme`_
organization of the *Paris-Saclay* University at the `DaSciM`_ laboratory of LiX.

The project is currently carried out under the administrative supervision of the DaSciM Team Leader, `Michalis Vazirgiannis`_ as well as the scientific
supervision of Giannis Nikolentzos a postdoctoral researcher.
Other members of the lab aiding with the construction of this library are: Xristos Giatsidis, Stratis Limnios and Konstantinos Skianis.

The project BSD license - a prerequisite of *sk-learn* compatibility - is accepted by the editor as it does not offend the virtues of the free-software community.

.. _DaSciM: http://www.lix.polytechnique.fr/dascim/
.. _Michalis Vazirgiannis: http://www.lix.polytechnique.fr/~mvazirg/
.. _Labex DigiCosme: https://digicosme.lri.fr/tiki-index.php
