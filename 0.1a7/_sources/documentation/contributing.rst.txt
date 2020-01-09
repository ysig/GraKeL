.. _contributing:

============
Contributing
============

All contributions are welcome! If you are not sure about how you can contribute, please contact the authors of the library.

Areas you can contribute
------------------------
Curious about how you can contribute to *GraKeL*? Here are a few ideas!

* **Implementing a kernel**: The number of graph kernels that have been proposed in the past years is very large. *GraKeL* contains implementations of several of these kernels, but still, there are many kernels that are not contained in the library. You can help making *GraKeL* more complete by implementing new graph kernels.

* **Optimizing kernel computation**: A lot of kernels are implemented in this library, which will increase in the future
  if users get interested contributing to this project. Python is a language that has balanced high-level and good-looking
  code of an imperative and object-oriented manner with computational efficiency. Although this is true, usually it cannot
  achieve the same computational achievements of other low-level languages, such as C and C++ or either VM optimizations
  such as those done in Java, being restricted to the script/interpreter execution scheme, which in some cases is outperformed
  by the pre-referenced programming languages. So except than computationally optimizing kernels through algorithm analysis
  which is generally done by the users or through the scientific community in time, it can be important to either re-implement
  existing kernel functions using wrapped C++ packages, that we either write ourselves or import from existing libraries
  like numpy, scipy, sklearn etc.

* **Improving the** :code:`Graph` **class**: As discussed in the :ref:`core_concepts` section, the :code:`Graph` class supports both adjacency matrix and edgelist representations of graphs. There are also methods that allow *GraKeL* to read graphs in various formats (e.g., NetworkX graphs). Furthermore, there are methods that implement graph algorithm (e.g., computation of shortest paths). These operations have to be efficient, both in terms of time and space complexity. Therefore, the :code:`Graph` class needs to be optimized. Altenatively, the project may benefit a lot from a *Cython* implementation.

* **Redesigning the** :code:`Kernel` **class**: The :code:`Kernel` class was designed in such a way that it satisfies some constraints (e.g., compatibility with *scikit-learn*) and is as simple as possible. This class can be extended to support families of methods that are not currently  *deep graph kernels*.

* **Unit-Testing**: As far as the kernel module is concerned, we have not managed to come up with any methodology for testing if the kernels are correctly implemented. We could use some "reference" code to check if our kernels produce identical results on some datasets, however, in most cases, this is not practical. Our tests check if the kernel matrices produced by the kernels are positive semidefinite, however, this can be true even if a kernel is not correctly implemented. We would like to design new tests that can verify the validity of implemeted kernels.

* **Utilizing Parallel/Concurrent computation**: It is important that the :code:`GraphKernel` generic-wrapper supports a parallel computation scheme,
  through the objects of the :code:`Kernel` class, for executing procedures faster if the user wants so. This feature appears as currently implemented,
  but actually its implementation is wrong. It is a major problem to tackle that on the other hand, should not be so difficult.

* **Examples and tutorials**: Have you created an example or tutorial that makes use of the *GraKeL* library? Please let us know. We would be more than happy to include it in our list of examples or tutorials.

.. _master: https://github.com/ysig/GraKeL
.. _develop: https://github.com/ysig/GraKeL/tree/develop


Who to Blame for the GraKeL Project
-----------------------------------
The *GraKeL* project started in 2018 as part of a one year project funded by `Labex DigiCosme`_. The main contributor to *GraKeL*'s development is Ioannis Siglidis. Ioannis is also responsible for its maintenance. Giannis Nikolentzos is also an active contributor. The project was carried out under the supervision of Professor `Michalis Vazirgiannis`_ at the LIX laboratory of École Polytechnique. The following people have also contributed to the project: Christos Giatsidis, Stratis Limnios and Konstantinos Skianis.

License
-------
GraKeL is distributed under the **BSD 3-clause** license. The library makes use of the C++ source code of BLISS_ (a tool for computing automorphism groups and canonical labelings of graphs) which is **LGPL** licensed. Futhermore, the cvxopt_ package (a software package for convex optimization) which is an optional dependency of GraKeL is **GPL** licensed.

.. _Labex DigiCosme: https://digicosme.lri.fr/tiki-index.php
.. _Michalis Vazirgiannis: http://www.lix.polytechnique.fr/~mvazirg/
.. _BLISS: http://www.tcs.hut.fi/Software/bliss
.. _cvxopt: https://cvxopt.org/
