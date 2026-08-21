.. _contributing:

============
Contributing
============

All contributions are welcome! If you are not sure about how you can contribute, please contact the authors of the library.

Areas you can contribute
------------------------
Curious about how you can contribute to *GraKeL*? Here are a few ideas!

* **Implementing a kernel**: The number of graph kernels that have been proposed in the past years is very large. *GraKeL* contains implementations of several of these kernels, but still, there are many kernels that are not contained in the library. You can help making *GraKeL* more complete by implementing new graph kernels.

* **Optimizing kernel computation**: We have done our best to write maintainable and efficient Python code. However, this does not mean that the code cannot be further optimized. For even higher efficiency, some graph kernels can be re-implemented using wrapped C++ packages. Furthermore, most kernels solve combinatorial problems for which more efficient algorithms (than the employed ones) may exist. 

* **Improving the** :class:`grakel.Graph` **class**: As discussed in the :ref:`core_concepts` section, the :class:`grakel.Graph` class supports both adjacency matrix and edgelist representations of graphs. There are also methods that allow *GraKeL* to read graphs in various formats (e.g., NetworkX graphs). Furthermore, there are methods that implement graph algorithms (e.g., shortest path distances). These operations have to be efficient, both in terms of time and space complexity. Therefore, the :class:`grakel.Graph` class needs to be optimized. The project may benefit a lot from a *Cython* implementation of the class.

* **Redesigning the** :class:`grakel.Kernel` **class**: The :class:`grakel.Kernel` class was designed to satisfy some constraints (e.g., compatibility with *scikit-learn*) and to be as simple as possible. This class can be extended to support families of kernels or frameworks that are not currently developed such as *deep graph kernels*.

* **Unit-Testing**: As far as the kernel module is concerned, we have not managed to come up with any methodology for testing if the kernels are correctly implemented. We could use some "reference" code to check if our kernels produce identical results on some datasets, however, in most cases, this is not practical. Our tests check if the kernel matrices produced by the kernels are positive semidefinite, however, this can be true even if a kernel is not correctly implemented. We would like to design new tests that can verify the validity of implemeted kernels.

* **Parallel/Concurrent execution**: The :class:`grakel.GraphKernel` class supports a parallel computation scheme (i.e., using the :code:`n_jobs` attribute), but this has not been implemented for all the kernels or the current implementation is not optimal. Implementations that allow parallel computation of the kernels are of high importance since they can lead to significant speed-ups of kernel computations.

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
