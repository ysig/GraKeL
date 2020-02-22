.. _myok:

========================
Creating your Own Kernel
========================

As mentioned in the :ref:`core_concepts` subsection, each kernel imported from the :code:`GraphKernel` generic wrapper
and found inside :code:`grakel.kernels` sub-package, inherits the :code:`Kernel` class found there.
In order to write any kernel that will be integrated (see :ref:`contributing`), in our package we would like the
user to inherit that class. This may seem restricting but is not. We will try to demonstrate this in the
following sections and we would like to mention the fact that the whole kernel package **is not final**, but must
be always considered final when deploying a new package, because all the kernels should have a unified common interface.

-------------------------------------------
Overriding the :code:`Kernel` class methods
-------------------------------------------
In order to start we will present the current interface of the :code:`Kernel` class
(public methods) and guide the user how he can write a simple base :code:`Kernel`, such as the *vertex-histogram-kernel*.

The *vertex-histogram-kernel*, defined in :cite:`sugiyama2015halting` p.4 section 2.3 is a simple kernel, that
calculates label histograms for each graph, that is: *counts the number of occurrences for each label
value and as a kernel between two graphs calculates the sum of products between frequencies of common
occurrences*.

To design this kernel let's first learn some things about the kernel class.

Which methods should I implement
--------------------------------
Each kernel should have the following methods:

* :code:`__init__` method implemented. This method should be always overrided.

* :code:`fit`. Calculate features of the reference dataset.

* :code:`fit_transform`. Calculate features and the kernel matrix on the reference dataset.

* :code:`transform`. Calculate the kernel matrix between the fitted dataset and the transformed.
* :code:`diagonal`. Calculate the diagonal of the kernel matrix matrix produced between all elements of transform.

So let's start with our example and define an :code:`__init__` method for our *vertex-histogram-kernel*.

Let's import some library elements we will use

.. literalinclude:: code_for_examples/vertex_kernel.py
   :language: python
   :lines: 1-3

and let's define the *vertex-histogram-kernel*

.. literalinclude:: code_for_examples/vertex_kernel.py
   :language: python
   :lines: 6-59

As you see there a few things we should note. Firstly the :code:`verbose`, :code:`normalize`,
:code:`n_jobs` parameters should be always included in the kernel definition and passed as is,
to the father method by calling his :code:`__init__` method. All :code:`__init__` parameters should be stored
to class method parameters with attributes of the exactly same name, as this is important for the sci-kit learn
transformer **Pipeline** and any needed initialization registered in :code:`initialized_` attribute must be registered
inside :code:`initialize_` method as an update of the father object. If the user decides to overwrite :code:`fit` or :code:`fit_transform` he
must call :code:`self.initialize_()` as the first command of any process inside the kernel. This allows a valid
parameter *resetting* through :code:`set_params`, satisfying an important aspect for **Pipeline** consistency. 

Irrelevant from the interface the :code:`Kernel` class offers two methods, overriding which, a user can write
his own kernel. These methods are :code:`parse_input` and :code:`pairwise_operation`.

The :code:`parse_input`, :code:`pairwise_operation` methods
-----------------------------------------------------------
The :code:`parse_input` method actually produces all the features used from the kernel.
The default procedure is that both :code:`fit` and :code:`transform` use :code:`parse_input`
to create features for each graph producing a list where each element corresponds to the index
indicating that graph, from the corresponding iterable (while ignoring empty elements).
Those lists are then used feeding single list elements to the :code:`pairwise_operation`
method, calculating a final kernel value.

.. note::
    In order not to repeat again-again the word *kernel*, when defining a kernel class
    this word is omitted on the kernel name definition.

So to continue or example we first define :code:`parse_input` inside the :code:`VertexHistogram` class:

.. literalinclude:: code_for_examples/vertex_kernel.py
   :language: python
   :lines: 61-109

The procedure of reading the input is really standard, but must be specified for each kernel
according to its minimum acceptable input, that is the least elements with which it can be computed.
After reading the labels, we calculate frequencies with a :code:`Counter` on the label values, while
adding them to a list for all the non empty element of the iterable. Finally on those list
elements that are produced from :code:`parse_input` the pairwise operation that calculates
the kernel value is

.. literalinclude:: code_for_examples/vertex_kernel.py
   :language: python
   :lines: 111-125

where we count the :math:`\sum_{l \in L_{i} \cap L_{j}} f^{i}_{l}*f^{j}_{l}`, using the property of a :code:`Counter`
object, returning a **0** value on empty occurrences.

Now going back to the fascinating water-example

.. code-block:: python

    >>> H2O = Graph([[0, 1, 1], [1, 0, 0], [1, 0, 0]], {0: 'O', 1: 'H', 2: 'H'})
    >>> H3O = Graph([[0, 1, 1, 1], [1, 0, 0, 0], [1, 0, 0, 0], [1, 0, 0, 0]], {0: 'O', 1: 'H', 2: 'H', 3:'H'})
    >>> vh_kernel = vertex_histogram(normalize=True)
    >>> vh_kernel.fit_transform([H2O])[0, 0]
    1.0
    >>> vh_kernel.transform([H3O])[0, 0]
    0.9899494936611665

So this kernel works!

.. Note::
    The :code:`normalize` argument uses :code:`diagonal`, which makes it a good way of checking if this method is implemented also (correctly).


Why should I follow the list format on parse input?
---------------------------------------------------
You can still avoid changing the list format result on :code:`parse_input`, by instead of overriding pairwise operation, doing so with overriding :code:`_calculate_kernel_matrix`
method. This method must be able to receive one element whose default value is :code:`None`, where if it is, the method should calculate the kernel matrix between elements of :code:`self.X` (where the output of :code:`parse_input` is stored) and otherwise calculate between the input of :code:`Y` (which also comes from :code:`parse_input`) and :code:`self.X`.
By doing so, the developer needs to override the :code:`diagonal` method, because it normally uses :code:`pairwise_operation`. The diagonal method can also see **the last transformed
input as resulted from :code:`parse_input`** inside the attribute :code:`self._Y`.

But why not always use a pairwise operation?
--------------------------------------------
There are cases where the matrix calculation as a batch is much more computationally efficient.
This can be demonstrated, by utilizing libraries such as :code:`numpy`.

Let's go back to the vertex histogram kernel and do so. The class definition of :code:`__init__` stands as above.
We will redefine the :code:`parse_input` method, in order for utilizing the :code:`numpy` package.

Let's start with some new import statements

.. literalinclude:: code_for_examples/vertex_kernel_advanced.py
   :language: python
   :lines: 1-4

and define the new :code:`parse_input`

.. literalinclude:: code_for_examples/vertex_kernel_advanced.py
   :language: python
   :lines: 7-139

Let's stop here at two points.

* What is :code:`self._method_calling`? The parse input is a global method used by three basic functions (:code:`fit`, :code:`fit_transform` and :code:`transform`)
  to extract certain features, from the input while parsing it. Although in most cases the procedure is unique there are minor points where the execution
  of the algorithm should differentiate, according to if the information is intended for the fit method or not. :code:`self._method_calling` is a class attribute
  that should be initialized in *1* if any *method* is called from :code:`fit`, to *2* if it is called from :code:`fit_transform` and to *3* if it is called from
  :code:`transform`. The ordinary :code:`fit_transform`, calls :code:`fit` to parse the input in which case :code:`self._method_calling` will never be set to *2*, but it is added in favor of *completeness*.

* Why are you doing something else on :code:`fit` and on :code:`transform`? Generally :code:`fit` stands for setting the reference dataset, that is, it holds information that
  act as reference for any further transformation. In our case the :code:`labels`, appearing in all graphs are enumerated in order to be added on a feature matrix.
  The labels corresponding to the dataset of :code:`fit` should occupy the first position of the feature matrix, where any new labels of the matrix should be
  added on new positions coherently and forgotten (because there is no need to be remembered).

Now we would like to define the :code:`_calculate_kernel_matrix` method, calculating simple the kernel matrix as a result of a simple dot product

.. literalinclude:: code_for_examples/vertex_kernel_advanced.py
   :language: python
   :lines: 141-165

as well as the diagonal method (using Einstein summation convention which is a really fast method for speeding up
the diagonal calculation of the final matrix, as found `here`_)

.. literalinclude:: code_for_examples/vertex_kernel_advanced.py
   :language: python
   :lines: 167-193

Now let's solve a classification problem on the **"DD"** dataset, by following a standard procedure.
First import and download the dataset

.. code-block:: python

    >>> from grakel import GraphKernel, datasets
    >>> DD = datasets.fetch_dataset("DD", verbose=False)
    >>> DD_data = DD.data
    >>> DD_classes = DD.target

afterwards split train/test in a pseudo-random way holding 

.. code-block:: python

    >>> from sklearn.model_selection import train_test_split
    >>> X, Y, y_train, y_test = train_test_split(DD_data, DD_classes, test_size=0.1, random_state=42)

.. note::
    For understanding why 42, was chosen as a :code:`random_state`, check `this`_.

now let's calculate the :code:`K_train` and :code:`K_test` kernels as

.. code-block:: python

    >>> vh_kernel = vertex_histogram(normalize=True)
    >>> K_train = vh_kernel.fit_transform(X)
    >>> K_test = vh_kernel.transform(Y)

and finally classify, by using/finding the best C (that is a parameter emphasizing to the SVM,
`how much you want to avoid misclassifying each training example`_)

.. code-block:: python

    >>> from sklearn.metrics import accuracy_score
    >>> from sklearn.svm import SVC
    >>> from numpy import arange

    >>> C_grid = 10. ** arange(4,10,1) / len(DD.data)

    >>> acc_score = 0
    >>> for i in range(C_grid.shape[0]):
    >>>     clf = SVC(kernel='precomputed', C = C_grid[i])
    >>>     clf.fit(K_train, y_train)
    >>>     y_pred = clf.predict(K_test)
    >>>     acc_score = max(acc_score, round(accuracy_score(y_test, y_pred)*100, 2))

which produces an accuracy score close to the maximum 78.2% documented on :cite:`kriege2016valid`

.. code-block:: python

    >>> print(acc_score, "%")
    76,27 %

.. _here: https://stackoverflow.com/questions/14758283/is-there-a-numpy-scipy-dot-product-calculating-only-the-diagonal-entries-of-the
.. _this: https://www.google.fr/search?q=the%20answer%20to%20life%20the%20universe%20and%20everything
.. _how much you want to avoid misclassifying each training example: https://stats.stackexchange.com/questions/31066/what-is-the-influence-of-c-in-svms-with-linear-kernel


What if I don't want to follow, this structure at all?
------------------------------------------------------
Although the basic methods (:code:`fit`, :code:`fit_transform`, :code:`transform`, :code:`diagonal`) are needed for the kernel
pipeline, the kernel structure is not imposed, but proposed. The user can always write her/his methods, as long as they are coherent with some design specifications needed for a *valid* kernel, mainly:

* Unified input support between all structures, while expecting input in a relevant way, that is accepting a :code:`Graph` type-object
  or an iterable of at least one element, which position **0** should signify the graph Object, position **1** should signify the
  node labels and position **2** the edge labels.
* Support of :code:`normalization` in a valid way is imperative.
* The method :code:`diagonal` is used in order for nested-graph-kernels (such as *Weisfeiler-Lehman*) to be able to calculate a correct
  *normalized* output, so it is also imperative, to be valid. Its output should be supported both as a scalar and as a vector.


---------------------
Getting more advanced
---------------------
.. note::
    This documentation entry is a preamble to the :ref:`contributing` section and will evolve much more as users contribute to this project.

What if you do not want just to write your kernel, but fork the project in order to integrate an optimized version of your new kernel.

Wrapping C++ functions for the kernels to use
---------------------------------------------
We have implemented a simple `Cython`_ extension for the use of **C++** inside python kernel implementation.
This package can be found inside :code:`grakel/kernels/_c_functions`. There the user can find the following files.

::

    grakel/kernels/_c_functions
                        ├── functions.cpp
                        ├── functions.pxd
                        ├── functions.pyx
                        ├── include
                        │      └── functions.hpp
                        └── src
                             └── *.cpp

In order to add a new function, what he/she should do is to add its source file inside :code:`src` while including/implementing
its definition found inside the :code:`include/functions.hpp`. Afterwards by externing the function definition inside the library
file :code:`functions.pxd`, by importing it from :code:`include/functions.hpp` you should add it inside the functions.pyx, file
where the essential IO wraping of the C++ function is done inside Python. Finally you must add the address of new *C++* you wrote on
the :code:`_c_functions` :code:`Extension` found on :code:`setup.py`, by adding it to the list given to the :code:`sources` argument.

.. _Cython: http://cython.org/

Parallelization
---------------
A basic infrastracture for utilizing parallelization possibilities of the kernel computation, was created inside the kernel class.
The user define n_jobs for each kernel and the indexes of the kernel matrix are *linearized* and splitted to the effective number of jobs
on which the pairwise operation is applied in parallel. We use the :code:`joblib` library, as found inside the :code:`sklearn.externals`.
On certain frameworks (Weisfeiler Lehman, Hadamard Code) the strategy for applying parallel kernel calculation is different, where the
task is the single kernel calculation of a kernel matrix from the base kernel. Parallelization is a feature we wanted to include in our
library because it can give the opportunity to the user to increase the efficiency of kernel matrix computation on large dataset, **although
the induced overhead in the most cases seems not to worth it**. The developer can either follow our road and call the initialization
of the father method creating a :code:`joblib.Parallel` object with a :code:`threading` *backend api* (stored in :code:`_parallel`) and use it inside
her/his code anyway he wants or follow a different strategy that seems to have a computational advantage. In case someone wants to contribute
in redisigning and extending the parallelization proccess he/she can see how in :ref:`contributing`.


Bibliography
------------
.. bibliography:: ../biblio.bib
   :filter: docname in docnames
