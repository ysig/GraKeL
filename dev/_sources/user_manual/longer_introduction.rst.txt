.. _longer_introduction:

=====================
A longer Introduction
=====================

What is the `GraphKernel` class
-------------------------------
`GraphKernel` is a class decorator, which means that it takes a collection of classes, in our case Graph-Kernels and creates a uniform interface for all them, while providing the user with a way of adding various features.

These features can be listed as follows:

* :code:`kernel` : The kernel can be either a :code:`base_kernel` or a list of :code:`general_kernels` that end in a :code:`base_kernel`.

    - :code:`base_kernel` : 
        To initialize a :code:`base_kernel` kernel the procedure is simple. Any :code:`base_kernel` is a dictionary containing its name under the :code:`'name'` field and its parameterization on separate fields that signify kernel parameters and their values. The :code:`shortest_path` kernel we show on the introduction is such a kernel.

        .. note::
            The decorator sometimes wraps two kernels in one (as in the :code:`MultiscaleLaplacian` and :code:`MultiscaleLaplacianFast`) and in order to learn
            the meaning of each parameter the user is suggested to read the documentation found on :ref:`kernels`.

    - :code:`frameworks` : 
        This type of kernels takes as a :code:`base_kernel` another kernel object. This kernel is also a dictionary containing its name under the :code:`'name'` field and its 
        parameterization on separate fields that signify kernel parameters and their values. The kernel produced from all the rest kernels in the list is considered as a :code:`base_kernel` and its parametrization, can be applied on the list consecutive elements. The :code:`weisfeiler_lehman` we show on the introduction is such a kernel.

    If no parameters are given at parametrization, default values are assigned.

* :code:`normalize` : Any kernel has the ability to provide a normalized output. This is important because otherwise the user should calculate the whole kernel matrix at the :code:`fit_transform` stage in order to make a valid normalize representation of the kernel matrix. Code example:
    Say we have a :code:`G_train`, :code:`G_test`. This

    .. code-block:: python

        >>> gk = GraphKernel(kernel=[{"name": "weisfeiler_lehman", "niter": 5}, {"name": "subtree_wl"}], normalize=True)
        >>> # Calculate the kernel matrix.
        >>> K_train = gk.fit_transform(G_train)
        >>> K_test = gk.transform(G_test)

    should be equivalent as process (set aside a non-deterministic kernel or split), with this

    .. code-block:: python

        >>> gk = GraphKernel(kernel=[{"name": "weisfeiler_lehman", "niter": 5}, {"name": "subtree_wl"}], normalize=False)
        >>> # Calculate the kernel matrix.
        >>> K = gk.fit_transform(G)
        >>> K_diag = K.diagonal()
        >>> K_train_diag, K_test_diag = K_diag[train_diag], K_diag[test_diag]
        >>> K_train = K[train_indices, :][:, train_indices] / np.sqrt(np.outer(K_train_diag, K_train_diag))
        >>> K_test = K[test_indices, :][:, train_indices] / np.sqrt(np.outer(K_test_diag, K_train_diag))

    but in the second case we make some more computations, in reward that fit samples are drown from both
    train and test datasets, producing a different result in some kernels by unifying intuitively
    train and test data when creating features, which may be desired (e.g. on the :code:`MultiscaleLaplacianFast`).

* :code:`Nystroem` : Nystroem is very well known method, for approximating kernel matrices on huge datasets.
    The kernel matrix is calculated only in random drown subsets of graphs, whose size can be defined by the user 
    by setting an int value inside the :code:`Nystroem` parameter. After that, all the kernel values are calculated
    in the sample space. An example indicating the power of nystroem can be indicated below:

    | Example: We will perform a simple classification task.

    Download the dataset and split to train and test

    .. code-block:: python

        >>> from grakel import datasets
        >>> MUTAG = datasets.fetch_dataset("MUTAG", verbose=False)
        >>> MUTAG_data, y = MUTAG.data, MUTAG.target
        >>> split_point = int(len(MUTAG_data) * 0.9)
        >>> X, Y = MUTAG_data[:split_point], MUTAG_data[split_point:]

    Initialise a :code:`GraphKernel`, using :code:`Nystroem` of 20 samples

    .. code-block:: python

        >>> from grakel import GraphKernel
        >>> wl_kernel = GraphKernel(kernel = [{"name": "weisfeiler_lehman", "niter": 5}, {"name": "subtree_wl"}], Nystroem=20)
        >>> K_train = wl_kernel.fit_transform(X)
        >>> K_test = wl_kernel.transform(Y)
        >>> print(K_train.shape)
        (169, 10)
        >>> print(K_test.shape)
        (19, 10)


    Classify using a standard SVC

    .. code-block:: python

        >>> y_train, y_test = y[:split_point], y[split_point:]
        >>> from sklearn import svm
        >>> clf = svm.SVC()
        >>> clf.fit(K_train, y_train)
        SVC(C=1.0, cache_size=200, class_weight=None, coef0=0.0,
          decision_function_shape='ovr', degree=3, gamma='auto', kernel='rbf',
          max_iter=-1, probability=False, random_state=None, shrinking=True,
          tol=0.001, verbose=False)
        >>> y_pred = clf.predict(K_test)

    finnaly calculate accuracy score

    .. code-block:: python
        
        >>> from sklearn.metrics import accuracy_score
        >>> print(str(round(accuracy_score(y_test, y_pred)*100, 2)), "%")
        78.95 %

    .. note::
        | As we see the accuracy of the classification is the same, allthough instead of doing
        | ~ 169 * (169-1) /2 + 19 * 169 = 17,407 computations we did
        | ~ 20 * (20-1)/ 2 + 20 * 169 + 20* 19 = 3,950 computations.

* :code:`n_jobs` : Some kernels have operations that can be executed concurrently, making computation faster 
    when user uses a significant amount of data, to overcome the parallelization overhead. :code:`n_jobs` follows
    the same formulation as in scikit-learn where giving as input 0 or -1 :code:`n_jobs` signifies initializing all the
    possible workers and if given a positive number, initializes that amount of workers if possible. There are kernels
    where this feature is not implemented either from already using low level parallelization from other libraries (such as numpy)
    or there was not a way that applying parallelization seemed to *worth it*. In such case the kernel pops a specified warning.

    .. note::
        The efficiency of parallelization is generally revealed when doing kernel computation on large datasets where the
        final operation that calculates the kernel value (generally between features extracted from graphs) is the one of the
        computation bottlenecks of the whole operation.

    To extend these feature to more kernels or to propose new computational strategies see how you canc **contribute** in :ref:`contributing`.

* :code:`random_seed` : We would in generally want to satisfy the need of the user to provide
    a :code:`random_seed` either to kernels that are probabilistic, or to randomize accordingly
    procedures of the :code:`GraphKernel` that need randomization such as :code:`Nystroem`, where the
    decorator draws probabilistically a number of components from the number of fitted samples.

    Let's give an example of a probabilistic kernel using our old water example. We will use a very well known kernel called *Graphlet-Sampling*, where we will
    sample probabilistically 5 subgraph samples from each graph either :math:`\mathbf{H}_{2}\mathbf{O}` or :math:`\mathbf{H}_{3}\mathbf{O}^{+}`.

    After initializing the input

    .. code-block:: python

        >>> from grakel import GraphKernel
        >>> H2O = [[[[0, 1, 1], [1, 0, 0], [1, 0, 0]], {0: 'O', 1: 'H', 2: 'H'}]]
        >>> H3O = [[[[0, 1, 1, 1], [1, 0, 0, 0], [1, 0, 0, 0], [1, 0, 0, 0]], {0: 'O', 1: 'H', 2: 'H', 3:'H'}]]

    let's calculate a default kernel value

    .. code-block:: python

        >>> gs_kernel = GraphKernel(kernel=dict(name="graphlet_sampling", n_samples=5))
        >>> gs_kernel.fit(H2O)
        >>> gs_kernel.transform(H3O)
        20.0

    Note that if a random seed is not given as an argument either to the :code:`GraphKernel` or to the kernel parameters
    a default will be used. Now let's try to give one as a parameter of the kernel (say 3)

    .. code-block:: python

        gs_kernel = GraphKernel(kernel=dict(name="graphlet_sampling", n_samples=5, random_seed=3))
        gs_kernel.fit(H2O)
        gs_kernel.transform(H3O)
        10.0

    As we see a new value has been calculated because the default seed is now not used. If know a :code:`random_seed`
    is initialized inside the decorator and no parameter is given signifying a :code:`random_seed` to the :code:`kernel`
    argument then if the kernel has such parameter the default will be used. This is demonstrated in what follows

    .. code-block:: python

        gs_kernel = GraphKernel(kernel=dict(name="graphlet_sampling", n_samples=5), random_seed=3)
        gs_kernel.fit(H2O)
        gs_kernel.transform(H3O)
        10.0

    where we get the same result. Now if both a :code:`GraphKernel` has a :code:`random_seed` and the :code:`kernel` is provided
    with one inside parametrization, the second will be used inside the :code:`kernel` and the first outside, in the rest code area
    covered by the decorator, as expected. To demonstrate we will show is the following:

    .. code-block:: python

        >>> gs_kernel = GraphKernel(kernel=dict(name="graphlet_sampling", n_samples=5, random_seed=3), random_seed=10)
        >>> gs_kernel.fit(H2O)
        >>> gs_kernel.transform(H3O)
        10.0

    where

    .. code-block:: python

        >>> gs_kernel = GraphKernel(kernel=dict(name="graphlet_sampling", n_samples=5, random_seed=10))
        >>> gs_kernel.fit(H2O)
        >>> gs_kernel.transform(H3O)
        15.0


* :code:`verbose` : 
    .. note::
        Verbose is an argument that is currently unsupported (has no impact), but is set for future implementation of some output messages.

To understand what the :code:`GraphKernel` object is doing, we must see inherently what objects it decorates.

The `Kernel` class
------------------
This :code:`Object` is any object inherited from the :ref:`kernel` (which can be imported from :code:`grakel`).

Normally a kernel function, between graphs should be considered as a function with to arguments,
such as :math:`k \; : \; \mathcal{G} \times \mathcal{G} \rightarrow \mathbb{R}`.
This raises two issues, namely one of efficiency and one of compatibility:

1. The first one has to do with the fact, that there are major computational advantages if instead of calculating the kernel pairwise, we calculate the whole kernel matrix.

2. The second has to do with the fact, that we wanted our project to be integrable inside the `sk learn template`_. From this template the most relevant structure was the sci-kit transformer, which consists of three inherent methods: :code:`fit`, :code:`fit_transform`, :code:`transform`.

So the way we conceptually attached the kernel definition to that design pattern was:

- The :code:`fit` part should fix a graph dataset as the base of comparison calculating necessary features.

- The :code:`fit_transform` should fit and calculate the kernel matrix on the fitted dataset.

- The :code:`transform` should calculate the matrix produced between a new dataset (namely the *test*) and the fitted dataset.

The deconstruction of the kernel matrix calculation from a function :math:`\mathcal{K}: \mathcal{G}^{\text{train}} \times \mathcal{G}^{\text{test}} \rightarrow \mathbb{R}^{n_{\text{test}}} \times \mathbb{R}^{n_{\text{train}}}`
to a `currying`_ scheme :math:`\mathcal{K}: \mathcal{G}^{\text{train}} \rightarrow \mathcal{G}^{\text{test}} \rightarrow \mathbb{R}^{n_{\text{test}}} \times \mathbb{R}^{n_{\text{train}}}` is not always equivalent in the
result, if some of the data of :math:`\mathcal{G}^{\text{train}}`, must be combined with data of :math:`\mathcal{G}^{\text{test}}` to produce the fit reference-features. In such cases
as mentioned above, namely in the case of :code:`multiscale_laplacian`, if the user wants :math:`\mathcal{G}^{\text{train}} \rightarrow \mathcal{G}^{\text{test}}` to be concerned
before fit we advise him to use the :code:`fit_transform`, function in the whole of the train and test data and separate the kernel matrices on the result.

.. note::
    The very idea that lies before fitting concerns holding a reference dataset. This means a collections of features should be stored into memory and **not** get corrupted throughout various applications of :code:`transform`. This however - the need of copying and protecting the reference data - produces a computational overhead in kernels such as the :code:`odd_sth` where the user will may prefer the computational advantages of applying a sole :code:`fit_transform`.

Using a :code:`kernel` type object through the decorator, should be equivalent with doing so without the decorator, if the correct parametrization is given.
The decorator **does not** restrict any *user-oriented* interface of the kernels, except if the user wants to write a kernel of his own.
If you want to know more about the kernel structure in order to write your own see :ref:`myok`.

To demonstrate a small example of the above we will construct our own a WL-subtree kernel instead of using the decorator.
To do so first import the :code:`WeisfeilerLehman` and :code:`VertexHistogram` (where :code:`vertex_histogram` is equivalent
with the :code:`subtree_kernel`) kernels as

.. code-block:: python

    >>> from grakel import WeisfeilerLehman
    >>> from grakel import VertexHistogram

If we see the documentation of :ref:`weisfeiler_lehman`, we can see that it accepts two arguments upon initialization: a :code:`niter` and a :code:`base_kernel`. The :code:`base_kernel` is a tuple consisting of a :code:`kernel` type object and a dictionary of arguments. To initialize a Weisfeiler-Lehman with 5 iterations and a subtree base-kernel.

.. code-block:: python

    >>> wl_kernel = WeisfeilerLehman(niter=5, base_kernel=(VertexHistogram, {}))

This is also equivalent with doing (as long as we have no arguments)

.. code-block:: python

    >>> wl_kernel = WeisfeilerLehman(niter=5, base_kernel=VertexHistogram)

Now let's go back again to our favorite MUTAG problem.

.. code-block:: python

    >>> from grakel import datasets
    >>> MUTAG = datasets.fetch_dataset("MUTAG", verbose=False)
    >>> MUTAG_data, y = MUTAG.data, MUTAG.target
    >>> split_point = int(len(MUTAG_data) * 0.9)
    >>> X, Y = MUTAG_data[:split_point], MUTAG_data[split_point:]

If what we said till now is correct, the :code:`GraphKernel` object should produce the same kernel matrix output on the MUTAG train/test data split.

.. code-block:: python

    >>> from grakel import GraphKernel
    >>> wl_graph_kernel = GraphKernel(kernel = [{"name": "weisfeiler_lehman", "niter": 5}, {"name": "subtree_wl"}])
    >>> # The alias "subtree_wl" is supported inside the decorator
    >>> from numpy import array_equal
    >>> array_equal(wl_graph_kernel.fit_transform(X), wl_kernel.fit_transform(X))
    True
    >>> array_equal(wl_graph_kernel.transform(Y), wl_kernel.transform(Y))
    True

.. _currying: https://en.wikipedia.org/wiki/Currying
.. _sk learn template: https://github.com/scikit-learn-contrib/project-template

Why not a more structured input for Graphs?
-------------------------------------------
The flattened input type provided for all kernels (graph-dictionary/adjacency, node-labels, edge-labels) may raise the question,
why does not this library, accept a well known type of Graph input as the one constructed from `networkx`_ or `igraph`_.
Networkx library is known for producing a very big memory overhead, which seems unimportant when the user wants to use
very basic graph methods such calculating shortest paths or getting a vertex neighbor. Because what we wanted to wrap
around a graph class was really simple: conversion between dictionary and adjacency formats, format agnostic - format imposing
methods and very basic graph oriented supplementary methods, such as *Shortest-Path matrix* calculation, we designed
a Graph class of our own, used inside most of our kernels, in order to resolve to a common object - graph format reference.
This specificity of kernel format, as well as the absence of a need for complex calculations concerning the field of graphs
lead us to the creation of :ref:`Graph`.

Let's go back to the H2O example:
First we will import the :code:`Graph` object from :code:`Grakel`

.. code-block:: python

    >>> from grakel import Graph

Firstly let's collect all the dictionary formats and show that they are equivalent.
We start by calculating a graph object for the native format of graph dictionary which corresponds to the following:

.. code-block:: python

    >>> H2Od = dict()
    >>> H2Od[0] = {'a': {'b': 1., 'c': 1.}, 'b': {'a': 1}, 'c': {'a': 1}}

Now let's initialize all the other

.. code-block:: python

    >>> H2Od[1] = {'a': ['b', 'c'], 'b': ['a'], 'c':['b']}
    >>> H2Od[2] = {('a', 'b'): 1., ('a', 'c'): 1., ('c', 'a'): 1., ('b', 'a'): 1.}
    >>> H2Od[3] = [('a', 'b'), ('a', 'c'), ('b', 'a'), ('c', 'a')]
    >>> H2Od[4] = [('a', 'b', 1.), ('a', 'c', 1.), ('b', 'a', 1.), ('c', 'a', 1.)]

and compute the result

.. code-block:: python

    >>> any(Graph(H2Od[i]).get_edge_dictionary() == H2Od[0] for i in range(1, 5))
    True

Now let's do the same for adjacency matrix type formats. The numpy array is the native adjacency-matrix format:

.. code-block:: python

    >>> from numpy import array
    >>> H2O = dict()
    >>> H2O[0] = array([[0, 1, 1], [1, 0, 0], [1, 0, 0]])

and with the conversion of other input type formats

.. code-block:: python

    >>> H2O[1] = [[0, 1, 1], [1, 0, 0], [1, 0, 0]]
    >>> from scipy.sparse import csr_matrix
    >>> H2O[2] = csr_matrix(([1, 1, 1, 1], ([0, 0, 1, 2], [1, 2, 0, 0])), shape=(3, 3))

we can demonstrate equality as

.. code-block:: python

    >>> from numpy import array_equal
    >>> all(array_equal(Graph(H2O[i]).get_adjacency_matrix(), H2O[0]) for i in range(1, 3))
    True

Now we would like to initialize two :code:`Graph` type objects one for adjacency_matrix and one for edge_dictionary and show that they are equivalent (using also labels).
First initialize the graph object, created from an adjacency matrix:

.. code-block:: python

    >>> H2O_labels = {0: 'O', 1: 'H', 2: 'H'}
    >>> H2O_edge_labels = {(0, 1): 'pcb', (1, 0): 'pcb', (0, 2): 'pcb', (2, 0): 'pcb'}
    >>> adj_graph = Graph(H2O[0], H2O_labels, H2O_edge_labels, "all")

and one from an edge dictionary:

.. code-block:: python

    >>> H2Od_labels = {'a': 'O', 'b': 'H', 'c': 'H'}
    >>> H2Od_edge_labels = {('a', 'b'): 'pcb', ('b', 'a'): 'pcb', ('a', 'c'): 'pcb', ('c', 'a'): 'pcb'}
    >>> edge_dict_graph = Graph(H2Od[0], H2Od_labels, H2Od_edge_labels, "all")

Firstly we will demonstrate equality of graph type formats:

.. code-block:: python

    >>> array_equal(adj_graph.get_adjacency_matrix(), edge_dict_graph.get_adjacency_matrix())
    True

and

.. code-block:: python

    >>> adj_graph.get_edge_dictionary() == edge_dict_graph.get_edge_dictionary()
    True

and afterwards between labels for :code:`"adjacency"` object formats, defined by the :code:`purpose` argument of the :code:`get_labels` method from the :code:`Graph` type object and for both vertices or edges defined by the :code:`label_type` format of the same method, as

.. code-block:: python

    >>> all((adj_graph.get_labels(purpose="adjacency", label_type=lt), edge_dict_graph.get_labels(purpose="adjacency", label_type=lt)) for lt in ["vertex", "edge"])
    True

Checking equality of the inverse ("edge_dictionary") want hold, because the adjacency matrix, when initialized does not have information about the vertex symbols.
Here we should emphasize that **vertex symbols should be a :code:`sortable` in order for an indexing to be possible**.

.. note::
    When initializing a :code:`Graph` object the 4th argument (named :code:`graph_format`), corresponds to the format the :code:`Graph` will be stored to. The default value of this argument is :code:`"auto"`, which stores the graph in the given format, if it is valid. Explicit format "choices" such as :code:`"adjacency"` or :code:`"dictionary"`, will (covert if needed and) store the :code:`Graph` in this format type. By initializing the :code:`Graph` format as all in the above example, we simply make sure that the :code:`Graph` instance will contain both adjacency and dictionary graph representations and their corresponding edge and adjacency labels for both nodes and edges. Although the methods :code:`get_adjacency_matrix` and `get_edge_dictionary`, construct and return such a graph representation if non existent, the :code:`get_labels` method will change the graph format if the requested labels are not in the desired format and pop a certain warning. If the user wants to avoid doing so he can either set the explicit format afterwards by executing

    .. code-block:: python

        >>> adj_graph = Graph(H2O[0], H2O_labels, H2O_edge_labels)
        >>> adj_graph.set_format("all")

    or declare which is the desired format format he wants the graph to support and it will be included automatically by executing

    .. code-block:: python

        >>> adj_graph.desired_format("dictionary")

    which in that case will set the :code:`Graph` instance format from :code:`"adjacency"` to :code:`"all"`, in order to include the specified format.

After this long introduction of what the :code:`Graph` Object is, the way this can interest the user is by utilizing as input for :code:`GraphKernel`.
Because this Object will act as a mutable-object, any necessary format conversion inside a dataset will happen only ones and the user can execute
multiple kernels on a single dataset with repeating conversions again and again. An important thing to mention here is that a kernel Object **should
not** cause information loss concerning a the :code:`Graph` data Object given as input.

Now let's demonstrate the simple water example on a Shortest-Path kernel, using :code:`Graph` type objects.
First initialize those objects:

.. code-block:: python

    >>> H2O = Graph([[0, 1, 1], [1, 0, 0], [1, 0, 0]], {0: 'O', 1: 'H', 2: 'H'})
    >>> H3O = Graph([[0, 1, 1, 1], [1, 0, 0, 0], [1, 0, 0, 0], [1, 0, 0, 0]], {0: 'O', 1: 'H', 2: 'H', 3:'H'})

And calculate fit transform

.. code-block:: python

    >>> from grakel import GraphKernel
    >>> sp_kernel = GraphKernel(kernel = {"name": "shortest_path"}, normalize=True)
    >>> sp_kernel.fit_transform([H2O])
    1.0

and finally the normalized kernel value, between :math:`\mathbf{H}_{2}\mathbf{O}` and :math:`\mathbf{H}_{3}\mathbf{O}^{+}`

.. code-block:: python

    >>> sp_kernel.transform([H3O])
    0.9428090415820634

which is equivalent with the originally computation, we did on introduction.

.. _networkx: https://networkx.github.io/
.. _igraph: http://igraph.org/python/
