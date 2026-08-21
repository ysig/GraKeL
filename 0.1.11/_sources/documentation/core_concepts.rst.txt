.. _core_concepts:

=============
Core Concepts
=============

*GraKeL* is organized as shown in the Figure below.

.. figure:: ../_figures/grakel_schema.svg
  :scale: 100%
  :align: center

Specifically, the package consists of two subpackages:

- :class:`grakel.datasets`
- :class:`grakel.kernels`

The :class:`grakel.datasets` subpackage contains functions for automatically downloading graph datasets. On the other hand, the :class:`grakel.kernels` subpackage contains the :class:`grakel.Kernel` class, and the implementations of all graph kernels and frameworks. Note that all graph kernels and frameworks such as the :class:`grakel.ShortestPath` kernel or the :class:`grakel.WeisfeilerLehman` framework inherit from the :class:`grakel.Kernel` class.

GraKeL also contains the following two classes:

- :class:`grakel.Graph`
- :class:`grakel.GraphKernel`

The :class:`grakel.Graph` class is used to represent graphs and also provides functions for manipulating them. The :class:`grakel.GraphKernel` class is a generic wrapper. This class provides a uniform interface for all the implemented graph kernels and frameworks.

We next present some core components of *GraKeL*.

What is the :class:`grakel.GraphKernel` Class?
----------------------------------------------
The :class:`grakel.GraphKernel` class is a *generic wrapper class*. This class provides a uniform interface for all the implemented graph kernels and frameworks. A graph kernel can be described by an instance of this class, and it holds the attributes listed below:

* :code:`kernel` : Specifies the graph kernel to be computed. It can be either a :code:`base_graph_kernel` or a list that contains one or more :code:`framework` along with exactly one :code:`base_graph_kernel`. The :code:`base_graph_kernel` needs to be the last element in the list.
    - :code:`base_graph_kernel` : Α :code:`base_graph_kernel` is a kernel that compares graphs to each other. It is represented by a dictionary which contains a key :code:`'name'` whose value  corresponds to the name of the kernel. The dictionary can also contain other keys that specify the parameters of the kernel and their values. For instance, we can initialize a shortest path kernel as follows.

    .. code-block:: python

        >>> from grakel import GraphKernel
        >>> gk = GraphKernel(kernel={"name": "shortest_path", "with_labels": False})

    - :code:`framework` : A :code:`framework` works on top of graph kernels. It takes a :code:`base_graph_kernel` as input. Frameworks correspond to dictionaries that contain their name as the value of the key :code:`'name'`, and their parameters. A :code:`framework` combined with a :code:`base_graph_kernel` corresponds to a :code:`base_graph_kernel` and can be passed on to another :code:`framework`. For example, a kernel that applies the Weisfeiler-Lehman framework on top of the shortest path kernel is initialized as follows.

    .. code-block:: python

        >>> from grakel import GraphKernel
        >>> gk = GraphKernel(kernel=[{"name": "weisfeiler_lehman", "n_iter": 5}, {"name": "shortest_path"}])

* :code:`normalize` : A kernel can provide either an unnormalized or a normalized output.
    The normalized kernel value between two graphs :math:`G_1` and :math:`G_2` is computed as follows: :math:`k(G_1, G_2)/\sqrt{k(G_1, G_1) k(G_2, G_2)}`. This normalization ensures that the kernel value between a graph and itself is equal to 1, while the kernel value between a graph and any other graph takes values between 0 and 1.

    | **Example**
    
    Suppose we have a set of training graphs :code:`G_train`, and a set of test graphs :code:`G_test`. We compute the normalized kernel matrices using the Weisfeiler-Lehman subtree kernel as follows.

    .. code-block:: python

        >>> gk = GraphKernel(kernel=[{"name": "weisfeiler_lehman", "n_iter": 5}, {"name": "subtree_wl"}], normalize=True)
        >>> # Calculate the normalized kernel matrices
        >>> K_train = gk.fit_transform(G_train)
        >>> K_test = gk.transform(G_test)

    The above is equivalent (for deterministic kernels) to the code below.

    .. code-block:: python

        >>> gk = GraphKernel(kernel=[{"name": "weisfeiler_lehman", "n_iter": 5}, {"name": "subtree_wl"}], normalize=False)
        >>> K = gk.fit_transform(G)
        >>> K_diag = K.diagonal()
        >>> K_train_diag, K_test_diag = K_diag[idx_train], K_diag[idx_test]
        >>> # Calculate the normalized kernel matrices
        >>> K_train = K[idx_train, :][:, idx_train] / np.sqrt(np.outer(K_train_diag, K_train_diag))
        >>> K_test = K[idx_test, :][:, idx_train] / np.sqrt(np.outer(K_test_diag, K_train_diag))

    Note that in the second case, we perform more computations since we also compare the graphs of the test set to each other.

* :code:`Nystroem` : The Nyström method is a well-established approach for approximating kernel matrices on large datasets.
    If :math:`n` is the number of samples, computing and storing the kernel matrix requires :math:`\mathcal{O}(n^2)` time and memory, respectively. Therefore, applying kernel methods will become unfeasible when :math:`n` is large. The Nyström approximation can allow a significant speed-up of the calculations by computing an approximation :math:`\tilde{\mathbf{K}}` of rank :math:`q` of the kernel matrix. The method uses a subset of the training data as basis and reduces the storage and complexity requirements to :math:`\mathcal{O}(n q)`. The value of :math:`q` is specified by the user by setting :code:`Nystroem` equal to an integer value. An example demonstrating the power of the Nyström method is given below.

    | **Example**

    We first download the MUTAG dataset and split it into a training and a test set.

    .. doctest:: 

        >>> from grakel.datasets import fetch_dataset
        >>> from sklearn.model_selection import train_test_split
        >>> MUTAG = fetch_dataset("MUTAG", verbose=False)
        >>> G = MUTAG.data
        >>> y = MUTAG.target
        >>> G_train, G_test, y_train, y_test = train_test_split(G, y, test_size=0.1, random_state=42)

    We next initialize a Weisfeiler-Lehman subtree kernel using :code:`GraphKernel`, and we also make use of :code:`Nystroem` with :math:`q=20` to approximate the kernel matrix.

    .. doctest:: 

        >>> from grakel import GraphKernel
        >>> gk = GraphKernel(kernel=[{"name": "weisfeiler_lehman", "n_iter": 5}, "subtree_wl"], Nystroem=20)
        >>> K_train = gk.fit_transform(G_train)
        >>> K_test = gk.transform(G_test)
        >>> print(K_train.shape)
        (169, 20)
        >>> print(K_test.shape)
        (19, 20)

    Then, we train a standard SVM classifier with linear kernel, and use the classifier to make predictions.

    .. doctest:: 

        >>> from sklearn.svm import SVC
        >>> clf = SVC(kernel='linear')
        >>> clf.fit(K_train, y_train)
        SVC(kernel='linear')
        >>> y_pred = clf.predict(K_test)

    Finally, we calculate the classification accuracy.

    .. code-block:: python

        >>> from sklearn.metrics import accuracy_score
        >>> print(str(round(accuracy_score(y_test, y_pred)*100, 2)), "%")
        94.74 %

    .. note::
        | To compute the full kernel matrices, we needed to perform :math:`~ 169 * (169-1) /2 + 19 * 169 = 17,407` kernel computations. Instead, we performed :math:`~ 20 * (20-1)/ 2 + 20 * 169 + 20* 19 = 3,950` kernel computations. As we can see, the approximation also led to an increase in performance.

* :code:`n_jobs` : Some kernels consist of operations that can be executed in parallel, leading to a reduction in the running time.
    The :code:`n_jobs` attribute has the same functionality as that of scikit-learn. It determines the number of jobs that will run in parallel. If :code:`n_jobs` is set equal to -1, all the processors will be utilized. Note that this attribute will not have an impact on the computation of some kernels whose code is not parallelized. These kernels either take advantage of the parallelization inherent in other libraries (e.g., NumPy) or their code is only partially parallelizable or not parallelizable at all. In such scenarios, a warning is issued.

    If you are interested in parallelizing any of the implemented kernels, you can *contribute* to the *GraKeL* project. To find out how you can contribute, please have a look at :ref:`contributing`.

* :code:`random_state` : This attribute is used for initializing the internal random number generator.
    It has no effect on deterministic graph kernels, but only on kernels that involve some random process (e.g., those that perform sampling). It also applies to the :code:`Nystroem` function of the :code:`GraphKernel` class which also performs sampling. If int, :code:`random_state` is the seed used by the random number generator. Otherwise, it can be a :code:`RandomState` instance. If :code:`None`, the random number generator is the :code:`RandomState` instance used by :code:`np.random`. The use of :code:`random_state` is illustrated in the following example.

    | **Example**

    We first create the graph representations of the following two molecules: (1) water :math:`\mathbf{H}_{2}\mathbf{O}` and (2) hydronium :math:`\mathbf{H}_{3}\mathbf{O}^{+}`, an ion of water produced by protonation.

    .. doctest::

       >>> from grakel import Graph
       >>>
       >>> H2O_adjacency = [[0, 1, 1], [1, 0, 0], [1, 0, 0]]
       >>> H2O_node_labels = {0: 'O', 1: 'H', 2: 'H'}
       >>> H2O = Graph(initialization_object=H2O_adjacency, node_labels=H2O_node_labels)
       >>>
       >>> H3O_adjacency = [[0, 1, 1, 1], [1, 0, 0, 0], [1, 0, 0, 0], [1, 0, 0, 0]]
       >>> H3O_node_labels = {0: 'O', 1: 'H', 2: 'H', 3:'H'}
       >>> H3O = Graph(initialization_object=H3O_adjacency, node_labels=H3O_node_labels)

    We will then compute the *graphlet kernel* between the two molecules. The graphlet kernel counts the number of common graphlets (i.e., small subgraphs) in two graphs. Instead of exaustively enumerating all the graphlets, it usually samples a number of them. In this example, we will sample 5 graphlets from each graph.

    .. doctest::

        >>> from grakel import GraphKernel
        >>> gk = GraphKernel(kernel=dict(name="graphlet_sampling", sampling=dict(n_samples=5)))
        >>> gk.fit([H2O])
        GraphKernel(kernel={'name': 'graphlet_sampling', 'sampling': {'n_samples': 5}})    

    .. code-block:: python

        >>> gk.transform([H3O])
        array([[10.]])

    Note that we did not set :code:`random_state` to some value, and therefore it took its default :code:`None` value. We will now set :code:`random_state` equal to 42.

    .. doctest:: 

        >>> gk = GraphKernel(kernel=dict(name="graphlet_sampling", sampling=dict(n_samples=5), random_state=20))
        >>> gk.fit([H2O])
        GraphKernel(kernel={'name': 'graphlet_sampling', 'random_state': 20,
                            'sampling': {'n_samples': 5}})

        >>> gk.transform([H3O])
        array([[20.]])

    As you can see, the new kernel value is not equal to the previous one. If we re-run the above code, we will obtain the same kernel value since the algorithm will sample exactly the same graphlets from both graphs. As shown below, we can also obtain the same kernel value if :code:`random_state` is initialized as an attribute of :code:`GraphKernel` instead of the graphlet kernel itself.

    .. doctest::

        >>> gk = GraphKernel(kernel=dict(name="graphlet_sampling", sampling=dict(n_samples=5)), random_state=20)
        >>> gk.fit([H2O])
        GraphKernel(kernel={'name': 'graphlet_sampling', 'sampling': {'n_samples': 5}},
                    random_state=20)
    
        >>> gk.transform([H3O])
        array([[20.]])

    If we provide a :code:`random_state` value to both :code:`GraphKernel` and :code:`kernel`, then each one will have an effect only on the corresponding instances.

    .. doctest::

        >>> gk = GraphKernel(kernel=dict(name="graphlet_sampling", sampling=dict(n_samples=5, random_state=0)), random_state=20)
        >>> gk.fit([H2O])
        GraphKernel(kernel={'name': 'graphlet_sampling',
                            'sampling': {'n_samples': 5, 'random_state': 0}},
                    random_state=20)
    
        >>> gk.transform([H3O])
        array([[20.]])

    while

    .. doctest::

        >>> gk = GraphKernel(kernel=dict(name="graphlet_sampling", sampling=dict(n_samples=5)), random_state=0)
        >>> gk.fit([H2O]).transform([H3O])
        array([[10.]])


* :code:`verbose` : Currently not supported.
    .. note::
        :code:`verbose` is an attribute that is currently not supported, but may be supported in the future for printing progress messages.

We will next focus on the :class:`grakel.Kernel` class. Instances of this class are wrapped in an instance of the :class:`grakel.GraphKernel` class that was presented above.

The :class:`grakel.Kernel` class
--------------------------------
All graph kernels inherit from this class.

A graph kernel is a function :math:`k` between two graphs. That is, :math:`k \; : \; \mathcal{G} \times \mathcal{G} \rightarrow \mathbb{R}` where :math:`\mathcal{G}` is the space of graphs. We usually do not have just two graphs, but a large set of graphs, and we are interested to compare these graphs to each other using some kernel. In almost all cases, it is more computationally efficient to compute all the kernel values in one step than computing the kernel value for each pair individaully. Therefore, we designed our kernels to take sets of graphs as input instead of just two graphs.

The *GraKeL* package had also to be compatible with *scikit-learn*. From the different scikit-learn structures, the one that fitted best to our setting was the :code:`TransformerMixin` class, which consists of the following three methods: :code:`fit`, :code:`fit_transform` and :code:`transform`. The three methods are designed to perform the following tasks in our package:

- The :code:`fit` method extracts kernel dependent features from an input graph collection.

- The :code:`fit_transform` method does the same job as :code:`fit`, but also computes the kernel matrix emerging from the input graph collection.

- The :code:`transform` method calculates the kernel matrix between a new collection of graphs and the one given as input to :code:`fit` or to :code:`fit_transform`.

.. note::
    The :code:`fit` and :code:`fit_transform` methods usually extract some features from the set of graphs that is given as input. These features are stored into memory and are not modified by the applications of the :code:`transform` method. This (the need to copy and protect the extracted data) however adds some overhead to the computation of some kernels such as the ODD-STh kernel. In such cases, the user may prefer to use the :code:`fit_transform` method once and then manually retrieve the two kernel matrices.

The Figure below illustrates how the :class:`grakel.Kernel` class is organized.

.. figure:: ../_figures/kernel_schema.svg
  :scale: 100%
  :align: center

Besides the three methods discussed above, there also exist some other methods as shown in the Figure. As can be seen, these methods are called by the :code:`fit`, :code:`fit_transform` and :code:`transform` methods. The :code:`diagonal` method is used for normalizing kernel matrices. It returns the self-kernel values of all the graphs given as input to :code:`fit` along with those given as input to :code:`transform`, provided that this method has been called. The :code:`parse_input` method extracts features from the collection of graphs that is given as input either to :code:`fit` or to :code:`transform`. The :code:`pairwise_operation` method computes the kernel between two graphs. This method is used by the :code:`calculate_kernel_matrix` method which generates kernel matrices from collections of graphs. Finally, the :code:`initialize_` function is used just for initialization purposes.


A kernel initialized as an instance of the :class:`grakel.Kernel` class is equivalent to an instance of the :class:`grakel.GraphKernel` generic wrapper corresponding to the same kernel if the attributes of the two kernels are identical to each other. To illustrate this, we will employ a deterministic graph kernel (the Wesfeiler-Lehman subtree kernel) and we will investigate if the kernel values produced by the two instances of the kernel are equal to each other.

We first initialize the instance of the :class:`grakel.Kernel` class. This corresponds to the Weisfeiler-Lehman framework on top of the vertex histogram kernel.

.. code-block:: python

    >>> from grakel import WeisfeilerLehman, VertexHistogram
    >>> gk_1 = WeisfeilerLehman(n_iter=5, base_graph_kernel=VertexHistogram)

We have set the :code:`base_graph_kernel` attribute equal to the :class:`grakel.kernels.VertexHistogram` class. Note that the :code:`base_graph_kernel` attribute can also be set equal to a tuple consisting of a :class:`grakel.kernel` class and a dictionary containing the attributes of the corresponding kernel and their values. Above, we have set the attributes of the vertex histogram kernel to their default values. Therefore, the above code is equivalent to the following.

.. code-block:: python

    >>> gk_1 = WeisfeilerLehman(n_iter=5, base_graph_kernel=(VertexHistogram, {}))

We will perform our experiment on the MUTAG dataset.

.. code-block:: python

    >>> from grakel.datasets import fetch_dataset
    >>> MUTAG = fetch_dataset("MUTAG", verbose=False)
    >>> G = MUTAG.data
    >>> y = MUTAG.target
    >>> K_1 = gk_1.fit_transform(G)

We will now test if the kernel matrix produced by the instance of the :class:`grakel.GraphKernel` class is equal to the one produced by the instance of the :class:`grakel.Kernel` class.

.. code-block:: python

    >>> from grakel import GraphKernel
    >>> from numpy import array_equal
    >>> gk_2 = GraphKernel(kernel = [{"name": "weisfeiler_lehman", "n_iter": 5}, {"name": "subtree_wl"}]) # The alias "subtree_wl" is supported inside the generic wrapper
    >>> K_2 = gk_2.fit_transform(G)
    >>> array_equal(K_1, K_2)
    True

As we can see, the two matrices are indeed equal to each other.

Why Not a More Advanced Graph Representation?
---------------------------------------------
As already discussed, the graph objects in *GraKeL* are instances of the :class:`grakel.Graph` class. The :class:`grakel.Graph` class is very simple, and this may raise the question why *GraKeL* does not utilize the graph structures of well-established graph libraries such as `networkx`_ and `igraph`_. The answer is that the operations that most kernels perform on graphs are relatively simple and easily implementable. For instance, a kernel may need to retrieve the neighbors of a vertex or to compute the shortest paths between all pairs of nodes. Standard graph libraries provide many more functions, and they are specially designed such that all these functions are computed efficiently. Since *GraKeL* would only utilize a small fraction of these functions, introducing an extra dependency to some large library seemed not to be a good idea.

We will again experiment with the two molecules: (1) water :math:`\mathbf{H}_{2}\mathbf{O}` and (2) hydronium :math:`\mathbf{H}_{3}\mathbf{O}^{+}`.

We will first initialize five water molecules using the different edgelist representations and show that they are equivalent to each other.

.. code-block:: python

    >>> from grakel import Graph
    >>> H2Od = list()
    >>> H2Od.append(Graph({'a': {'b': 1., 'c': 1.}, 'b': {'a': 1}, 'c': {'a': 1}}))
    >>> H2Od.append(Graph({'a': ['b', 'c'], 'b': ['a'], 'c':['b']}))
    >>> H2Od.append(Graph({('a', 'b'): 1., ('a', 'c'): 1., ('c', 'a'): 1., ('b', 'a'): 1.}))
    >>> H2Od.append(Graph([('a', 'b'), ('a', 'c'), ('b', 'a'), ('c', 'a')]))
    >>> H2Od.append(Graph([('a', 'b', 1.), ('a', 'c', 1.), ('b', 'a', 1.), ('c', 'a', 1.)]))

Then, we compare the first representation against all the other.

.. code-block:: python

    >>> any(H2Od[i].get_edge_dictionary() == H2Od[0].get_edge_dictionary() for i in range(1, 5))
    True

Now, we will do the same for the case of the adjacency matrix representations.

.. code-block:: python

    >>> import numpy as np
    >>> from scipy.sparse import csr_matrix
    >>> H2O = list()
    >>> H2O.append(Graph(np.array([[0, 1, 1], [1, 0, 0], [1, 0, 0]])))
    >>> H2O.append(Graph([[0, 1, 1], [1, 0, 0], [1, 0, 0]]))
    >>> H2O.append(Graph(csr_matrix(([1, 1, 1, 1], ([0, 0, 1, 2], [1, 2, 0, 0])), shape=(3, 3))))

Then, we again compare the first representation against all the other.

.. code-block:: python

    >>> from numpy import array_equal
    >>> all(array_equal(H2O[i].get_adjacency_matrix(), H2O[0].get_adjacency_matrix()) for i in range(1, 3))
    True

Next, we will create two instances of the :code:`grakel.Graph` class, the first using the adjacency_matrix representation and the second using the edgelist representation. We will also assign labels to the nodes and edges of the two graphs. Then, we will show that the two representations are equivalent to each other.

We create the adjacency matrix and use this matrix to create the first object.

.. code-block:: python

    >>> H2O_adj = np.array([[0, 1, 1], [1, 0, 0], [1, 0, 0]])
    >>> H2O_labels = {0: 'O', 1: 'H', 2: 'H'}
    >>> H2O_edge_labels = {(0, 1): 'pcb', (1, 0): 'pcb', (0, 2): 'pcb', (2, 0): 'pcb'}
    >>> adj_graph = Graph(H2O_adj, H2O_labels, H2O_edge_labels, "all")

We then create the second graph object.

.. code-block:: python

    >>> H2Od_edg = {'a': {'b': 1., 'c': 1.}, 'b': {'a': 1}, 'c': {'a': 1}}
    >>> H2Od_labels = {'a': 'O', 'b': 'H', 'c': 'H'}
    >>> H2Od_edge_labels = {('a', 'b'): 'pcb', ('b', 'a'): 'pcb', ('a', 'c'): 'pcb', ('c', 'a'): 'pcb'}
    >>> edge_dict_graph = Graph(H2Od_edg, H2Od_labels, H2Od_edge_labels, "all")

We test if the adjacency matrices of the two objects are equal to each other.

.. code-block:: python

    >>> array_equal(adj_graph.get_adjacency_matrix(), edge_dict_graph.get_adjacency_matrix())
    True

and

.. code-block:: python

    >>> adj_graph.get_edge_dictionary() == edge_dict_graph.get_edge_dictionary()
    True

Finally, we also compare the labels of the nodes and the edges of the two objects.

.. code-block:: python

    >>> all((adj_graph.get_labels(purpose="adjacency", label_type=lt), edge_dict_graph.get_labels(purpose="adjacency", label_type=lt)) for lt in ["vertex", "edge"])
    True

Above, we showed that the adjacency matrices of the two objects are equal to each other. The same does not hold for their edge dictionaries (i.e., :code:`edge_dictionary`) since the adjacency matrix contains no information about the names of the nodes. Note that these names have to be instances of some **sortable** datatype such that indexing can be performed.

.. note::
    The fourth attribute of the constructor of the :code:`grakel.Graph` class (i.e., :code:`graph_format`) corresponds to the format into which the graph object will be stored. The default value of this attribute is :code:`"auto"` which maintains the format that is passed on to the constructor. This attribute can also take the values :code:`"adjacency"`, :code:`"dictionary"`, and :code:`all`. The last value ensures that the :code:`grakel.Graph` instance will contain both representations and their corresponding node and edge labels. Note that the :code:`get_adjacency_matrix` and :code:`get_edge_dictionary` methods create and return the corresponding graph representation if it does not exist. On the other hand, the :code:`get_labels` method will modify the graph format if the labels are not in the proper format and a warning will also be issued. Note that the user can set the :code:`graph_format` attribute to some value later on as follows.

    .. code-block:: python

        >>> adj_graph = Graph(H2O_adj, H2O_labels, H2O_edge_labels)
        >>> adj_graph.change_format("all")

    Alternatively, the user can specify which is his/her desired format, and it will be created if it does not exist.

    .. code-block:: python

        >>> adj_graph.desired_format("dictionary")

The methods of the graph kernels take lists of :class:`grakel.Graph` objects as input, extract the necessary features and may return some matrices. It should be mentioned that the :class:`grakel.Kernel` objects are not allowed to modify the graphs that they take as input.

.. _networkx: https://networkx.github.io/
.. _igraph: http://igraph.org/python/
