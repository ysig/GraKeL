.. _introduction:

.. include:: ../.special.rst

====================
A short Introduction
====================

What is grakel?
---------------
GraKeL is a library for the study, use and integration of an upcoming collection
of techniques, inside the field of Machine Learning known as graph kernels. These
techniques utilize information derived from a conceived structure of the data, in
order to apply conventional machine learning techniques for achieving tasks as
classification, ranking, etc. Graph Kernels have been widely used in fields such
as chemistry, bio-informatics, social networks and malware detection and are starting
to be considered as the state-of-the-art solution for various problems inside the field of ML.

What is a Graph Kernel?
-----------------------
A graph kernel is a measure of similarity between two graphs, that obeys a certain
mathematical constraint, which is *that the calculation of each similarity measure between two graphs implies a representation of this two graphs in a* `hilbert space`_, where those to graphs are represented as vectors. This can be notated as the 
following:

    .. math:: 
        k \; : \; \mathcal{G} \times \mathcal{G} \rightarrow \mathbb{R} \text{ is a 
        valid kernel, if there exists a map } \phi \; :\; \mathcal{G} \rightarrow 
        \mathbb{H}, \\\text{for a Hilbert space } \mathbb{H} \text{ where each kernel 
        value can be computed as }\\ k(G_{i}, G_{j}) = \langle G_{i}, G_{j} \rangle
        \text{ where } \langle \;.\; ,\; .\;\rangle \text{ signifies an inner product inside this space}.

The above definition is satisfied in the literature, by proving that the produced kernel matrix from any collection of graphs :math:`\{G_{i}, \text{for } i\in [N]\}`, where :math:`[K]_{ij} = k(G_{i}, G_{j})`, is `Positive Semi Definite`_.

.. _hilbert space: https://en.wikipedia.org/wiki/Hilbert_space
.. _Positive Semi Definite: https://en.wikipedia.org/wiki/Positive-definite_matrix

Initializing a graph kernel
---------------------------
A very well known graph kernel found in literature is the Shortest Path Kernel first introduced by Karsten M. Borgwardt
and Hans-Peter Kriegel on a 2005 article [see :cite:`Borgwardt2005ShortestpathKO`] titled **"Shortest-Path Kernels on Graphs"**,
really essential as an origin of the Graph Kernel field.

After following the instructions found on :ref:`installation`, in order to initialize a *Shortest Path* kernel
using the **grakel** library, you just need to do the following:

.. code-block:: python

    >>> from grakel import GraphKernel
    >>> sp_kernel = GraphKernel(kernel={"name": "shortest_path"})

Kernels as the above are considered as *base kernels*, meaning that they can be computed onto the sets
of Graphs needed only a minor parametrization. A second type of kernels appear in the literature which
we will call *meta-kernels*, which apply transformation operations upon graph objects in order to apply
kernel calculations on certain steps using a *base kernel*, aggregating their result in a certain way.
A kernel like *Weisfeiler-Lehman* introduced by Nino Shervashidze at 2011, published on a journal with the title
"Weisfeiler Lehman Kernels" [see :cite:`Shervashidze2011WeisfeilerLehmanGK`], used a method for approximating a
solution to the graph isomorphism problem, in order to generate a graph refinement scheme that would imply
bigger expressiveness to the base_kernel calculations (an interesting `post`_ explaining the intuition of this kernel).

.. _post: http://blog.smola.org/post/33412570425/the-weisfeiler-lehman-algorithm-and-estimation-on

To initialize such a kernel, using the default subtree kernel, found originally on the paper's
page 9, eq. 2, you can do the following:

.. code-block:: python

    >>> from grakel import GraphKernel
    >>> wl_kernel = GraphKernel(kernel=[{"name": "weisfeiler_lehman"}, {"name": "subtree_wl"}])

Calculate a kernel
------------------
Let's consider a toy example, comparing water :math:`\mathbf{H}_{2}\mathbf{O}` with hydronium
:math:`\mathbf{H}_{3}\mathbf{O}^{+}`, an ion of water produced by protonation.

For start we would calculate the kernel value of water with itself:

.. code-block:: python

    >>> H2O = [[[[0, 1, 1], [1, 0, 0], [1, 0, 0]], {0: 'O', 1: 'H', 2: 'H'}]]
    >>> sp_kernel.fit_transform(H2O)
    array([[12.]])

Now to calculate the graph similarity to hydronium based on the shortest path
graph kernel

.. code-block:: python

    >>> H3O = [[[[0, 1, 1, 1], [1, 0, 0, 0], [1, 0, 0, 0], [1, 0, 0, 0]], {0: 'O', 1: 'H', 2: 'H', 3:'H'}]]
    >>> sp_kernel.transform(H3O)
    array([[24.]])

This result seems like the water molecule is more similar to hydronium, than with itself.
This is a false assumption derived from the fact that the kernel calculation is not normalized.

To apply normalization we add such an argument on the GraphKernel method
initialization and continue

.. code-block:: python

    >>> sp_kernel = GraphKernel(kernel={"name": "shortest_path"}, normalize=True)
    >>> sp_kernel.fit_transform(H2O)
    array([[1.]])
    >>> sp_kernel.transform(H3O)
    array([[0.94280904]])


The input type
--------------
On the above example concerning water and hydronium, we provided a very strange input object
without saying anything about it. The input concerns the user mostly when dealing with an
API, so we will examine it in detail although this is a small introduction.

The input of any kernel - either on the stage of fit or of transform - (to learn more about the
kernel design see :ref:`longer_introduction`) is an iterable, where each element
contains in following order, the next 3 basic elements:

1. The first element is a valid graph object. Valid graph objects can be separated in two major categories (both :red:`weighted` and :blue:`un-weighted`):
    
    * Dictionary representations: This is an *edge* oriented approach, where the input can have one of the following formats:
        - | :red:`2-level nested dictionaries from edge symbols to weights.`  
          | Example: :code:`H2O = {'a': {'b': 1., 'c': 1.}, 'b': {'a': 1}, 'c': {'a': 1.}}`
    
        - | :blue:`Dictionary of symbols to list of symbols.`  
          | Example: :code:`H2O = {'a': ['b', 'c'], 'b': ['a'], 'c': ['b']}`
    
        - | :red:`Dictionary of tuples to weights.`  
          | Example: :code:`H2O = {('a', 'b'): 1., ('a', 'c'): 1., ('c', 'a'): 1., ('b', 'a'): 1.}`
    
        - | :blue:`Iterable of tuples of lenght 2.`  
          | Example: :code:`H2O = [('a', 'b'), ('a', 'c'), ('b', 'a'), ('c', 'a')]`
    
        - | :blue:`Iterable of tuples of length 3.`  
          | Example: :code:`H2O = [('a', 'b', 1.), ('a', 'c', 1.), ('b', 'a', 1.), ('c', 'a', 1.)]`
    
      As seen above all the graph objects are considered **directed** graphs.

    * Array representations: This is a *vertex* oriented approach, where the input can have on of the following formats:
    
        - | :red:`array-like lists of lists`  
          | Example: :code:`H2O = [[0, 1, 1], [1, 0, 0], [1, 0, 0]]`
    
        - | :red:`np.array`  
          | Example: :code:`H2O = numpy.array([[0, 1, 1], [1, 0, 0], [1, 0, 0]])`
    
        - | :red:`sparse matrix (scipy.sparse)`  
          | Example: :code:`H2O = scipy.sparse.csr_matrix(([1, 1, 1, 1], ([0, 0, 1, 2], [1, 2, 0, 0])), shape=(3, 3))`

2. The second *optional* element is a graph labeling of vertices (or nodes):
    
    * | Dictionary representations: Dictionary between vertex symbols and label symbols.  
      | Example: :code:`H2O_labels = {'a': 'O', 'b': 'H', 'c': 'H'}`
    
    * | Array representations: Dictionary between numbers with int keys from :math:`0 \cdots |V|-1` to label symbols.  
      | Example: :code:`H2O_labels = {0: 'O', 1: 'H', 2:'H'}`
    
    .. note::
      Normally in the literature *labels* correspond to scalars or single symbols and not to vector-like objects, which are defined as *attributes*.
      As far as the input representation is concerned the second object is either labels or attributes for graph **vertices** and the distinction
      between attributed kernels or labeled once is specified for each kernel. Here we have made the assumption that never a kernel uses both labels
      or attributes for *vertices* and if it does so, a label representation can be applied such that the kernel can use only the one kind of labels.

3. The third *optional* element is a graph labeling of edges:
    
    * | Dictionary representations: Dictionary between tuples of vertex symbols for all edges and label symbols.  
      | Example: :code:`H2O_edge_labels = {('a', 'b'): 'pcb', ('b', 'a'): 'pcb', ('a', 'c'): 'pcb', ('c', 'a'): 'pcb'}`
    
    * | Array representations: Dictionary between numbers tuples of int keys from :math:`0 \cdots |V|-1` for all matrix entries considered as edges and label symbols.  
      | Example: :code:`H2O_edge_labels = {(0, 1): 'pcb', (1, 0): 'pcb', (0, 2): 'pcb', (2, 0): 'pcb'}`
    
    .. note::
      As soon as the same unification of node labels and attributes is valid, a second distinction should be made here.
      Labels between edges are not weight values. This means that if the user wants to apply such approach to a kernel, that uses weight values in such a way as the
      Random Walk Kernel, she/he should enrich the graph-type input with weights between all edges that correspond to the edge-labels he/she is intended to use.

As defined above the input should be an iterable of any iterable producing at most one and at least three (or more for certain kernels elements).
To signify absence of node labels if the elements produced by each iterable are more than 2 or edge labels if the labels produced are more than 3
the user can provide the empty list or a None Object.

Fitting on a dataset
--------------------
The next important step and final for our short introduction is to see how to apply a kernel on dataset of graphs and labels.

To do so we will utilize the :code:`fetch_dataset` function found on :ref:`datasets`.
Firstly download the dataset:

.. code-block:: python

    >>> from grakel import GraphKernel, datasets
    >>> mutag = datasets.fetch_dataset("MUTAG", verbose=False)
    >>> mutag_data = mutag.data

Note that the :code:`fetch_dataset` function returns a sklearn.utils.Bunch object, where
the graph-data can be found in the data class member of the result.
    
Now let's initialize a Weisfeiler-Lehman Kernel with 5 iterations:

.. code-block:: python

    >>> wl_kernel = GraphKernel(kernel = [{"name": "weisfeiler_lehman", "niter": 5}, {"name": "subtree_wl"}], normalize=True)

Now let's split the dataset in a train/test manner and calculate fit on the train set.

.. code-block:: python

    >>> split_point = int(len(mutag_data) * 0.9)
    >>> X_train, X_test = mutag_data[:split_point], mutag_data[split_point:]
    >>> wl_kernel = GraphKernel(kernel=[{"name": "weisfeiler_lehman", "niter": 5}, {"name": "subtree_wl"}], normalize=True)

In order to apply classification on a dataset based on the calculation of a kernel matrix, one generally needs 
the matrix between all the training data. Namely given :math:`\mathcal{G}^{\text{train}}` a collection of graphs, calculate
the kernel values with a function :math:`\mathcal{K}: \mathcal{G}^{\text{train}} \rightarrow \mathbb{R}^{n_{\text{train}}} \times \mathbb{R}^{n_{\text{train}}}`,
where :math:`n_{\text{train}}` is the number of graphs inside the training set. This function that simply outputs the kernel
matrix between all graphs of the graphs of the training set is equivalent with

.. code-block:: python

    >>> K_train = wl_kernel.fit_transform(X_train)

The :code:`wl_kernel` is now fitted with the train data and we would like given a collection of graphs :math:`\mathcal{G}^{\text{test}}`
to calculate all the kernel values with a function :math:`\mathcal{K}: \mathcal{G}^{\text{train}} \times \mathcal{G}^{\text{test}} \rightarrow \mathbb{R}^{n_{\text{test}}} \times \mathbb{R}^{n_{\text{train}}}` where :math:`n_{\text{test}}` is the number of graphs inside the test set. This function can be calculated as:

.. code-block:: python

    >>> K_test = wl_kernel.transform(X_test)

which is equivalent to calculating:

.. code-block:: python

    >>> K_test = wl_kernel.fit(X_train).transform(X_test)

except the case where the kernel is not deterministic and aside the fact that fitting in the most
cases takes the majority of the overall kernel computation time.

Finally to demonstrate a classification task using a standard SVM, with a precomputed kernel
(a very well known process in the field's literature) we first take the targets (which are the
class labels) as follows:

.. code-block:: python

    >>> y = MUTAG.target
    >>> y_train, y_test = y[:split_point], y[split_point:]

:code:`import` and initialize a sk-learn :code:`SVC`

.. code-block:: python

    >>> from sklearn.svm import SVC
    >>> clf = SVC(kernel='precomputed')

classify

.. code-block:: python

    >>> clf.fit(K_train, y_train)
    SVC(C=1.0, cache_size=200, class_weight=None, coef0=0.0,
	  decision_function_shape='ovr', degree=3, gamma='auto',
	  kernel='precomputed', max_iter=-1, probability=False, random_state=None,
	  shrinking=True, tol=0.001, verbose=False)
    >>> y_pred = clf.predict(K_test)

and print the accuracy score

.. code-block:: python

    >>> from sklearn.metrics import accuracy_score
    >>> print("%2.2f %%" %(round(accuracy_score(y_test, y_pred)*100)))
    78.95 %
