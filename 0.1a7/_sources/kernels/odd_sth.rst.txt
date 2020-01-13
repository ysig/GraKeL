.. _odd_sth:

ODD-STh Kernel
==============
The ODD-STh kernel is a kernel between labeled graphs. Its approach derives from the idea of utilizing tree-based kernels, i.e. kernels that take as input graphs that are trees. Such kernels are in general more computationally efficient as trees are constrained to interesting properties.  
The idea behind the ODD-STh kernel proposed in :cite:`Martino2012ATK`, has to do with decomposition of two graph to ordered DAGs and adding the kernel values between all pairs of DAGs of the original graphs as:

.. math::

    K_{K_{DAG}}(G_{1}, G_{2}) = \sum_{\substack{D_{1} \in DD(G_{1}) \\ D_{2} \in DD(G_{2})}} K_{DAG}(D1, D2)

where :math:`DD(G_{i})` corresponds to a graph decomposition of this graph and :math:`K_{DAG}` is a kernel between DAGs. As a DAG decomposition of each graph they considered the set of all directed BFS explorations starting from each node inside the graph, as follows in the picture:

.. figure:: ../_static/odd_sth_1.png
    :scale: 100 %

    A simple DAG decomposition of a single graph
   
       

Now in order to move from DAGs to trees each :math:`K_{DAG}` kernel was calculated as the sum of tree kernel between derived trees between each of the two DAGs:

.. math::

    K_{DAG} = \sum_{\substack{v_{1} \in V(D_{1}) \\ v_{2} \in V(D_{2})}} C(root(T(v_{1})), root(T(v_{2})))

where :math:`T()` corresponds to the tree-visits on DAGs (which preserve an essence of \textit{ordering} as found in (:cite:`Martino2012ATK`, section 5.2). An example of such tree visits follows:

.. figure:: ../_static/odd_sth_2.png
    :scale: 100 %

    Ordered tree visits on a DAG decomposed from a graph

:math:`C()` is a kernel between trees, where in our case it will be the **S**\ ub-\ **T**\ ree Kernel (as found in :cite:`STKernel`).

.. note::
    Tree isomorphism can `be decided in linear time on the sum of the number of nodes and the number of edges <https://www.geeksforgeeks.org/tree-isomorphism-problem/>`_

For increasing the efficiency of this algorithm for the new set of DAG decomposition, known as ODD (*Ordered Dag Decomposition*), an aggregation of all the decomposition in a single DAG was proposed notated as :math:`BigDAG`. This method introduced in (:cite:`Martino2006`, MinimalDAG: Figure 2, p. 3), aggregates nodes having same labels with frequencies if they correspond to the same path on each DAG, while conserves the existence of nodes that cannot be aggregated.

.. figure:: ../_static/odd_sth_3.png
    :scale: 100 %

    Construction of a :math:`BigDAG` from two DAGs

Doing so allows as to replace the kernel computation:

.. math::

    K_{K_{DAG}}(G_{1}, G_{2}) = \sum_{\substack{D_{1} \in DD(G_{1}) \\ D_{2} \in DD(G_{2})}} K_{DAG}(D1, D2)

with:

.. math::

    K_{BigDAG}(G_{1}, G_{2}) = \sum_{\substack{u_{1} \in V(BigDAG(G_{1}))\\ u_{2} \in V(BigDAG(G_{2})}} f_{u_{1}}f_{u_{2}}C(u_{1}, u_{2})

where :math:`f_{u}` is the frequency counter of the node :math:`u` and :math:`C(u, v)` is the number of matching proper subtrees from :math:`u` and :math:`v`. An even more abstract idea they followed was to created a :math:`Big^{2}DAG` where all the :math:`BigDAGs` created from each graph, would be aggregated to a single one, in the same way as in trees, but instead of incrementing frequencies on common nodes a frequency vector of appended frequencies for each DAG, was constructed.

.. figure:: ../_static/odd_sth_4.png
    :scale: 100 %

    Construction of a :math:`Big^{2}DAG` from two :math:`BigDAGs`

In the final :math:`Big^{2}DAG` graph, the computation of the kernel matrix is all about calculating the following formula:

.. math::

    K_{Big^{2}DAG}(G_{i}, G_{j}) = \sum_{u_{1}, u_{2} \in V(Big^{2}DAG)} F_{u_{1}}[i] * F_{u_{2}}[j] C(u_{1}, u_{2})

which is equivalent to:

.. math::

    K_{Big^{2}DAG}(G_{i}, G_{j}) = \sum_{u \in V(Big^{2}DAG)} F_{u}[i] * F_{u}[j] C(u, u)

because the subtree kernel will have a match only between identical subtrees, that is:

.. math::

    C(u_{1}, u_{2}) \not= 0 \leftrightarrow T(u_{1}) = T(u_{2})

Finally in order to construct the :math:`Big^{2}DAG` each vertex would be represented by a tuple containing a unique hash (whose uniqueness has to do with the ordering) a frequency vector and a depth, which where utilized for calculating the kernel value. In order to restrict the size of the produced graphs a parameter :math:`h` was introduced which restricts the maximum depth of the BFS exploration when doing the graph decomposition.

The Ordered Dag Decomposition - Sub-Tree :math:`h` (ODD-STh) kernel can be found implemented below:

.. currentmodule:: grakel

.. autosummary::

   OddSth

.. note::

    Because the :math:`Big^{2}DAG` graph should be preserved through consequent transformations, the cost of copying it may make :code:`fit_transform` calculation between all the graphs of *train* and *test* faster than fitting on *train* graphs and transforming on *test* graphs.

Bibliography
------------
.. bibliography:: graph_kernels.bib
   :filter: docname in docnames