.. _neighborhood_hash:

Neighborhood Hash Kernel
========================


The neighborhood hash kernel assumes node-labeled graphs :cite:`Hido2009ALG`.
It compares graphs by updating their node labels and counting the number of common labels.
The kernel replaces the discrete node labels with binary arrays of fixed length, and it then employs logical operations to update the labels so that they contain information regarding the neighborhood structure of each vertex.

Let :math:`\ell : \mathcal{V} \rightarrow \Sigma` be a function that maps vertices to an alphabet :math:`\Sigma` which is the set of possible discrete node labels.
Hence, given a vertex :math:`v`, :math:`\ell(v) \in \Sigma` is the label of vertex :math:`v`.
The algorithm first transform each discrete node label to a bit label.
A bit label is a binary array consisting of :math:`d` bits as

.. math::

    s = \{ b_1, b_2, \ldots, b_d \}

where the constant :math:`d` satisfies :math:`2^d - 1 \gg |\Sigma|` and :math:`b_1, b_2, \ldots, b_d \in \{0, 1\}`.

The most important step of the algorithm involves a procedure that updates the labels of the vertices.
To achieve that, the kernel makes use of two very common bit operations: (1) the exclusive or (:math:`XOR`) operation, and (2) the bit rotation (:math:`ROT`) operation.
Let :math:`XOR(s_i, s_j) = s_i \oplus s_j` denote the :math:`XOR` operation between two bit labels :math:`s_i` and :math:`s_j` (\ie the :math:`XOR` operation is applied to all their components).
The output of the operation is a new binary array whose components represent the :math:`XOR` value between the corresponding components of the :math:`s_i` and :math:`s_j` arrays.
The :math:`ROT_o` operation takes as input a bit array and shifts its last :math:`o` bits to the left by :math:`o` bits and moves the first :math:`o` bits to the right end as shown below  

.. math::

    ROT_o(s) = \{ b_{o+1}, b_{o+2}, \ldots, b_d, b_1, \ldots, b_o \}


Below, we present in detail two procedures for updating the labels of the vertices: (1) the simple neighborhood hash, and (2) the count-sensitive neighborhood hash.

Simple Neighborhood Hash
------------------------

Given a graph :math:`G=(V,E)` with bit labels, the simple neighborhood hash update procedure computes a neighborhood hash for each vertex using the logical operations :math:`XOR` and :math:`ROT`.
More specifically, given a vertex :math:`v \in V`, let :math:`\mathcal{N}(v)=\{ u_1,\ldots,u_d \}` be the set of neighbors of :math:`v`.
Then, the kernel computes the neighborhood hash as

.. math::

    NH(v) = ROT_1 \big( \ell(v) \big) \oplus \big( \ell(u_1) \oplus \ldots \oplus \ell(u_d) \big)

The resulting hash :math:`NH(v)` is still a bit array of length :math:`d`, and we regard it as the new label of :math:`v`.
This new label represents the distribution of the node labels around :math:`v`.
Hence, if :math:`v_i` and :math:`v_j` are two vertices that have the same label (\ie :math:`\ell(v_i) = \ell(v_j)`) and the label sets of their neighborhors are also identical, their hash values will be the same (\ie :math:`NH(v_i) = NH(v_j))`.
Otherwise, they will be different except for accidental hash collisions.
The main idea behind this update procedure is that the hash value is independent of the order of the neighborhood values due to the properties of the :math:`XOR` operation.
Hence, one can check whether or not the distributions of neighborhood labels of two vertices are equivalent without sorting or matching these two label sets.

Count-sensitive Neighborhood Hash
---------------------------------

The simple neighborhood hash update procedure described above suffers from some problematic hash collisions.
Specifically, the neighborhood hash values for two independent nodes have a small probability of being the same even if there is no accidental hash collision.
Such problematic hash collisions may affect the positive semidefiniteness of the kernel.
To address that problem, the count-sensitive neighborhood hash update procedure counts the number of occurences of each label in the set.
More specifically, it first uses a sorting algorithm (\eg radix sort) to align the bit labels of the neighbors, and then, it extracts the unique labels (set :math:`\{ \ell_1, \ldots, \ell_l \}` in the case of :math:`l` unique labels) and for each label counts its number of occurences.
Then, it updates each unique label based on its number of occurences as follows

.. math::

    \ell'_i = ROT_o \big( \ell_i \oplus o \big)

where :math:`\ell_i, \ell'_i` is the initial and updated label respectively, and :math:`o` is the number of occurences of that label in the set of neighbors.
The above operation makes the hash values unique by depending on the number of label occurrences.
Then, the count-sensitive neighborhood hash is computed as

.. math::

    CSNH(v) = ROT_1 \big( \ell(v) \big) \oplus \big( \ell'_1 \oplus \ldots \oplus \ell'_l \big)

Both the simple and the count-sensitive neighborhood hash can be seen as general approaches for enriching the labels of vertices based on the label distribution of their neighborhood vertices.


Kernel Calculation
------------------

The neighborhood hash update procedures presented above aggregate the information of the neighborhood vertices to each vertex.
Then, given two graphs :math:`G` and :math:`G'`, the updated labels of their vertices are compared using the following function

.. math::

    \kappa(G, G') = \frac{c}{|V| + |V'| - c}

where :math:`c` is the number of labels the two graphs have in common.
This function is equivalent to the Tanimoto coefficent which is commonly used as a similarity measure between sets of discrete values and which has been proven to be positive semidefinite :cite:`gower1971general`.

The label-update procedures is not necessary to be applied once, but they can be applied iteratively.
By updating the bit labels several times, the new labels can capture high-order relationships between vertices.
For instance, if the procedure is performed :math:`h` times in total, the updated label :math:`\ell(v)` of a vertex :math:`v` represents the label distribution of its :math:`h`-neighbors.
Hence, two vertices :math:`v_i, v_j` with identical labels and connections among their :math:`r`-neighbors will be assigned the same label.
Given a graph :math:`G=(V,E)`, let :math:`G_1, \ldots, G_h` denote the updated graphs where the node labels have been updated :math:`1,\ldots,h` times, respectively.
Then, given two graphs :math:`G` and :math:`G'`, the neighborhood hash kernel is defined as

.. math::

    k(G, G') = \frac{1}{h} \sum_{i=1}^h \kappa(G_i, G'_i)

The computational complexity of the neighborhood hash kernel is :math:`\mathcal{O}(dhn\bar{D})` where :math:`n=|V|` is the number of vertices of the graphs and :math:`\bar{D}` is the average degree of their vertices.

The implementation of the neighborhood hash kernel can be found below

.. currentmodule:: grakel

.. autosummary::

   NeighborhoodHash

Bibliography
------------
.. bibliography:: graph_kernels.bib
   :filter: docname in docnames
