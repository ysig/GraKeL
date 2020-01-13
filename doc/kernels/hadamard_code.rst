.. _hadamard_code:

Hadamard Code Kernel
====================

A similar framework to the neighborhood-hashing kernel and the Weisfeiler-Lehman kernel was introduced by Tetsuya Kataoka and Akihito Inokuchi in :cite:`icpram16`, known as Hadamard-code kernel.
Given a collection of **labeled** graphs :math:`\mathbf{G}=[G]^{N}_{i=1}` collect the set :math:`\Sigma` of all distinct labels inside :math:`\mathbf{G}`. A :math:`2^{k}`-th Hadamard code matrix :math:`H_{2^{k}}` is defined as follows:

.. math::

    H_{2^{k+1}}= \begin{cases}
    \begin{pmatrix}
        1 & 1\\
        1 & -1
    \end{pmatrix},\text{ if }k = 0
    \\\\
    \begin{pmatrix}
        H_{2^{k}} & H_{2^{k}}\\
        H_{2^{k}} & -H_{2^{k}}
    \end{pmatrix},\text{if } k > 0
    \end{cases}

Now by defining a Hadamard matrix :math:`\mathbb{H} = H_{2^{\lceil \log_{2}|\Sigma|\rceil}}`, then we initially label each node inside a graph:

.. math::

    l^{(0)}(v) = \mathtt{row}_{i}\mathbb{H},\text{ }\textbf{iff}\text{ }label(v) = \Sigma_{i}

Based on this initial labeling the following relabeling rule:

.. math::

    l^{(k+1)}(v) = l^{(k)}(v) + \sum_{u \in N(v)}l^{(k)}(u)

was used. :math:`N(v)` is used to denote the neighborhood of a node :math:`v`.
Following the above scheme, relabeling is applied iteratively for a fixed number of iterations, while each kernel matrix (calculated from a give *base-kernel*) between relabeled graphs is aggregated to a total one, through summation.

.. figure:: ../_static/kataoka1.png

    An example of the relabeling procedure of the Hadamard code kernel for a single graph.

    

The implementation of the hadamard code kernel framework can be found below. Note that use can use :code:`base_kernel` to attach as a base kernel any kernel for **labeled** graphs.

.. currentmodule:: grakel

.. autosummary::

   HadamardCode
   
Bibliography
------------
.. bibliography:: graph_kernels.bib
   :filter: docname in docnames
