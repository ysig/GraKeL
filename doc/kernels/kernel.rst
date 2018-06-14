.. _kernel:

Kernel (general class)
======================

In the literature, a graph kernel appears as a function: :math:`k \; : \; \mathcal{G} \times \mathcal{G} \rightarrow \mathbb{R}` for which there exists a map:

.. math::

    \phi \; :\; \mathcal{G} \rightarrow \mathbb{H}, \text{for a Hilbert space } \mathbb{H}

where each kernel value can be computed as :math:`k(G_{i}, G_{j}) = \langle \phi(G_{i}), \phi(G_{j}) \rangle` where :math:`\langle \;.\; ,\; .\;\rangle` signifies an inner product inside this space. The matrix :math:`\mathbf{K}_{ij} = k(G_{i}, G_{j})` produced is known as the kernel matrix. For a kernel matrix to be valid it is required to be *positive semi definite* (i.e. :math:`\lambda_{min} \ge 0`).

.. note::

   The kernels implemented inside this package have all a **Polynomial** time complexity.

What we found out was that if instead of designing the kernel to be calculated pairwise, we designed it for a collection of graphs :math:`[G]_{i=1}^{N}` we have a significant computational advantage. Given two collections of graphs: :math:`G^{n}, G^{m}`, the full kernel matrix :math:`\mathcal{K}`, should look like:

.. math::

    \mathcal{K} =
    \left[
    \begin{array}{c||c}
    \mathcal{K}^{n\times n} & \mathcal{K}^{n\times m} \\
    \hline
    \hline
    \mathcal{K}^{m\times n} & \mathcal{K}^{m\times m}
    \end{array}
    \right]

Any :code:`Kernel` object inherited class should match the following behavior:

- :math:`\mathcal{K}^{n\times n}=\texttt{<KernelName>.fit\_transform}(\mathcal{G}^{n})`
- :math:`\mathcal{K}^{m\times n}=\texttt{<KernelName>.fit}(\mathcal{G}^{\text{n}}).\texttt{transform}(\mathcal{G}^{\text{m}})`
- :math:`\mathcal{K}=\texttt{<KernelName>.fit\_transform}([\mathcal{G}^{n}\; \mathcal{G}^{m}])`

In a graph classification task, a very common area for using graph kernels what we need to calculate in general are the matrices :math:`\mathcal{K}^{n\times n}` and :math:`\mathcal{K}^{m\times n}`.

.. currentmodule:: grakel

Parametrization
---------------

Any :code:`Kernel` inherited object comes with three parameters:

- :code:`verbose` is a :code:`bool` parameter in order for the kernel to print messages that are related with the progress of the execution.
- :code:`normalize` parameter is a :code:`bool` parameter which ensures that the kernel output will be normalized, that is :math:`[\mathcal{\hat{K}}]_{ij} = \frac{[\mathcal{K}]_{ij}}{\sqrt[]{[\mathcal{K}]_{ii}*[\mathcal{K}]_{jj}}}`.
- :code:`n_jobs` is an :code:`int` parameter that defines the amount of parallel jobs on which parts of the kernel calculation will be executed (if a parallelization has been implemented).

.. note::

    In order for normalization to happen even in a framework scheme, its :code:`Kernel` should have an implemented :code:`diagonal` method.

The :code:`Kernel` class as discussed above can be found below:

.. autosummary::

   Kernel
