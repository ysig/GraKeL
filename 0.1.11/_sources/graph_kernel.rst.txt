.. _graph_kernel:

GraphKernel (class)
===================
Documentation for the graph-kernel generic wrapper.

The `GraphKernel` *Decorator*
-----------------------------

A specification follows for GraphKernel where each of the supported kernel is presented
with whether it utilizes information from node-labels, edge-labels, node-attributes, edge-attributes
inside the kernel calculation.

+----------------------------------------------------------+----------------+---------------+----------------+----------------+
|                                                          |      Labels                    |           Attributes            |
|                    Kernels                               +----------------+---------------+----------------+----------------+
|                                                          |  node          | edge          | node           | edge           |
+==========================================================+================+===============+================+================+
| random_walk, [:code:`with_labels=False`]                 | Ignores        | Ignores       | Ignores        | Ignores        |
+----------------------------------------------------------+----------------+---------------+----------------+----------------+
| random_walk, :code:`with_labels=True`                    | Needs          | Ignores       | Ignores        | Ignores        |
+----------------------------------------------------------+----------------+---------------+----------------+----------------+
| shortest_path, [:code:`with_labels=False`]               | Ignores        | Ignores       | Ignores        | Ignores        |
+----------------------------------------------------------+----------------+---------------+----------------+----------------+
| shortest_path, :code:`with_labels=True`                  | Needs          | Ignores       | Nonce-Behavior | Ignores        |
+----------------------------------------------------------+----------------+---------------+----------------+----------------+
| shortest_path, :code:`as_attributes=True`                | Slow           | Ignores       | Needs          | Ignores        |
+----------------------------------------------------------+----------------+---------------+----------------+----------------+
| graphlet_sampling                                        | Ignores        | Ignores       | Ignores        | Ignores        |
+----------------------------------------------------------+----------------+---------------+----------------+----------------+
| multiscale_laplacian                                     | Forbids        | Ignores       | Needs          | Ignores        |
+----------------------------------------------------------+----------------+---------------+----------------+----------------+
| subgraph_matching, [:code:`type(kv)=type(ke)='labels'`]  | type(kv)       | type(ke)      | type(kv)       | type(ke)       |
+----------------------------------------------------------+----------------+---------------+----------------+----------------+
| lovasz_theta                                             | Ignores        | Ignores       | Ignores        | Ignores        |
+----------------------------------------------------------+----------------+---------------+----------------+----------------+
| svm_theta                                                | Ignores        | Ignores       | Ignores        | Ignores        |
+----------------------------------------------------------+----------------+---------------+----------------+----------------+
| neighborhood_hash                                        | Needs          | Ignores       | Nonce-Behavior | Ignores        |
+----------------------------------------------------------+----------------+---------------+----------------+----------------+
| neighborhood_subgraph_pairwise_distance, NSPD            | Needs          | Needs         | Nonce-Behavior | Nonce-Behavior |
+----------------------------------------------------------+----------------+---------------+----------------+----------------+
| odd_sth                                                  | Needs          | Ignores       | Nonce-Behavior | Ignores        |
+----------------------------------------------------------+----------------+---------------+----------------+----------------+
| propagation, [:code:`with_attributes=False`]             | Needs          | Ignores       | Nonce-Behavior | Ignores        |
+----------------------------------------------------------+----------------+---------------+----------------+----------------+
| propagation, :code:`with_attributes=True`                | Nonce-Behavior | Ignores       | Needs          | Ignores        |
+----------------------------------------------------------+----------------+---------------+----------------+----------------+
| pyramid_match                                            | Needs          | Ignores       | Nonce-Behavior | Ignores        |
+----------------------------------------------------------+----------------+---------------+----------------+----------------+
| vertex_histogram, subtree_wl                             | Needs          | Ignores       | Nonce-Behavior | Ignores        |
+----------------------------------------------------------+----------------+---------------+----------------+----------------+
| edge_histogram                                           | Ignores        | Needs         | Ignores        | Nonce-Behavior |
+----------------------------------------------------------+----------------+---------------+----------------+----------------+
| graph_hopper                                             | Nonce-Behavior | Ignores       | Needs          | Ignores        |
+----------------------------------------------------------+----------------+---------------+----------------+----------------+
| weisfeiler_lehman                                        | Needs          | *base_kernel* | *base_kernel*  | *base_kernel*  |
+----------------------------------------------------------+----------------+---------------+----------------+----------------+
| hadamard_code                                            | Needs          | *base_kernel* | *base_kernel*  | *base_kernel*  |
+----------------------------------------------------------+----------------+---------------+----------------+----------------+
| core_framework                                           | *base_kernel*  | *base_kernel* | *base_kernel*  | *base_kernel*  |
+----------------------------------------------------------+----------------+---------------+----------------+----------------+

where **[**, **]** denounce the default value and the *kernel-names* and *arguments* are as given in the :code:`__init__` of the
:code:`GraphKernel` object and **base_kernel** signifies means that given input with this type of *meta-information* the behavior
of the kernel should **only** depend on the behavior of the base kernel to this type of *meta-information*.

.. currentmodule:: grakel

.. autosummary::
   :toctree: generated/
   :template: class.rst

   GraphKernel
