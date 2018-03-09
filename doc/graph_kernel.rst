.. _graph_kernel:

GraphKernel (class)
===================
Documentation for the graph-kernel decorator.

The `GraphKernel` *Decorator*
-----------------------------

A specification follows for GraphKernel where each of the supported kernel is presented
with whether it utilizes information from node-labels, edge-labels, node-attributes, edge-attributes
inside the kernel calculation.

+--------------------------------------------------+---------+---------+----------------+----------------+
|                                                  |      Labels       |           Attributes            |
|                    Kernels                       +---------+---------+----------------+----------------+
|                                                  |  node   | edge    | node           | edge           |
+==================================================+=========+=========+================+================+
| random_walk                                      | Ignores | Ignores | Ignores        | Ignores        |
+--------------------------------------------------+---------+---------+----------------+----------------+
| shortest_path, with_labels=False                 | Ignores | Ignores | Ignores        | Ignores        |
+--------------------------------------------------+---------+---------+----------------+----------------+
| shortest_path, with_labels=True                  | Needs   | Ignores | Nonce-Behavior | Ignores        |
+--------------------------------------------------+---------+---------+----------------+----------------+
| shortest_path, as_attributes=True                | Slow    | Ignores | Needs          | Ignores        |
+--------------------------------------------------+---------+---------+----------------+----------------+
| graphlet_sampling                                | Ignores | Ignores | Ignores        | Ignores        |
+--------------------------------------------------+---------+---------+----------------+----------------+
| multiscale_laplacian                             | Forbids | Ignores | Needs          | Ignores        |
+--------------------------------------------------+---------+---------+----------------+----------------+
| subgraph_matching                                | Forbids | Ignores | Needs          | Ignores        |
+--------------------------------------------------+---------+---------+----------------+----------------+
| lovasz_theta                                     | Ignores | Ignores | Ignores        | Ignores        |
+--------------------------------------------------+---------+---------+----------------+----------------+
| svm_theta                                        | Ignores | Ignores | Ignores        | Ignores        |
+--------------------------------------------------+---------+---------+----------------+----------------+
| neighborhood_hash                                | Needs   | Ignores | Nonce-Behavior | Ignores        |
+--------------------------------------------------+---------+---------+----------------+----------------+
| neighborhood_subgraph_pairwise_distance, NSPD    | Needs   | Needs   | Nonce-Behavior | Nonce-Behavior |
+--------------------------------------------------+---------+---------+----------------+----------------+
| odd_sth                                          | Needs   | Ignores | Nonce-Behavior | Ignores        |
+--------------------------------------------------+---------+---------+----------------+----------------+
| propagation                                      | Needs   | Ignores | Nonce-Behavior | Ignores        |
+--------------------------------------------------+---------+---------+----------------+----------------+
| pyramid_match                                    | Needs   | Ignores | Nonce-Behavior | Ignores        |
+--------------------------------------------------+---------+---------+----------------+----------------+
| vertex_histogram, subtree_wl                     | Needs   | Ignores | Nonce-Behavior | Ignores        |
+--------------------------------------------------+---------+---------+----------------+----------------+
| edge_histogram                                   | Ignores | Needs   | Ignores        | Nonce-Behavior |
+--------------------------------------------------+---------+---------+----------------+----------------+

.. currentmodule:: grakel

.. autosummary::
   :toctree: generated/
   :template: class.rst

   grakel.GraphKernel
