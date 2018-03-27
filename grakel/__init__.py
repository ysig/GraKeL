"""Init file for the whole grakel project."""
from grakel import datasets

from grakel.graph import Graph

from grakel.graph_kernels import GraphKernel

from grakel.kernels import kernel
from grakel.kernels import vertex_histogram
from grakel.kernels import edge_histogram
from grakel.kernels import graphlet_sampling
from grakel.kernels import random_walk
from grakel.kernels import random_walk_labeled
from grakel.kernels import shortest_path
from grakel.kernels import shortest_path_attr
from grakel.kernels import weisfeiler_lehman
from grakel.kernels import neighborhood_hash
from grakel.kernels import pyramid_match
from grakel.kernels import subgraph_matching
from grakel.kernels import neighborhood_subgraph_pairwise_distance
from grakel.kernels import lovasz_theta
from grakel.kernels import svm_theta
from grakel.kernels import odd_sth
from grakel.kernels import propagation
from grakel.kernels import hadamard_code
from grakel.kernels import multiscale_laplacian
from grakel.kernels import multiscale_laplacian_fast

__all__ = [
    "datasets",
    "GraphKernel",
    "Graph",
    "kernel",
    "graphlet_sampling",
    "random_walk",
    "random_walk_labeled",
    "shortest_path",
    "shortest_path_attr",
    "weisfeiler_lehman",
    "neighborhood_hash",
    "pyramid_match",
    "subgraph_matching",
    "neighborhood_subgraph_pairwise_distance",
    "lovasz_theta",
    "svm_theta",
    "odd_sth",
    "propagation",
    "hadamard_code",
    "multiscale_laplacian",
    "multiscale_laplacian_fast",
    "vertex_histogram",
    "edge_histogram"
]

# Generic release markers:
#   X.Y
#   X.Y.Z   # For bugfix releases
#
# Admissible pre-release markers:
#   X.YaN   # Alpha release
#   X.YbN   # Beta release
#   X.YrcN  # Release Candidate
#   X.Y     # Final release
#
# Dev branch marker is: 'X.Y.dev' or 'X.Y.devN' where N is an integer.
# 'X.Y.dev0' is the canonical version of 'X.Y.dev'
#
__version__ = '0.1a2'
