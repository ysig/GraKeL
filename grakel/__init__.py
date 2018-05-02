"""Init file for the whole grakel project."""
from grakel import datasets

from grakel.graph import Graph

from grakel.graph_kernels import GraphKernel


from grakel.kernels import Kernel

from grakel.kernels import GraphletSampling
from grakel.kernels import RandomWalk
from grakel.kernels import RandomWalkLabeled
from grakel.kernels import ShortestPath
from grakel.kernels import ShortestPathAttr
from grakel.kernels import WeisfeilerLehman
from grakel.kernels import NeighborhoodHash
from grakel.kernels import PyramidMatch
from grakel.kernels import SubgraphMatching
from grakel.kernels import NeighborhoodSubgraphPairwiseDistance
from grakel.kernels import LovaszTheta
from grakel.kernels import SvmTheta
from grakel.kernels import OddSth
from grakel.kernels import Propagation
from grakel.kernels import PropagationAttr
from grakel.kernels import HadamardCode
from grakel.kernels import MultiscaleLaplacian
from grakel.kernels import MultiscaleLaplacianFast
from grakel.kernels import VertexHistogram
from grakel.kernels import EdgeHistogram
from grakel.kernels import GraphHopper
from grakel.kernels import CoreFramework


__all__ = [
    "datasets",
    "GraphKernel",
    "Graph",
    "Kernel",
    "GraphletSampling",
    "RandomWalk",
    "RandomWalkLabeled",
    "ShortestPath",
    "ShortestPathAttr",
    "WeisfeilerLehman",
    "NeighborhoodHash",
    "PyramidMatch",
    "SubgraphMatching",
    "NeighborhoodSubgraphPairwiseDistance",
    "LovaszTheta",
    "SvmTheta",
    "OddSth",
    "Propagation",
    "PropagationAttr",
    "HadamardCode",
    "MultiscaleLaplacian",
    "MultiscaleLaplacianFast",
    "VertexHistogram",
    "EdgeHistogram",
    "GraphHopper",
    "CoreFramework"
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
__version__ = '0.1a3'
