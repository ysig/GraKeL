"""General tests, concerning sci-kit compatibility."""
# Author: Ioannis Siglidis <y.siglidis@gmail.com>
# License: BSD 3 clause
# Currently checking only for picklability.

import warnings
import numpy as np
import pickle

from grakel.datasets import generate_dataset

from grakel import GraphKernel
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
from grakel.kernels import VertexHistogram
from grakel.kernels import EdgeHistogram
from grakel.kernels import GraphHopper
from grakel.kernels import CoreFramework
from grakel.kernels import WeisfeilerLehmanOptimalAssignment

verbose, normalize = False, True
default_eigvalue_precision = float("-1e-5")
rs = np.random.RandomState(42)
warnings.filterwarnings("ignore")

cvxopt = True
try:
    import cvxopt
except ImportError:
    cvxopt = False


def is_picklable(obj):
    try:
        pickle.dumps(obj)
    except Exception as ex:
        return False, ex
    return True


def test_random_walk():
    """Picklability test for the Simple Random Walk kernel."""
    train, _ = generate_dataset(n_graphs=100,
                                r_vertices=(10, 20),
                                r_connectivity=(0.4, 0.8),
                                r_weight_edges=(0.01, 12.0),
                                n_graphs_test=40,
                                random_state=rs,
                                features=None)

    rw_kernel = RandomWalk(verbose=verbose, normalize=normalize)
    rw_kernel.fit(train)
    assert is_picklable(rw_kernel)


def test_random_walk_labels():
    """Picklability test for the Labelled Random Walk kernel."""
    train, _ = generate_dataset(n_graphs=100,
                                r_vertices=(10, 20),
                                r_connectivity=(0.4, 0.8),
                                r_weight_edges=(0.01, 12.0),
                                n_graphs_test=40,
                                random_state=rs,
                                features=('nl', 3))

    rw_kernel = RandomWalkLabeled(verbose=verbose, normalize=normalize)
    rw_kernel.fit(train)
    assert is_picklable(rw_kernel)


def test_shortest_path():
    """Picklability test for the Shortest Path kernel."""
    train, _ = generate_dataset(n_graphs=100,
                                r_vertices=(10, 20),
                                r_connectivity=(0.4, 0.8),
                                r_weight_edges=(1, 1),
                                n_graphs_test=40,
                                random_state=rs,
                                features=('nl', 3))

    sp_kernel = ShortestPath(verbose=verbose, normalize=normalize)

    sp_kernel.fit(train)
    assert is_picklable(sp_kernel)

    train, _ = generate_dataset(n_graphs=50,
                                r_vertices=(5, 10),
                                r_connectivity=(0.4, 0.8),
                                r_weight_edges=(1, 1),
                                n_graphs_test=20,
                                random_state=rs,
                                features=('na', 5))

    sp_kernel = ShortestPathAttr(verbose=verbose, normalize=normalize)

    sp_kernel.fit(train)
    assert is_picklable(sp_kernel)


def test_graphlet_sampling():
    """Picklability test for the Graphlet Sampling Kernel [+ generic-wrapper]."""
    train, _ = generate_dataset(n_graphs=100,
                                r_vertices=(10, 20),
                                r_connectivity=(0.4, 0.8),
                                r_weight_edges=(1, 1),
                                n_graphs_test=40,
                                random_state=rs,
                                features=('nl', 3))

    gs_kernel = GraphletSampling(verbose=verbose, normalize=normalize, sampling=dict(n_samples=50))
    gk = GraphKernel(kernel={"name": "graphlet_sampling",
                             "sampling": {"n_samples": 50}},
                     verbose=verbose, normalize=normalize)
    gs_kernel.fit(train)
    assert is_picklable(gs_kernel)
    gk.fit(train)
    assert is_picklable(gk)


def test_weisfeiler_lehman():
    """Picklability test for the Weisfeiler Lehman kernel."""
    train, _ = generate_dataset(n_graphs=100,
                                r_vertices=(10, 20),
                                r_connectivity=(0.4, 0.8),
                                r_weight_edges=(1, 1),
                                n_graphs_test=40,
                                random_state=rs,
                                features=('nl', 3))

    wl_st_kernel = WeisfeilerLehman(verbose=verbose, normalize=normalize,
                                    base_graph_kernel=VertexHistogram)
    wl_st_kernel.fit(train)
    assert is_picklable(wl_st_kernel)


def test_weisfeiler_lehman_optimal_assignment():
    """Picklability test for the Weisfeiler Lehman kernel."""
    train, _ = generate_dataset(n_graphs=100,
                                r_vertices=(10, 20),
                                r_connectivity=(0.4, 0.8),
                                r_weight_edges=(1, 1),
                                n_graphs_test=40,
                                random_state=rs,
                                features=('nl', 3))

    wl_oa_kernel = WeisfeilerLehmanOptimalAssignment(verbose=verbose, normalize=normalize)
    wl_oa_kernel.fit(train)
    assert is_picklable(wl_oa_kernel)


def test_pyramid_match():
    """Picklability test for the Pyramid Match kernel [+ generic-wrapper]."""
    train, _ = generate_dataset(n_graphs=100,
                                r_vertices=(10, 20),
                                r_connectivity=(0.4, 0.8),
                                r_weight_edges=(1, 1),
                                n_graphs_test=40,
                                random_state=rs,
                                features=('nl', 3))

    pm_kernel = PyramidMatch(verbose=verbose, normalize=normalize)
    gk = GraphKernel(kernel={"name": "pyramid_match"}, verbose=verbose,
                     normalize=normalize)
    pm_kernel.fit(train)
    assert is_picklable(pm_kernel)
    gk.fit(train)
    assert is_picklable(gk)


def test_pyramid_match_no_labels():
    """Picklability test for the Pyramid Match kernel with no labels [+ generic-wrapper]."""
    train, _ = generate_dataset(n_graphs=100,
                                r_vertices=(10, 20),
                                r_connectivity=(0.4, 0.8),
                                r_weight_edges=(1, 1),
                                n_graphs_test=40,
                                random_state=rs,
                                features=None)

    pm_kernel = PyramidMatch(verbose=verbose, normalize=normalize, with_labels=False)
    gk = GraphKernel(kernel={"name": "pyramid_match", "with_labels": False},
                     verbose=verbose, normalize=normalize)
    pm_kernel.fit(train)
    assert is_picklable(pm_kernel)
    gk.fit(train)
    assert is_picklable(gk)


def test_neighborhood_hash():
    """Picklability test for the Neighborhood Hash kernel."""
    train, _ = generate_dataset(n_graphs=100,
                                r_vertices=(10, 20),
                                r_connectivity=(0.4, 0.8),
                                r_weight_edges=(1, 1),
                                n_graphs_test=40,
                                random_state=rs,
                                features=('nl', 3))

    nh_kernel = NeighborhoodHash(verbose=verbose, normalize=normalize)
    nh_kernel.fit(train)
    assert is_picklable(nh_kernel)


def test_subgraph_matching():
    """Picklability test for the Subgraph Matching kernel."""
    # node-label/edge-label
    train, _ = generate_dataset(n_graphs=100,
                                r_vertices=(10, 20),
                                r_connectivity=(0.4, 0.8),
                                r_weight_edges=(1, 1),
                                n_graphs_test=40,
                                random_state=rs,
                                features=('nl', 3, 'el', 4))

    sm_kernel = SubgraphMatching(verbose=verbose, normalize=normalize)
    sm_kernel.fit(train)
    assert is_picklable(sm_kernel)

    # node-label/edge-attribute
    train, _ = generate_dataset(n_graphs=50,
                                r_vertices=(5, 10),
                                r_connectivity=(0.4, 0.8),
                                r_weight_edges=(1, 1),
                                n_graphs_test=20,
                                random_state=rs,
                                features=('nl', 3, 'ea', 5))

    sm_kernel = SubgraphMatching(verbose=verbose, normalize=normalize, ke=np.dot)
    sm_kernel.fit(train)
    assert is_picklable(sm_kernel)

    # node-attribute/edge-label
    train, _ = generate_dataset(n_graphs=50,
                                r_vertices=(5, 10),
                                r_connectivity=(0.4, 0.8),
                                r_weight_edges=(1, 1),
                                n_graphs_test=20,
                                random_state=rs,
                                features=('na', 4, 'el', 3))

    sm_kernel = SubgraphMatching(verbose=verbose, normalize=normalize, kv=np.dot)
    sm_kernel.fit(train)
    assert is_picklable(sm_kernel)

    # node-attribute/edge-attribute
    train, _ = generate_dataset(n_graphs=50,
                                r_vertices=(5, 10),
                                r_connectivity=(0.4, 0.8),
                                r_weight_edges=(1, 1),
                                n_graphs_test=20,
                                random_state=rs,
                                features=('na', 4, 'ea', 6))

    sm_kernel = SubgraphMatching(verbose=verbose, normalize=normalize, ke=np.dot, kv=np.dot)
    sm_kernel.fit(train)
    assert is_picklable(sm_kernel)


def test_neighborhood_subgraph_pairwise_distance():
    """Picklability test for the Neighborhood Subgraph Pairwise Distance kernel [+ generic-wrapper]."""
    train, _ = generate_dataset(n_graphs=100,
                                r_vertices=(5, 10),
                                r_connectivity=(0.4, 0.8),
                                r_weight_edges=(1, 1),
                                n_graphs_test=40,
                                random_state=rs,
                                features=('nl', 5, 'el', 4))

    nspd_kernel = NeighborhoodSubgraphPairwiseDistance(verbose=verbose, normalize=normalize)
    gk = GraphKernel(kernel={
        "name": "neighborhood_subgraph_pairwise_distance"},
                     verbose=verbose, normalize=normalize)
    nspd_kernel.fit(train)
    assert is_picklable(nspd_kernel)

    gk.fit(train)
    assert is_picklable(gk)


if cvxopt:
    def test_lovasz_theta():
        """Picklability test for the Lovasz-theta distance kernel."""
        train, _ = generate_dataset(n_graphs=50,
                                    r_vertices=(5, 10),
                                    r_connectivity=(0.4, 0.8),
                                    r_weight_edges=(1, 1),
                                    n_graphs_test=20,
                                    random_state=rs,
                                    features=None)

        lt_kernel = LovaszTheta(verbose=verbose, normalize=normalize)
        lt_kernel.fit(train)
        assert is_picklable(lt_kernel)


def test_svm_theta():
    """Picklability test for the SVM-theta distance kernel."""
    train, _ = generate_dataset(n_graphs=100,
                                r_vertices=(10, 20),
                                r_connectivity=(0.4, 0.8),
                                r_weight_edges=(1, 1),
                                n_graphs_test=40,
                                random_state=rs,
                                features=None)

    svm_kernel = SvmTheta(verbose=verbose, normalize=normalize)
    svm_kernel.fit(train)
    assert is_picklable(svm_kernel)


def test_odd_sth():
    """Picklability test for the ODD-STh kernel [+ generic-wrapper]."""
    train, _ = generate_dataset(n_graphs=100,
                                r_vertices=(10, 20),
                                r_connectivity=(0.4, 0.8),
                                r_weight_edges=(1, 1),
                                n_graphs_test=40,
                                random_state=rs,
                                features=('nl', 4))

    odd_sth_kernel = OddSth(verbose=verbose, normalize=normalize)
    gk = GraphKernel(kernel={"name": "odd_sth"},
                     verbose=verbose, normalize=normalize)

    odd_sth_kernel.fit(train)
    assert is_picklable(odd_sth_kernel)
    gk.fit(train)
    assert is_picklable(gk)


def test_propagation():
    """Picklability test for the Propagation kernel."""
    train, _ = generate_dataset(n_graphs=100,
                                r_vertices=(10, 20),
                                r_connectivity=(0.4, 0.8),
                                r_weight_edges=(float("1e-5"), 10),
                                n_graphs_test=40,
                                random_state=rs,
                                features=('nl', 4))

    propagation_kernel = Propagation(verbose=verbose, normalize=normalize)
    propagation_kernel.fit(train)
    assert is_picklable(propagation_kernel)

    train, _ = generate_dataset(n_graphs=100,
                                r_vertices=(10, 20),
                                r_connectivity=(0.4, 0.8),
                                r_weight_edges=(float("1e-5"), 10),
                                n_graphs_test=40,
                                random_state=rs,
                                features=('na', 5))

    propagation_kernel_attr = PropagationAttr(verbose=verbose, normalize=normalize)
    propagation_kernel_attr.fit(train)
    assert is_picklable(propagation_kernel_attr)


def test_hadamard_code():
    """Picklability test for the Hadamard Code kernel."""
    train, _ = generate_dataset(n_graphs=100,
                                r_vertices=(10, 20),
                                r_connectivity=(0.4, 0.8),
                                r_weight_edges=(1, 1),
                                n_graphs_test=40,
                                random_state=rs,
                                features=('nl', 5))

    hadamard_code_kernel = HadamardCode(verbose=verbose, normalize=normalize,
                                        base_graph_kernel=VertexHistogram)
    hadamard_code_kernel.fit(train)
    assert is_picklable(hadamard_code_kernel)


def test_multiscale_laplacian():
    """Picklability test for the Fast Multiscale Laplacian kernel."""
    # Initialise kernel
    train, _ = generate_dataset(n_graphs=100,
                                r_vertices=(10, 20),
                                r_connectivity=(0.4, 0.8),
                                r_weight_edges=(1, 1),
                                n_graphs_test=40,
                                random_state=rs,
                                features=('na', 5))

    mlf_kernel = MultiscaleLaplacian(verbose=verbose, normalize=normalize)

    mlf_kernel.fit(train)
    assert is_picklable(mlf_kernel)


def test_vertex_histogram():
    """Picklability test for the Vertex Histogram Kernel [+ generic-wrapper]."""
    train, _ = generate_dataset(n_graphs=100,
                                r_vertices=(10, 20),
                                r_connectivity=(0.4, 0.8),
                                r_weight_edges=(1, 1),
                                n_graphs_test=40,
                                random_state=rs,
                                features=('nl', 5))

    vh_kernel = VertexHistogram(verbose=verbose, normalize=normalize)
    gk = GraphKernel(kernel={"name": "vertex_histogram"},
                     verbose=verbose, normalize=normalize)

    vh_kernel.fit(train)
    assert is_picklable(vh_kernel)
    gk.fit(train)
    assert is_picklable(gk)


def test_edge_histogram():
    """Picklability test for the Edge Histogram kernel [+ generic-wrapper]."""
    train, _ = generate_dataset(n_graphs=100,
                                r_vertices=(10, 20),
                                r_connectivity=(0.4, 0.8),
                                r_weight_edges=(1, 1),
                                n_graphs_test=40,
                                random_state=rs,
                                features=('el', 4))

    eh_kernel = EdgeHistogram(verbose=verbose, normalize=normalize)
    gk = GraphKernel(kernel={"name": "edge_histogram"},
                     verbose=verbose, normalize=normalize)

    eh_kernel.fit(train)
    assert is_picklable(eh_kernel)
    gk.fit(train)
    assert is_picklable(gk)


def test_graph_hopper():
    """Picklability test for the Graph Hopper kernel."""
    train, _ = generate_dataset(n_graphs=100,
                                r_vertices=(10, 20),
                                r_connectivity=(0.4, 0.8),
                                r_weight_edges=(1, 1),
                                n_graphs_test=40,
                                random_state=rs,
                                features=('na', 4))

    gh_kernel = GraphHopper(verbose=verbose, normalize=normalize)
    gh_kernel.fit(train)
    assert is_picklable(gh_kernel)


def test_core_framework():
    """Picklability test for the Core kernel Framework [+ generic-wrapper]."""
    train, _ = generate_dataset(n_graphs=100,
                                r_vertices=(10, 20),
                                r_connectivity=(0.4, 0.8),
                                r_weight_edges=(1, 1),
                                n_graphs_test=40,
                                random_state=rs,
                                features=('nl', 4))

    base_graph_kernel = (WeisfeilerLehman, dict(base_graph_kernel=VertexHistogram))
    core_framework = CoreFramework(verbose=verbose, normalize=normalize, base_graph_kernel=base_graph_kernel)

    kernel = [{"name": "core_framework"}, {"name": "weisfeiler_lehman"}, {"name": "vertex_histogram"}]
    gk = GraphKernel(kernel=kernel, verbose=verbose, normalize=normalize)
    core_framework.fit(train)
    assert is_picklable(core_framework)
    gk.fit(train)
    assert is_picklable(gk)


if __name__ == "__main__":
    warnings.filterwarnings("once")
    verbose = True

    test_random_walk()
    test_random_walk_labels()
    test_shortest_path()
    test_graphlet_sampling()
    test_weisfeiler_lehman()
    test_weisfeiler_lehman_optimal_assignment()
    test_pyramid_match()
    test_pyramid_match_no_labels()
    test_neighborhood_hash()
    test_subgraph_matching()
    test_neighborhood_subgraph_pairwise_distance()
    if cvxopt:
        test_lovasz_theta()
    test_svm_theta()
    test_odd_sth()
    test_propagation()
    test_hadamard_code()
    test_multiscale_laplacian()
    test_vertex_histogram()
    test_edge_histogram()
    test_graph_hopper()
    test_core_framework()
