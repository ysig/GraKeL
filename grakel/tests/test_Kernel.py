"""A random input test for Kernel Oblects similar to test_estimator."""
# Author: Ioannis Siglidis <y.siglidis@gmail.com>
# License: BSD 3 clause

import warnings
import numpy as np

from grakel.datasets import generate_dataset

from grakel import GraphKernel
from grakel.kernels import GraphletSampling
from grakel.kernels import RandomWalk
from grakel.kernels import RandomWalkLabeled
from grakel.kernels import ShortestPath
from grakel.kernels import ShortestPathAttr
from grakel.kernels import WeisfeilerLehman
from grakel.kernels import WeisfeilerLehmanOptimalAssignment
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

verbose, normalize = False, True
default_eigvalue_precision = float("-1e-5")
rs = np.random.RandomState(42)
warnings.filterwarnings("ignore")

cvxopt = True
try:
    import cvxopt
except ImportError:
    cvxopt = False


def test_random_walk():
    """Random input test for the Simple Random Walk kernel."""
    train, test = generate_dataset(n_graphs=100,
                                   r_vertices=(10, 20),
                                   r_connectivity=(0.4, 0.8),
                                   r_weight_edges=(0.01, 12.0),
                                   n_graphs_test=40,
                                   random_state=rs,
                                   features=None)

    rw_kernel = RandomWalk(verbose=verbose, normalize=normalize)
    try:
        rw_kernel.fit_transform(train)
        rw_kernel.transform(test)
        assert True
    except Exception as exception:
        assert False, exception


def test_random_walk_pd():
    """Random input test for the Simple Random Walk kernel [n_jobs=-1/generic-wrapper]."""
    train, test = generate_dataset(n_graphs=100,
                                   r_vertices=(10, 20),
                                   r_connectivity=(0.4, 0.8),
                                   r_weight_edges=(0.01, 12.0),
                                   n_graphs_test=40,
                                   random_state=rs,
                                   features=None)

    gk = GraphKernel(kernel="RW", verbose=verbose, normalize=normalize, n_jobs=-1)

    try:
        gk.fit_transform(train)
        gk.transform(test)
        assert True
    except Exception as exception:
        assert False, exception


def test_random_walk_labels():
    """Random input test for the Labelled Random Walk kernel."""
    train, test = generate_dataset(n_graphs=100,
                                   r_vertices=(10, 20),
                                   r_connectivity=(0.4, 0.8),
                                   r_weight_edges=(0.01, 12.0),
                                   n_graphs_test=40,
                                   random_state=rs,
                                   features=('nl', 3))

    rw_kernel = RandomWalkLabeled(verbose=verbose, normalize=normalize)

    try:
        rw_kernel.fit_transform(train)
        rw_kernel.transform(test)
        assert True
    except Exception as exception:
        assert False, exception


def test_random_walk_labels_pd():
    """Random input test for the Labelled Random Walk kernel [n_jobs=-1/generic-wrapper]."""
    train, test = generate_dataset(n_graphs=100,
                                   r_vertices=(10, 20),
                                   r_connectivity=(0.4, 0.8),
                                   r_weight_edges=(0.01, 12.0),
                                   n_graphs_test=40,
                                   random_state=rs,
                                   features=('nl', 3))

    gk = GraphKernel(kernel={"name": "RW", "with_labels": True},
                     verbose=verbose, normalize=normalize, )

    try:
        gk.fit_transform(train)
        gk.transform(test)
        assert True
    except Exception as exception:
        assert False, exception


def test_shortest_path():
    """Random input test for the Shortest Path kernel."""
    train, test = generate_dataset(n_graphs=100,
                                   r_vertices=(10, 20),
                                   r_connectivity=(0.4, 0.8),
                                   r_weight_edges=(1, 1),
                                   n_graphs_test=40,
                                   random_state=rs,
                                   features=('nl', 3))

    sp_kernel = ShortestPath(verbose=verbose, normalize=normalize)

    try:
        sp_kernel.fit_transform(train)
        sp_kernel.transform(test)
        assert True
    except Exception as exception:
        assert False, exception

    train, test = generate_dataset(n_graphs=50,
                                   r_vertices=(5, 10),
                                   r_connectivity=(0.4, 0.8),
                                   r_weight_edges=(1, 1),
                                   n_graphs_test=20,
                                   random_state=rs,
                                   features=('na', 5))

    sp_kernel = ShortestPathAttr(verbose=verbose, normalize=normalize)

    try:
        sp_kernel.fit_transform(train)
        sp_kernel.transform(test)
        assert True
    except Exception as exception:
        assert False, exception


def test_shortest_path_pd():
    """Random input test for the Shortest Path kernel [n_jobs=-1 (for attributed)/decorator]."""
    train, test = generate_dataset(n_graphs=100,
                                   r_vertices=(10, 20),
                                   r_connectivity=(0.4, 0.8),
                                   r_weight_edges=(1, 1),
                                   n_graphs_test=40,
                                   random_state=rs,
                                   features=('nl', 3))

    gk = GraphKernel(kernel="SP", verbose=verbose, normalize=normalize)

    try:
        gk.fit_transform(train)
        gk.transform(test)
        assert True
    except Exception as exception:
        assert False, exception

    train, test = generate_dataset(n_graphs=50,
                                   r_vertices=(5, 10),
                                   r_connectivity=(0.4, 0.8),
                                   r_weight_edges=(1, 1),
                                   n_graphs_test=20,
                                   random_state=rs,
                                   features=('na', 5))

    gk = GraphKernel(kernel={"name": "SP", "as_attributes": True},
                     verbose=verbose, normalize=normalize, n_jobs=-1)

    try:
        gk.fit_transform(train)
        gk.transform(test)
        assert True
    except Exception as exception:
        assert False, exception


def test_graphlet_sampling():
    """Random input test for the Graphlet Sampling Kernel [+ generic-wrapper]."""
    train, test = generate_dataset(n_graphs=100,
                                   r_vertices=(10, 20),
                                   r_connectivity=(0.4, 0.8),
                                   r_weight_edges=(1, 1),
                                   n_graphs_test=40,
                                   random_state=rs,
                                   features=('nl', 3))

    gs_kernel = GraphletSampling(verbose=verbose, normalize=normalize, sampling=dict(n_samples=50))
    gk = GraphKernel(kernel={"name": "GR", "sampling": {"n_samples": 50}},
                     verbose=verbose, normalize=normalize)

    try:
        gs_kernel.fit_transform(train)
        gs_kernel.transform(test)
        gk.fit_transform(train)
        gk.transform(test)
        assert True
    except Exception as exception:
        assert False, exception


def test_weisfeiler_lehman():
    """Random input test for the Weisfeiler Lehman kernel."""
    train, test = generate_dataset(n_graphs=100,
                                   r_vertices=(10, 20),
                                   r_connectivity=(0.4, 0.8),
                                   r_weight_edges=(1, 1),
                                   n_graphs_test=40,
                                   random_state=rs,
                                   features=('nl', 3))

    wl_st_kernel = WeisfeilerLehman(verbose=verbose, normalize=normalize,
                                    base_graph_kernel=VertexHistogram)

    try:
        wl_st_kernel.fit_transform(train)
        wl_st_kernel.transform(test)
        assert True
    except Exception as exception:
        assert False, exception


def test_weisfeiler_lehman_optimal_assignment():
    """Random input test for the Weisfeiler Lehman Optimal Assignment kernel."""
    train, test = generate_dataset(n_graphs=100,
                                   r_vertices=(10, 20),
                                   r_connectivity=(0.4, 0.8),
                                   r_weight_edges=(1, 1),
                                   n_graphs_test=40,
                                   random_state=rs,
                                   features=('nl', 3))

    wl_oa_kernel = WeisfeilerLehmanOptimalAssignment(
        verbose=verbose, normalize=normalize)

    try:
        wl_oa_kernel.fit_transform(train)
        wl_oa_kernel.transform(test)
        assert True
    except Exception as exception:
        assert False, exception


def test_weisfeiler_lehman_pd():
    """Random input test for the Weisfeiler Lehman kernel [n_jobs=-1/generic-wrapper]."""
    train, test = generate_dataset(n_graphs=100,
                                   r_vertices=(10, 20),
                                   r_connectivity=(0.4, 0.8),
                                   r_weight_edges=(1, 1),
                                   n_graphs_test=40,
                                   random_state=rs,
                                   features=('nl', 3))

    gk = GraphKernel(kernel="WL", verbose=verbose, normalize=normalize)

    try:
        gk.fit_transform(train)
        gk.transform(test)
        assert True
    except Exception as exception:
        assert False, exception


def test_pyramid_match():
    """Random input test for the Pyramid Match kernel [+ generic-wrapper]."""
    train, test = generate_dataset(n_graphs=100,
                                   r_vertices=(10, 20),
                                   r_connectivity=(0.4, 0.8),
                                   r_weight_edges=(1, 1),
                                   n_graphs_test=40,
                                   random_state=rs,
                                   features=('nl', 3))

    pm_kernel = PyramidMatch(verbose=verbose, normalize=normalize)
    gk = GraphKernel("PM", verbose=verbose, normalize=normalize)

    try:
        pm_kernel.fit_transform(train)
        pm_kernel.transform(test)
        gk.fit_transform(train)
        gk.transform(test)
        assert True
    except Exception as exception:
        assert False, exception


def test_pyramid_match_no_labels():
    """Random input test for the Pyramid Match kernel with no labels [+ generic-wrapper]."""
    train, test = generate_dataset(n_graphs=100,
                                   r_vertices=(10, 20),
                                   r_connectivity=(0.4, 0.8),
                                   r_weight_edges=(1, 1),
                                   n_graphs_test=40,
                                   random_state=rs,
                                   features=None)

    pm_kernel = PyramidMatch(verbose=verbose, normalize=normalize, with_labels=False)
    gk = GraphKernel(kernel={"name": "PM", "with_labels": False},
                     verbose=verbose, normalize=normalize)

    try:
        pm_kernel.fit_transform(train)
        pm_kernel.transform(test)
        gk.fit_transform(train)
        gk.transform(test)
        assert True
    except Exception as exception:
        assert False, exception


def test_neighborhood_hash():
    """Random input test for the Neighborhood Hash kernel."""
    train, test = generate_dataset(n_graphs=100,
                                   r_vertices=(10, 20),
                                   r_connectivity=(0.4, 0.8),
                                   r_weight_edges=(1, 1),
                                   n_graphs_test=40,
                                   random_state=rs,
                                   features=('nl', 3))

    nh_kernel = NeighborhoodHash(verbose=verbose, normalize=normalize)

    try:
        nh_kernel.fit_transform(train)
        nh_kernel.transform(test)
        assert True
    except Exception as exception:
        assert False, exception


def test_neighborhood_hash_pd():
    """Random input test for the Neighborhood Hash kernel [n_jobs=-1/generic-wrapper]."""
    train, test = generate_dataset(n_graphs=100,
                                   r_vertices=(10, 20),
                                   r_connectivity=(0.4, 0.8),
                                   r_weight_edges=(1, 1),
                                   n_graphs_test=40,
                                   random_state=rs,
                                   features=('nl', 3))

    gk = GraphKernel(kernel="NH", verbose=verbose, normalize=normalize, n_jobs=-1)

    try:
        gk.fit_transform(train)
        gk.transform(test)
        assert True
    except Exception as exception:
        assert False, exception


def test_subgraph_matching():
    """Random input test for the Subgraph Matching kernel."""
    # node-label/edge-label
    train, test = generate_dataset(n_graphs=100,
                                   r_vertices=(10, 20),
                                   r_connectivity=(0.4, 0.8),
                                   r_weight_edges=(1, 1),
                                   n_graphs_test=40,
                                   random_state=rs,
                                   features=('nl', 3, 'el', 4))

    sm_kernel = SubgraphMatching(verbose=verbose, normalize=normalize)

    try:
        sm_kernel.fit_transform(train)
        sm_kernel.transform(test)
        assert True
    except Exception as exception:
        assert False, exception

    # node-label/edge-attribute
    train, test = generate_dataset(n_graphs=50,
                                   r_vertices=(5, 10),
                                   r_connectivity=(0.4, 0.8),
                                   r_weight_edges=(1, 1),
                                   n_graphs_test=20,
                                   random_state=rs,
                                   features=('nl', 3, 'ea', 5))

    sm_kernel = SubgraphMatching(verbose=verbose, normalize=normalize, ke=np.dot)

    try:
        sm_kernel.fit_transform(train)
        sm_kernel.transform(test)
        assert True
    except Exception as exception:
        assert False, exception

    # node-attribute/edge-label
    train, test = generate_dataset(n_graphs=50,
                                   r_vertices=(5, 10),
                                   r_connectivity=(0.4, 0.8),
                                   r_weight_edges=(1, 1),
                                   n_graphs_test=20,
                                   random_state=rs,
                                   features=('na', 4, 'el', 3))

    sm_kernel = SubgraphMatching(verbose=verbose, normalize=normalize, kv=np.dot)

    try:
        sm_kernel.fit_transform(train)
        sm_kernel.transform(test)
        assert True
    except Exception as exception:
        assert False, exception

    # node-attribute/edge-attribute
    train, test = generate_dataset(n_graphs=50,
                                   r_vertices=(5, 10),
                                   r_connectivity=(0.4, 0.8),
                                   r_weight_edges=(1, 1),
                                   n_graphs_test=20,
                                   random_state=rs,
                                   features=('na', 4, 'ea', 6))

    sm_kernel = SubgraphMatching(verbose=verbose, normalize=normalize, ke=np.dot, kv=np.dot)

    try:
        sm_kernel.fit_transform(train)
        sm_kernel.transform(test)
        assert True
    except Exception as exception:
        assert False, exception


def test_subgraph_matching_pd():
    """Random input test for the Subgraph Matching kernel [n_jobs=-1/generic-wrapper]."""
    # node-label/edge-label
    train, test = generate_dataset(n_graphs=100,
                                   r_vertices=(10, 20),
                                   r_connectivity=(0.4, 0.8),
                                   r_weight_edges=(1, 1),
                                   n_graphs_test=40,
                                   random_state=rs,
                                   features=('nl', 3, 'el', 4))

    gk = GraphKernel(kernel={"name": "SM"}, verbose=verbose,
                     normalize=normalize, n_jobs=-1)

    try:
        gk.fit_transform(train)
        gk.transform(test)
        assert True
    except Exception as exception:
        assert False, exception

    # node-label/edge-attribute
    train, test = generate_dataset(n_graphs=50,
                                   r_vertices=(5, 10),
                                   r_connectivity=(0.4, 0.8),
                                   r_weight_edges=(1, 1),
                                   n_graphs_test=20,
                                   random_state=rs,
                                   features=('nl', 3, 'ea', 5))

    gk = GraphKernel(kernel={"name": "SM", "ke": np.dot},
                     verbose=verbose, normalize=normalize, n_jobs=-1)

    try:
        gk.fit_transform(train)
        gk.transform(test)
        assert True
    except Exception as exception:
        assert False, exception

    # node-attribute/edge-label
    train, test = generate_dataset(n_graphs=50,
                                   r_vertices=(5, 10),
                                   r_connectivity=(0.4, 0.8),
                                   r_weight_edges=(1, 1),
                                   n_graphs_test=20,
                                   random_state=rs,
                                   features=('na', 4, 'el', 3))

    gk = GraphKernel(kernel={"name": "SM", "kv": np.dot},
                     verbose=verbose, normalize=normalize, n_jobs=-1)

    try:
        gk.fit_transform(train)
        gk.transform(test)
        assert True
    except Exception as exception:
        assert False, exception

    # node-attribute/edge-attribute
    train, test = generate_dataset(n_graphs=50,
                                   r_vertices=(5, 10),
                                   r_connectivity=(0.4, 0.8),
                                   r_weight_edges=(1, 1),
                                   n_graphs_test=20,
                                   random_state=rs,
                                   features=('na', 4, 'ea', 6))

    gk = GraphKernel(kernel={"name": "SM", "kv": np.dot, "ke": np.dot},
                     verbose=verbose, normalize=normalize, n_jobs=-1)

    try:
        gk.fit_transform(train)
        gk.transform(test)
        assert True
    except Exception as exception:
        assert False, exception


def test_neighborhood_subgraph_pairwise_distance():
    """Random input test for the Neighborhood Subgraph Pairwise Distance kernel [+ generic-wrapper]."""
    train, test = generate_dataset(n_graphs=100,
                                   r_vertices=(5, 10),
                                   r_connectivity=(0.4, 0.8),
                                   r_weight_edges=(1, 1),
                                   n_graphs_test=40,
                                   random_state=rs,
                                   features=('nl', 5, 'el', 4))

    nspd_kernel = NeighborhoodSubgraphPairwiseDistance(verbose=verbose, normalize=normalize)
    gk = GraphKernel(kernel="NSPD", verbose=verbose, normalize=normalize)

    try:
        nspd_kernel.fit_transform(train)
        nspd_kernel.transform(test)
        gk.fit_transform(train)
        gk.transform(test)
        assert True
    except Exception as exception:
        assert False, exception


if cvxopt:
    def test_lovasz_theta():
        """Random input test for the Lovasz-theta distance kernel."""
        train, test = generate_dataset(n_graphs=50,
                                       r_vertices=(5, 10),
                                       r_connectivity=(0.4, 0.8),
                                       r_weight_edges=(1, 1),
                                       n_graphs_test=20,
                                       random_state=rs,
                                       features=None)

        lt_kernel = LovaszTheta(verbose=verbose, normalize=normalize)

        try:
            lt_kernel.fit_transform(train)
            lt_kernel.transform(test)
            assert True
        except Exception as exception:
            assert False, exception

    def test_lovasz_theta_pd():
        """Random input test for the Lovasz-theta distance kernel [n_jobs=-1/generic-wrapper]."""
        train, test = generate_dataset(n_graphs=50,
                                       r_vertices=(5, 10),
                                       r_connectivity=(0.4, 0.8),
                                       r_weight_edges=(1, 1),
                                       n_graphs_test=20,
                                       random_state=rs,
                                       features=None)

        gk = GraphKernel(kernel="lovasz_theta",
                         verbose=verbose, normalize=normalize, n_jobs=-1)

        try:
            gk.fit_transform(train)
            gk.transform(test)
            assert True
        except Exception as exception:
            assert False, exception


def test_svm_theta():
    """Random input test for the SVM-theta distance kernel."""
    train, test = generate_dataset(n_graphs=100,
                                   r_vertices=(10, 20),
                                   r_connectivity=(0.4, 0.8),
                                   r_weight_edges=(1, 1),
                                   n_graphs_test=40,
                                   random_state=rs,
                                   features=None)

    svm_kernel = SvmTheta(verbose=verbose, normalize=normalize)

    try:
        svm_kernel.fit_transform(train)
        svm_kernel.transform(test)
        assert True
    except Exception as exception:
        assert False, exception


def test_svm_theta_pd():
    """Random input test for the SVM-theta distance kernel [n_jobs=-1/generic-wrapper]."""
    train, test = generate_dataset(n_graphs=100,
                                   r_vertices=(10, 20),
                                   r_connectivity=(0.4, 0.8),
                                   r_weight_edges=(1, 1),
                                   n_graphs_test=40,
                                   random_state=rs,
                                   features=None)

    gk = GraphKernel(kernel="svm_theta",
                     verbose=verbose, normalize=normalize)

    try:
        gk.fit_transform(train)
        gk.transform(test)
        assert True
    except Exception as exception:
        assert False, exception


def test_odd_sth():
    """Random input test for the ODD-STh kernel [+ generic-wrapper]."""
    train, test = generate_dataset(n_graphs=100,
                                   r_vertices=(10, 20),
                                   r_connectivity=(0.4, 0.8),
                                   r_weight_edges=(1, 1),
                                   n_graphs_test=40,
                                   random_state=rs,
                                   features=('nl', 4))

    odd_sth_kernel = OddSth(verbose=verbose, normalize=normalize)
    gk = GraphKernel(kernel="ODD",
                     verbose=verbose, normalize=normalize)

    try:
        odd_sth_kernel.fit_transform(train)
        odd_sth_kernel.transform(test)
        gk.fit_transform(train)
        gk.transform(test)
        assert True
    except Exception as exception:
        assert False, exception


def test_propagation():
    """Random input test for the Propagation kernel."""
    train, test = generate_dataset(n_graphs=100,
                                   r_vertices=(10, 20),
                                   r_connectivity=(0.4, 0.8),
                                   r_weight_edges=(float("1e-5"), 10),
                                   n_graphs_test=40,
                                   random_state=rs,
                                   features=('nl', 4))

    propagation_kernel = Propagation(verbose=verbose, normalize=normalize)

    try:
        propagation_kernel.fit_transform(train)
        propagation_kernel.transform(test)
        assert True
    except Exception as exception:
        assert False, exception

    train, test = generate_dataset(n_graphs=100,
                                   r_vertices=(10, 20),
                                   r_connectivity=(0.4, 0.8),
                                   r_weight_edges=(float("1e-5"), 10),
                                   n_graphs_test=40,
                                   random_state=rs,
                                   features=('na', 5))

    propagation_kernel_attr = PropagationAttr(verbose=verbose, normalize=normalize)

    try:
        propagation_kernel_attr.fit_transform(train)
        propagation_kernel_attr.transform(test)
        assert True
    except Exception as exception:
        assert False, exception


def test_propagation_pd():
    """Random input test for the Propagation kernel [n_jobs=-1/generic-wrapper]."""
    train, test = generate_dataset(n_graphs=100,
                                   r_vertices=(10, 20),
                                   r_connectivity=(0.4, 0.8),
                                   r_weight_edges=(float("1e-5"), 10),
                                   n_graphs_test=40,
                                   random_state=rs,
                                   features=('nl', 4))

    gk = GraphKernel(kernel="PR", verbose=verbose, normalize=normalize, n_jobs=-1)

    try:
        gk.fit_transform(train)
        gk.transform(test)
        assert True
    except Exception as exception:
        assert False, exception

    train, test = generate_dataset(n_graphs=100,
                                   r_vertices=(10, 20),
                                   r_connectivity=(0.4, 0.8),
                                   r_weight_edges=(float("1e-5"), 10),
                                   n_graphs_test=40,
                                   random_state=rs,
                                   features=('na', 5))

    gk = GraphKernel(kernel={"name": "PR", "with_attributes": True},
                     verbose=verbose, normalize=normalize, n_jobs=-1)

    try:
        gk.fit_transform(train)
        gk.transform(test)
        assert True
    except Exception as exception:
        assert False, exception


def test_hadamard_code():
    """Random input test for the Hadamard Code kernel."""
    train, test = generate_dataset(n_graphs=100,
                                   r_vertices=(10, 20),
                                   r_connectivity=(0.4, 0.8),
                                   r_weight_edges=(1, 1),
                                   n_graphs_test=40,
                                   random_state=rs,
                                   features=('nl', 5))

    hadamard_code_kernel = HadamardCode(verbose=verbose, normalize=normalize,
                                        base_graph_kernel=VertexHistogram)

    try:
        hadamard_code_kernel.fit_transform(train)
        hadamard_code_kernel.transform(test)
        assert True
    except Exception as exception:
        assert False, exception


def test_hadamard_code_pd():
    """Random input test for the Hadamard Code kernel [n_jobs=-1/generic-wrapper]."""
    train, test = generate_dataset(n_graphs=100,
                                   r_vertices=(10, 20),
                                   r_connectivity=(0.4, 0.8),
                                   r_weight_edges=(1, 1),
                                   n_graphs_test=40,
                                   random_state=rs,
                                   features=('nl', 5))

    gk = GraphKernel(kernel="HC",
                     verbose=verbose, normalize=normalize, n_jobs=-1)

    try:
        gk.fit_transform(train)
        gk.transform(test)
        assert True
    except Exception as exception:
        assert False, exception


def test_multiscale_laplacian():
    """Random input test for the Multiscale Laplacian kernel."""
    # Initialise kernel
    train, test = generate_dataset(n_graphs=100,
                                   r_vertices=(10, 20),
                                   r_connectivity=(0.4, 0.8),
                                   r_weight_edges=(1, 1),
                                   n_graphs_test=40,
                                   random_state=rs,
                                   features=('na', 5))

    mlf_kernel = MultiscaleLaplacian(verbose=verbose, normalize=normalize)

    try:
        mlf_kernel.fit_transform(train)
        mlf_kernel.transform(test)
        assert True
    except Exception as exception:
        assert False, exception


def test_multiscale_laplacian_pd():
    """Random input test for the Multiscale Laplacian kernel [n_jobs=-1/generic-wrapper]."""
    # Initialise kernel
    train, test = generate_dataset(n_graphs=100,
                                   r_vertices=(10, 20),
                                   r_connectivity=(0.4, 0.8),
                                   r_weight_edges=(1, 1),
                                   n_graphs_test=40,
                                   random_state=rs,
                                   features=('na', 5))

    gk = GraphKernel(kernel="ML", verbose=verbose, normalize=normalize, n_jobs=-1)

    try:
        gk.fit_transform(train)
        gk.transform(test)
        assert True
    except Exception as exception:
        assert False, exception


def test_vertex_histogram():
    """Random input test for the Vertex Histogram Kernel [+ generic-wrapper]."""
    train, test = generate_dataset(n_graphs=100,
                                   r_vertices=(10, 20),
                                   r_connectivity=(0.4, 0.8),
                                   r_weight_edges=(1, 1),
                                   n_graphs_test=40,
                                   random_state=rs,
                                   features=('nl', 5))

    vh_kernel = VertexHistogram(verbose=verbose, normalize=normalize)
    gk = GraphKernel(kernel="VH", verbose=verbose, normalize=normalize)

    try:
        vh_kernel.fit_transform(train)
        vh_kernel.transform(test)
        gk.fit_transform(train)
        gk.transform(test)
        assert True
    except Exception as exception:
        assert False, exception


def test_edge_histogram():
    """Random input test for the Edge Histogram kernel [+ generic-wrapper]."""
    train, test = generate_dataset(n_graphs=100,
                                   r_vertices=(10, 20),
                                   r_connectivity=(0.4, 0.8),
                                   r_weight_edges=(1, 1),
                                   n_graphs_test=40,
                                   random_state=rs,
                                   features=('el', 4))

    eh_kernel = EdgeHistogram(verbose=verbose, normalize=normalize)
    gk = GraphKernel(kernel="EH", verbose=verbose, normalize=normalize)

    try:
        eh_kernel.fit_transform(train)
        eh_kernel.transform(test)
        gk.fit_transform(train)
        gk.transform(test)
        assert True
    except Exception as exception:
        assert False, exception


def test_graph_hopper():
    """Random input test for the Graph Hopper kernel."""
    train, test = generate_dataset(n_graphs=100,
                                   r_vertices=(10, 20),
                                   r_connectivity=(0.4, 0.8),
                                   r_weight_edges=(1, 1),
                                   n_graphs_test=40,
                                   random_state=rs,
                                   features=('na', 4))

    gh_kernel = GraphHopper(verbose=verbose, normalize=normalize)

    try:
        gh_kernel.fit_transform(train)
        gh_kernel.transform(test)
        assert True
    except Exception as exception:
        assert False, exception


def test_graph_hopper_pd():
    """Random input test for the Graph Hopper kernel [n_jobs=-1/generic-wrapper]."""
    train, test = generate_dataset(n_graphs=100,
                                   r_vertices=(10, 20),
                                   r_connectivity=(0.4, 0.8),
                                   r_weight_edges=(1, 1),
                                   n_graphs_test=40,
                                   random_state=rs,
                                   features=('na', 4))

    gk = GraphKernel(kernel="GH",
                     verbose=verbose, normalize=normalize, n_jobs=-1)

    try:
        gk.fit_transform(train)
        gk.transform(test)
        assert True
    except Exception as exception:
        assert False, exception


def test_core_framework():
    """Random input test for the Core kernel Framework [+ generic-wrapper]."""
    train, test = generate_dataset(n_graphs=100,
                                   r_vertices=(10, 20),
                                   r_connectivity=(0.4, 0.8),
                                   r_weight_edges=(1, 1),
                                   n_graphs_test=40,
                                   random_state=rs,
                                   features=('nl', 4))

    base_graph_kernel = (WeisfeilerLehman, dict(base_graph_kernel=VertexHistogram))
    core_framework = CoreFramework(verbose=verbose, normalize=normalize, base_graph_kernel=base_graph_kernel)

    kernel = ["CORE", "WL"]
    gk = GraphKernel(kernel=kernel, verbose=verbose, normalize=normalize)
    try:
        core_framework.fit_transform(train)
        core_framework.transform(test)
        gk.fit_transform(train)
        gk.transform(test)
        assert True
    except Exception as exception:
        assert False, exception


if __name__ == "__main__":
    warnings.filterwarnings("once")
    verbose = True

    test_random_walk()
    test_random_walk_pd()
    test_random_walk_labels()
    test_random_walk_labels_pd()
    test_shortest_path()
    test_shortest_path_pd()
    test_graphlet_sampling()
    test_weisfeiler_lehman()
    test_weisfeiler_lehman_optimal_assignment()
    test_weisfeiler_lehman_pd()
    test_pyramid_match()
    test_pyramid_match_no_labels()
    test_neighborhood_hash()
    test_neighborhood_hash_pd()
    test_subgraph_matching()
    test_subgraph_matching_pd()
    test_neighborhood_subgraph_pairwise_distance()
    if cvxopt:
        test_lovasz_theta()
        test_lovasz_theta_pd()
    test_svm_theta()
    test_svm_theta_pd()
    test_odd_sth()
    test_propagation()
    test_propagation_pd()
    test_hadamard_code()
    test_hadamard_code_pd()
    test_multiscale_laplacian()
    test_multiscale_laplacian_pd()
    test_vertex_histogram()
    test_edge_histogram()
    test_graph_hopper()
    test_graph_hopper_pd()
    test_core_framework()
