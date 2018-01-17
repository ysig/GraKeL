"""Tests for the kernel sub-module."""
import numpy as np
import numpy.testing as npt

from grakel.kernels import dirac
from grakel.kernels import dirac_pair
from grakel.kernels import graphlet_sampling
from grakel.kernels import hadamard_code
from grakel.kernels import jsm
from grakel.kernels import lovasz_theta
from grakel.kernels import multiscale_laplacian
from grakel.kernels import neighborhood_hash
from grakel.kernels import neighborhood_subgraph_pairwise_distance
from grakel.kernels import odd_sth
from grakel.kernels import propagation
from grakel.kernels import pyramid_match
from grakel.kernels import random_walk
from grakel.kernels import shortest_path
from grakel.kernels import shortest_path_matrix
from grakel.kernels import subgraph_matching
from grakel.kernels import subtree_rg
from grakel.kernels import subtree_rg_pair
from grakel.kernels import svm_theta
from grakel.kernels import weisfeiler_lehman

global verbose, main, development

if __name__ == '__main__':
    import argparse
    # Create an argument parser for the installer of pynauty
    parser = argparse.ArgumentParser(description='A test file for all kernels')

    parser.add_argument(
        '--verbose',
        help='print kernels with their outputs on stdout',
        action="store_true")
    parser.add_argument(
        '--problematic',
        help='allow execution of problematic test cases in development',
        action="store_true")
    parser.add_argument(
        '--ignore_warnings',
        help='ignore warnings produced by kernel executions',
        action="store_true")

    meg = parser.add_mutually_exclusive_group()
    meg.add_argument(
        '--develop',
        help='execute only tests connected with current development',
        action="store_true")
    meg.add_argument(
        '--all',
        help='execute all tests',
        action="store_true")
    meg.add_argument(
        '--main',
        help='execute the main tests [default]',
        action="store_true")

    args = parser.parse_args()

    verbose = bool(args.verbose)
    if args.all:
        main, develop = True, True
    elif args.develop:
        main, develop = False, True
    else:
        main, develop = True, False
    problematic = bool(args.problematic)

    if bool(args.ignore_warnings):
        import warnings
        warnings.filterwarnings('ignore', category=UserWarning)
else:
    import warnings
    warnings.filterwarnings('ignore', category=UserWarning)
    main, develop, verbose = True, False, False

global X, Y, L, Le, phi

X = np.array([[0, 1, 2, 1, 0, 0.5, 0],
              [1, 0, 0, 0, 1, 2, 0.5],
              [2, 0, 0, 3, 0, 0, 2],
              [1, 0, 3, 0, 0, 0, 0],
              [0, 1, 0, 0, 0, 3, 1],
              [0.5, 2, 0, 0, 3, 0, 0],
              [0, 0.5, 2, 0, 1, 0, 0]])

# Route May 2001 map of Spirit Airlines
# Atlantic City, Chicago (O'Hare), Detroit, Fort Lauderdale, Fort Myers,
# Los Angeles, Melbourne, Myrtle Beach, Newark, New York (LaGuardia), Oakland,
# Orlando, Tampa, and Washington (Reagan National).
D = np.array([
    [0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0],
    [1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0],
    [0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0],
    [1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1],
    [1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
    [0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1],
    [1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0],
    [0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
    [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]])

LD = {0: 'rock', 1: 'post-punk', 2: 'rock', 3: 'indie', 4: 'classical',
      5: 'pop', 6: 'rock', 7: 'punk', 8: 'punk', 9: 'indie', 10: 'post-rock',
      11: 'post-punk', 12: 'jazz', 13: 'jazz', 14: 'jazz', 15: 'classical'}

L = {0: 'banana', 1: 'cherry', 2: 'banana',
     3: 'cherry', 4: 'peach', 5: 'cherry', 6: 'lime'}

phi = {0: [0, 1, 2], 1: [0, 0, 1], 2: [1, 1, 0], 3: [3, 0, 1],
       4: [0, 4, 0], 5: [1, 1, 1], 6: [0, 1,  0]}

Le = {(0, 0): 'b', (0, 1): 'c', (0, 2): 'd', (0, 3): 'b',
      (0, 4): 'c', (0, 5): 'a', (0, 6): 'd',
      (1, 0): 'b', (1, 1): 'c', (1, 2): 'd',
      (1, 3): 'b', (1, 4): 'c', (1, 5): 'a', (1, 6): 'd',
      (2, 0): 'd', (2, 1): 'c', (2, 2): 'y', (2, 3): 'b',
      (2, 4): 'r', (2, 5): 'q', (2, 6): 'a',
      (3, 0): 'b', (3, 1): 'a', (3, 2): 'd',
      (3, 3): 'w', (3, 4): 'c', (3, 5): 'q', (3, 6): 'w',
      (4, 0): 'a', (4, 1): 'c', (4, 2): 'u', (4, 3): 'q',
      (4, 4): 't', (4, 5): 'a', (4, 6): 'r',
      (5, 0): 'a', (5, 1): 'y', (5, 2): 'd', (5, 3): 'a',
      (5, 4): 't', (5, 5): 'c', (5, 6): 't',
      (6, 0): 'w', (6, 1): 'c', (6, 2): 'd', (6, 3): 'w',
      (6, 4): 'r', (6, 5): 'c', (6, 6): 'w'}


def test_dirac():
    """Test the dirac kernel."""
    if verbose:
        print("Dirac:", dirac(X, X, L, L))
    else:
        npt.assert_equal(15, dirac(X, X, L, L))


def test_random_walk_simple():
    """Test the simple random walk kernel."""
    if verbose:
        print("Random Walk [Simple]:",
              random_walk(X, X, lamda=0.1, method_type="simple"))
    else:
        npt.assert_almost_equal(
            -30.912616526802676,
            random_walk(X, X, lamda=0.1, method_type="simple"),
            decimal=3)


def test_random_walk_sylvester():
    """Test the sylvester random walk kernel."""
    if verbose:
        print("Random Walk [Sylvester]:",
              random_walk(X, X, lamda=0.1, method_type="sylvester"))
    else:
        npt.assert_almost_equal(
            -30.912616526802676,
            random_walk(X, X, lamda=0.1, method_type="sylvester"),
            decimal=3)


def test_shortest_path():
    """Test Shortest Path kernel."""
    if verbose:
        print("Shortest Path [Dijkstra]:",
              shortest_path(X, X, L, L, "dijkstra"))
        print("Shortest Path [Floyd Warshall]:",
              shortest_path(X, X, L, L, "floyd_warshall"))
        print("Shortest Path [Auto]:",
              shortest_path(X, X, L, L, "auto"))
        print("Shortest path Attributes [Dijkstra]:",
              shortest_path(X, X, phi, phi, as_attributes=True))
    else:
        npt.assert_equal(66, shortest_path(X, X, L, L, "dijkstra"))
        npt.assert_equal(66, shortest_path(X, X, L, L, "floyd_warshall"))
        npt.assert_equal(66, shortest_path(X, X, L, L, "auto"))


def test_subtree_RG():
    """Test the Subtree Ramon Gartner Kernel."""
    if verbose:
        print("Subtree [RG]:", subtree_rg(X, X, L, L, 2))
    else:
        bign = int(
            "154581500920690609988657432454780221337898410611171065400182321" +
            "827306812872340855281323188197412607643176769313264563764886960" +
            "20676136950620929717075828209329592328934916049")
        npt.assert_equal(bign, subtree_rg(X, X, L, L, 2))


def test_graphlet_sampling():
    """Test the Graphlet Sampling Kernel."""
    if verbose:
        print("Graphlets Sampling:",
              graphlet_sampling(X, X, k=5, delta=0.05, epsilon=0.05, a=-1))
    else:
        npt.assert_almost_equal(
            2575290,
            graphlet_sampling(X, X, 5, 0.05, 0.05, -1),
            decimal=3)


def test_weisfeiler_lehman():
    """Test the Weisfeiler Lehman kernel."""
    base_kernel = dict()
    base_kernel["dirac"] = lambda x, y: dirac_pair(x, y)
    base_kernel["shortest path"] = lambda x, y:\
        shortest_path_matrix({0: x}, {0: y})[0, 0]
    base_kernel["subtree"] = lambda x, y: subtree_rg_pair(x, y)

    if not verbose:
        npt_v = dict()
        npt_v["dirac"] = 50
        npt_v["shortest path"] = 276
        npt_v["subtree"] = float(
            "154581500920690609988657432454780221337898410611171065400182321" +
            "827306812872340855281323188197412607643176769313264563764886960" +
            "20676136950620929717075828209329592328934916049")
    for k in base_kernel.keys():
        if verbose:
            print("Weisfeiler_lehman - "+str(k)+":",
                  weisfeiler_lehman(X, X, L, L, base_kernel[k], 5))
        else:
            npt.assert_equal(npt_v[k],
                             weisfeiler_lehman(X, X, L, L, base_kernel[k], 5))


def test_multiscale_laplacian():
    """Test Multiscale Laplacian kernel."""
    if verbose:
        print("Multiscale Laplacian:", multiscale_laplacian(X, X, phi, phi))


def test_subgraph_matching():
    """Test Subgraph Matching kernel."""
    if verbose:
        print("Subgraph Matching:", subgraph_matching(X, X, L, L, Le, Le))


def test_lovasz_theta():
    """Test Lovasz-Theta kernel."""
    if verbose:
        print("Lovasz Theta:", lovasz_theta(D, D))


def test_svm_theta():
    """Test the SVM-theta kernel."""
    if verbose:
        print("SVM Theta:", svm_theta(D, D))


def test_neighborhood_subgraph_pairwise_distance_kernel():
    """Test the NSPDK-theta kernel."""
    if verbose:
        print("NSPDK:",
              neighborhood_subgraph_pairwise_distance(X, X, L, L, Le, Le))


def test_neighborhood_hash_kernel():
    """Test the neighborhood hash kernel."""
    if verbose:
        print("Neighborhood Hash - 'simple':",
              neighborhood_hash(X, X, L, L, nh_type='simple'))
        print("Neighborhood Hash - 'count-sensitive':",
              neighborhood_hash(X, X, L, L, nh_type='count-sensitive'))


def test_odd_sth():
    """Test the ODD-STh kernel."""
    # x = np.array([[0, 1, 1, 0], [1, 0, 1, 1], [1, 1, 0, 1], [0, 1, 1, 0]])
    # l = {0: 's', 1:'e', 2:'b', 3:'d'}
    if verbose:
        print("ODD-STh:", odd_sth(X, X, L, L, h=None))


def test_propagation():
    """Test the Propagation kernel."""
    if verbose:
        print("Propagation:", propagation(X, X, L, L))


def test_pyramid_match():
    """Test the Pyramid Match Kernel."""
    if verbose:
        print("Pyramid Match:", pyramid_match(D, D, LD, LD))


def test_hadamard_code():
    """Test the Hadamard Code kernel."""
    if verbose:
        print("Hadamard Code [simple]:",
              hadamard_code(X, X, L, L, dirac_pair))
        print("Hadamard Code [shortened]:",
              hadamard_code(X, X, L, L, dirac_pair, hc_type='shortened'))


def test_jsm():
    """Test the Jensen Shannon Representation Kernel."""
    if verbose:
        print("Jensen Shannon Representation Kernel:", jsm(D, D))


if verbose and main:
    test_dirac()
    test_random_walk_simple()
    test_random_walk_sylvester()
    test_shortest_path()
    test_subtree_RG()
    test_graphlet_sampling()
    test_weisfeiler_lehman()
    test_subgraph_matching()
    test_lovasz_theta()
    test_svm_theta()
    test_neighborhood_subgraph_pairwise_distance_kernel()
    test_neighborhood_hash_kernel()
    test_odd_sth()
    test_propagation()
    test_pyramid_match()
    test_hadamard_code()
    test_jsm()
if verbose and develop:
    if problematic:
        test_multiscale_laplacian()
