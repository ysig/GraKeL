"""Tests for the GraphKernel class."""
import argparse

import numpy as np
import numpy.testing as npt

from grakel.graph_kernels import GraphKernel

global verbose, main, development

if __name__ == '__main__':
    # Create an argument parser for the installer of pynauty
    parser = argparse.ArgumentParser(description='A test file for all kernels')

    parser.add_argument('--verbose', help='print kernels with their outputs on\
                        stdout', action="store_true")
    parser.add_argument('--problematic', help='allow execution of problematic\
                        test cases in development', action="store_true")
    parser.add_argument('--ignore_warnings', help='ignore warnings produced by\
                        kernel executions', action="store_true")

    meg = parser.add_mutually_exclusive_group()
    meg.add_argument('--develop', help='execute only tests connected with\
                     current development', action="store_true")
    meg.add_argument('--all', help='execute all tests', action="store_true")
    meg.add_argument('--main', help='execute the main tests [default]',
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

global X, L, k

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

L = {0: 'banana', 1: 'cherry', 2: 'banana', 3: 'cherry',
     4: 'peach', 5: 'cherry', 6: 'lime'}

phi = {0: [0, 1, 2], 1: [0, 0, 1], 2: [1, 1, 0], 3: [3, 0, 1],
       4: [0, 4, 0], 5: [1, 1, 1], 6: [0, 1, 0]}

Le = {(0, 0): 'b', (0, 1): 'c', (0, 2): 'd', (0, 3): 'b', (0, 4): 'c',
      (0, 5): 'a', (0, 6): 'd',
      (1, 0): 'b', (1, 1): 'c', (1, 2): 'd', (1, 3): 'b', (1, 4): 'c',
      (1, 5): 'a', (1, 6): 'd',
      (2, 0): 'd', (2, 1): 'c', (2, 2): 'y', (2, 3): 'b', (2, 4): 'r',
      (2, 5): 'q', (2, 6): 'a',
      (3, 0): 'b', (3, 1): 'a', (3, 2): 'd', (3, 3): 'w', (3, 4): 'c',
      (3, 5): 'q', (3, 6): 'w',
      (4, 0): 'a', (4, 1): 'c', (4, 2): 'u', (4, 3): 'q', (4, 4): 't',
      (4, 5): 'a', (4, 6): 'r',
      (5, 0): 'a', (5, 1): 'y', (5, 2): 'd', (5, 3): 'a', (5, 4): 't',
      (5, 5): 'c', (5, 6): 't',
      (6, 0): 'w', (6, 1): 'c', (6, 2): 'd', (6, 3): 'w', (6, 4): 'r',
      (6, 5): 'c', (6, 6): 'w'}

k = 5


def gk_test_dirac():
    """Test Dirac Kernel inside GraphKernel."""
    XX = list(zip(k*[X], k*[L]))
    gk = GraphKernel(
        verbose=verbose,
        kernel={"name": "dirac"})
    if verbose:
        print("Dirac:", gk.fit_transform(XX))
    else:
        XX_correct = np.full((k, k), 15)
        npt.assert_array_equal(XX_correct, gk.fit_transform(XX))


def gk_test_random_walk_simple():
    """Test Simple Random Walk Kernel inside GraphKernel."""
    XX = list(zip(k*[X], k*[L]))
    gk = GraphKernel(
        verbose=verbose,
        kernel={"name": "random_walk", "lamda": 0.1, "method_type": "simple"})
    if verbose:
        print("Simple:", gk.fit_transform(XX))
    else:
        XX_correct = np.full((k, k), -30.912616526802676)
        npt.assert_array_almost_equal(
            XX_correct, gk.fit_transform(XX), decimal=3)


def gk_test_random_walk_sylvester():
    """Test the Sylvester Random Walk Kernel inside GraphKernel."""
    XX = list(zip(k*[X], k*[L]))
    gk = GraphKernel(
        verbose=verbose,
        kernel={"name": "random_walk",
                "lamda": 0.1, "method_type": "sylvester"})
    if verbose:
        print("Sylvester:", gk.fit_transform(XX))
    else:
        XX_correct = np.full((k, k), -30.912616526802676)
        npt.assert_array_almost_equal(XX_correct,
                                      gk.fit_transform(XX), decimal=3)


def gk_test_shortest_path():
    """Test the Shortest Path Kernel inside GraphKernel."""
    XX = list(zip(k*[X], k*[L]))
    gk = GraphKernel(
        verbose=verbose,
        kernel={"name": "shortest_path", "algorithm_type": "dijkstra"})
    if verbose:
        print("Dijkstra:", gk.fit_transform(XX))
    else:
        XX_correct = np.full((k, k), 66)
        npt.assert_array_equal(XX_correct, gk.fit_transform(XX))

    gk = GraphKernel(
        verbose=verbose,
        kernel={"name": "shortest_path", "algorithm_type": "floyd_warshall"})
    if verbose:
        print("Floyd Warshall:", gk.fit_transform(XX))
    else:
        XX_correct = np.full((k, k), 66)
        npt.assert_array_equal(XX_correct, gk.fit_transform(XX))

    gk = GraphKernel(
        verbose=verbose,
        kernel={"name": "shortest_path", "algorithm_type": "auto"})
    if verbose:
        print("Auto:", gk.fit_transform(XX))
    else:
        XX_correct = np.full((k, k), 66)
        npt.assert_array_equal(XX_correct, gk.fit_transform(XX))


def gk_test_subtree_rg():
    """Test the Subtree Ramon Gartner Kernel inside GraphKernel."""
    XX = list(zip(k*[X], k*[L]))
    gk = GraphKernel(
        verbose=verbose,
        kernel={"name": "subtree_rg", "h": 5})
    if verbose:
        print("Subtree [RG]:", gk.fit_transform(XX))
    else:
        bign = float(
            "154581500920690609988657432454780221337898410611171065400182321" +
            "827306812872340855281323188197412607643176769313264563764886960" +
            "20676136950620929717075828209329592328934916049")
        XX_correct = np.full((k, k), bign)
        npt.assert_array_equal(XX_correct, gk.fit_transform(XX))


def gk_test_graphlet_sampling():
    """Test the Graphlet Sampling kernel inside GraphKernel."""
    XX = list(zip(k*[X], k*[L]))
    gk = GraphKernel(
        verbose=verbose,
        kernel={"name": "graphlets_sampling",
                "k": 5,
                "delta": 0.05,
                "epsilon": 0.05,
                "a": -1})
    if verbose:
        print("Graphlets Sampling:", gk.fit_transform(XX))
    else:
        XX_correct = np.array(
            [[2592593,  2575290,  2596677,  2619840,  2616037],
             [2575290,  2566830,  2583615,  2608620,  2603192],
             [2596677,  2583615,  2607152,  2628438,  2623008],
             [2619840,  2608620,  2628438,  2660086,  2647956],
             [2616037,  2603192,  2623008,  2647956,  2646356]])

        npt.assert_array_almost_equal(
            XX_correct,
            gk.fit_transform(XX),
            decimal=3)


def gk_test_weisfeiler_lehman():
    """Test the Weisfeiler Lehman kernel inside GraphKernel."""
    XX = list(zip(k*[X], k*[L]))
    base_kernel = dict()
    base_kernel["dirac"] = {"name": "dirac"}
    base_kernel["shortest path"] = {"name": "shortest_path"}
    base_kernel["subtree"] = {"name": "subtree_rg"}

    npt_v = dict()
    npt_v["dirac"] = 50
    npt_v["shortest path"] = 276
    npt_v["subtree"] = float(int(
        "154581500920690609988657432454780221337898410611171065400182321" +
        "827306812872340855281323188197412607643176769313264563764886960206" +
        "76136950620929717075828209329592328934916049"))
    for key in base_kernel.keys():
        gk = GraphKernel(
            verbose=verbose,
            kernel=[{"name": "weisfeiler_lehman", "niter": 5},
                    base_kernel[key]])
        if verbose:
            print("Weisfeiler Lehman - "+str(key)+":", gk.fit_transform(XX))
        else:
            XX_correct = np.full((k, k), npt_v[key])
            npt.assert_array_equal(XX_correct, gk.fit_transform(XX))


def gk_test_multiscale_laplacian():
    """Test the Multiscale Laplacian kernel inside GraphKernel."""
    XX = k*[[X, phi]]
    gk = GraphKernel(
        verbose=verbose,
        kernel={"name": "multiscale_laplacian"})

    if verbose:
        print("Multiscale Laplacian:", gk.fit_transform(XX))


def gk_test_subgraph_matching():
    """Test the Subgraph Matching kernel inside GraphKernel."""
    XX = list(k*[[X, L, Le]])
    gk = GraphKernel(
        verbose=verbose,
        kernel={"name": "subgraph_matching"})

    if verbose:
        print("Subgraph Matching:", gk.fit_transform(XX))


def gk_test_lovasz_theta():
    """Test the Lovasz Theta kernel inside GraphKernel."""
    XX = list(k*[[D]])
    gk = GraphKernel(
        verbose=verbose,
        kernel={"name": "lovasz_theta"})

    if verbose:
        print("Lovasz Theta:", gk.fit_transform(XX))


def gk_test_svm_theta():
    """Test the SVM Theta kernel inside GraphKernel."""
    XX = list(k*[[D]])
    gk = GraphKernel(
        verbose=verbose,
        kernel={"name": "svm_theta"})

    if verbose:
        print("SVM Theta:", gk.fit_transform(XX))


def gk_test_neighborhood_subgraph_pairwise_distance():
    """Test the NSPDK inside GraphKernel."""
    XX = list(k*[[X, L, Le]])
    gk = GraphKernel(verbose=verbose,
                     kernel={"name": "neighborhood_subgraph_pairwise_distance"}
                     )

    if verbose:
        print("NSPDK:", gk.fit_transform(XX))


def gk_test_neighborhood_hash():
    """Test the Neighborhood Hash Kernel inside GraphKernel."""
    XX = list(zip(k*[X], k*[L]))
    gk = GraphKernel(verbose=verbose,
                     kernel={"name": "neighborhood_hash", "nh_type": "simple"})

    if verbose:
        print("Neighborhood Hash - 'simple':", gk.fit_transform(XX))

    gk = GraphKernel(verbose=verbose,
                     kernel={"name": "neighborhood_hash",
                             "nh_type": "count-sensitive"})

    if verbose:
        print("Neighborhood Hash - 'count-sensitive':", gk.fit_transform(XX))


def gk_test_odd_sth():
    """Test the ODD-STh Kernel inside GraphKernel."""
    labels = {0: 'A', 1: 'B', 2: 'C', 3: 'D', 4: 'E', 5: 'F', 6: 'G'}

    XX = list(k*[[X, labels]])
    gk = GraphKernel(verbose=verbose,
                     kernel={"name": "odd_sth"})

    if verbose:
        print("ODD-STh:", gk.fit_transform(XX))


def gk_test_wl_nh():
    """Test the Weisfeiler Lehman - Neighborhood Hash Kernel inside GK."""
    XX = list(zip(k*[X], k*[L]))
    base_kernel = {"name": "neighborhood_hash",
                   "nh_type": "count-sensitive"}

    gk = GraphKernel(verbose=verbose,
                     kernel=[{"name": "weisfeiler_lehman",
                              "niter": 5},
                             base_kernel])

    if verbose:
        print("Weisfeiler Lehman - Neighboorhood Hash:", gk.fit_transform(XX))


def gk_test_propagation():
    """Test the Propagation Kernel inside GraphKernel."""
    XX = list(zip(k*[X], k*[L]))
    gk = GraphKernel(verbose=verbose,
                     kernel={"name": "propagation"})

    if verbose:
        print("Propagation:", gk.fit_transform(XX))


def gk_test_pyramid_match():
    """Test the Pyramid Match Kernel inside GraphKernel."""
    XX = list(zip(k*[X], k*[L]))
    gk = GraphKernel(verbose=verbose, kernel={"name": "pyramid_match"})

    if verbose:
        print("pyramid match:", gk.fit_transform(XX))


def gk_test_jsm():
    """Test the Jensen Shannon Representation Kernel inside GraphKernel."""
    XX = list(k*[[D]])
    gk = GraphKernel(verbose=verbose,
                     kernel={"name": "jsm"},
                     concurrency=2,
                     normalize=True)

    if verbose:
        print("Jensen Shannon Representation Kernel:", gk.fit_transform(XX))


if verbose and main:
    gk_test_dirac()
    gk_test_random_walk_simple()
    gk_test_random_walk_sylvester()
    gk_test_shortest_path()
    gk_test_subtree_rg()
    gk_test_graphlet_sampling()
    gk_test_weisfeiler_lehman()
    gk_test_subgraph_matching()
    gk_test_lovasz_theta()
    gk_test_svm_theta()
    gk_test_neighborhood_subgraph_pairwise_distance()
    gk_test_neighborhood_hash()
    gk_test_wl_nh()
    gk_test_odd_sth()
    gk_test_propagation()
    gk_test_pyramid_match()
    gk_test_jsm()

if verbose and develop:
    if problematic:
        gk_test_multiscale_laplacian()
