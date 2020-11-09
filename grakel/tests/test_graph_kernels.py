"""Tests for the GraphKernel class."""
# Author: Ioannis Siglidis <y.siglidis@gmail.com>
# License: BSD 3 clause
import os
from time import time
from warnings import warn

from sklearn.model_selection import train_test_split

from grakel.datasets import fetch_dataset
from grakel.datasets.base import read_data
from grakel.datasets import get_dataset_info
from grakel.graph_kernels import GraphKernel


fdir = os.path.dirname(__file__)

global verbose, main, development

if __name__ == '__main__':
    import argparse

    # Create an argument parser for the installer of pynauty
    parser = argparse.ArgumentParser(description='A test file for all kernels')

    parser.add_argument('--verbose', help='print kernels with their outputs' +
                        ' on stdout', action="store_true")
    parser.add_argument('--time', help='time the kernel computation (has ef' +
                        'fect only on verbose)', action="store_true")
    parser.add_argument('--problematic', help='allow execution of problemati' +
                        'c test cases in development', action="store_true")
    parser.add_argument('--slow', help='allow execution of slow test cases' +
                        ' in development', action="store_true")
    parser.add_argument('--normalize', help='normalize the kernel output',
                        action="store_true")
    parser.add_argument('--ignore_warnings', help='ignore warnings produced ' +
                        'by kernel executions', action="store_true")

    parser.add_argument(
        '--dataset',
        help='choose the dataset for tests requiring node/edge labels',
        type=str,
        default="MUTAG"
    )

    parser.add_argument(
        '--dataset_attr',
        help='choose the dataset for tests requiring node attributes',
        type=str,
        default="Cuneiform"
    )

    meg = parser.add_mutually_exclusive_group()
    meg.add_argument('--develop', help='execute only tests connected with ' +
                     'current development', action="store_true")
    meg.add_argument('--all', help='execute all tests', action="store_true")
    meg.add_argument('--main', help='execute the main tests [default]',
                     action="store_true")

    args = parser.parse_args()

    verbose = bool(args.verbose)
    time_kernel = bool(args.time)
    if args.all:
        main, develop = True, True
    elif args.develop:
        main, develop = False, True
    else:
        main, develop = True, False
    problematic = bool(args.problematic)
    slow = bool(args.slow)

    if bool(args.ignore_warnings):
        import warnings
        warnings.filterwarnings('ignore', category=UserWarning)

    normalize = bool(args.normalize)

    dataset_name = args.dataset
    dataset_attr_name = args.dataset_attr

else:
    import warnings
    warnings.filterwarnings('ignore', category=UserWarning)
    main, develop, problematic, slow = True, False, False, False
    normalize, verbose, time_kernel = False, False, False
    dataset_name = "MUTAG"
    dataset_attr_name = "Cuneiform"

global dataset_tr, dataset_te, dataset_attr_tr, dataset_attr_te

# consistency check for the dataset
dinfo = get_dataset_info(dataset_name)
if dinfo is None:
    raise TypeError('dataset not found')
elif not dinfo["nl"] and not dinfo["el"]:
    raise TypeError('dataset must have either node and edge labels')

# consistency check for the attribute dataset
dinfo_attr = get_dataset_info(dataset_attr_name)
if dinfo is None:
    raise TypeError('dataset for attributes not found')
elif not dinfo_attr["nl"] and not dinfo_attr["el"]:
    raise TypeError('dataset must have node attributes')


# The baseline dataset for node, edge_labels
global dataset, dataset_tr, dataset_te

try:
    dataset = fetch_dataset(dataset_name, with_classes=False, verbose=verbose).data
except Exception:
    # Offline testing
    warn('There was a problem fetching dataset for attributes: [' + dataset_name + ']')
    if dataset_name != 'MUTAG':
        warn('Switching back to baseline dataset MUTAG')
    warn('Using an offline version..')
    cwd = os.getcwd()
    os.chdir(os.path.join(fdir, 'data'))
    dataset = read_data('MUTAG', with_classes=False).data
    os.chdir(cwd)

dataset_tr, dataset_te = train_test_split(dataset,
                                          test_size=0.2,
                                          random_state=42)

# The baseline dataset for node/edge-attributes
global dataset_attr, dataset_attr_tr, dataset_attr_te

try:
    dataset_attr = fetch_dataset(dataset_attr_name, with_classes=False,
                                 prefer_attr_nodes=True,
                                 prefer_attr_edges=True,
                                 verbose=verbose).data
except Exception:
    # Offline testing
    warn('There was a problem fetching dataset for attributes: [' + dataset_attr_name + ']')
    if dataset_attr_name != 'Cuneiform':
        warn('Switching back to baseline dataset Cuneiform')
    warn('Using an offline version..')
    cwd = os.getcwd()
    os.chdir(os.path.join(fdir, 'data'))
    dataset_attr = read_data('Cuneiform', with_classes=False, prefer_attr_nodes=True,
                             prefer_attr_edges=True).data
    os.chdir(cwd)

dataset_attr_tr, dataset_attr_te = train_test_split(dataset_attr,
                                                    test_size=0.2,
                                                    random_state=42)

if verbose:
    # Guard for nosetests
    def test_random_walk():
        """Test the Simple random walk kernel [decorator]."""
        gk = GraphKernel(kernel={"name": "random_walk"}, verbose=verbose,
                         normalize=normalize)
        if verbose:
            print_kernel_decorator("Random Walk", gk, dataset_tr, dataset_te)

    def test_shortest_path():
        """Test the Shortest Path kernel [decorator]."""
        gk = GraphKernel(kernel={"name": "shortest_path"}, verbose=verbose,
                         normalize=normalize)
        if verbose:
            print_kernel_decorator("Shortest Path", gk, dataset_tr, dataset_te)

    def test_graphlet_sampling():
        """Test the Graphlet Sampling Kernel [decorator]."""
        gk = GraphKernel(kernel={"name": "graphlet_sampling",
                                 "sampling": {"n_samples": 200}},
                         verbose=verbose, normalize=normalize)

        if verbose:
            print_kernel_decorator("Graphlet Sampling", gk, dataset_tr, dataset_te)

    def test_weisfeiler_lehman():
        """Test the Weisfeiler Lehman kernel [decorator]."""
        gk = GraphKernel(kernel=[{"name": "weisfeiler_lehman"},
                         {"name": "vertex_histogram"}],
                         verbose=verbose, normalize=normalize)
        if verbose:
            print_kernel_decorator("WL/Subtree", gk, dataset_tr, dataset_te)

    def test_weisfeiler_lehman_optimal_assignment():
        """Test the Weisfeiler Lehman Optimal Assignment kernel [decorator]."""
        gk = GraphKernel(kernel="WL-OA", verbose=verbose, normalize=normalize)
        if verbose:
            print_kernel_decorator("WL-OA", gk, dataset_tr, dataset_te)

    def test_pyramid_match():
        """Test the Pyramid Match kernel [decorator]."""
        gk = GraphKernel(kernel={"name": "pyramid_match"}, verbose=verbose,
                         normalize=normalize)
        if verbose:
            print_kernel_decorator("Pyramid Match", gk, dataset_tr, dataset_te)

    def test_neighborhood_hash():
        """Test the Neighborhood Hash kernel [decorator]."""
        gk = GraphKernel(kernel={"name": "neighborhood_hash"}, verbose=verbose,
                         normalize=normalize)
        if verbose:
            print_kernel_decorator("Neighborhood Hash", gk, dataset_tr, dataset_te)

    def test_subgraph_matching():
        """Test the Subgraph Matching kernel [decorator]."""
        gk = GraphKernel(kernel={"name": "subgraph_matching"}, verbose=verbose,
                         normalize=normalize)
        if verbose:
            print_kernel_decorator("Subgraph Matching", gk, dataset_tr, dataset_te)

    def test_neighborhood_pairwise_distance():
        """Test the Neighborhood Subgraph Pairwise Distance kernel [decorator]."""
        gk = GraphKernel(kernel={
            "name": "neighborhood_subgraph_pairwise_distance"},
                         verbose=verbose, normalize=normalize)

        if verbose:
            print_kernel_decorator("NSPD", gk, dataset_tr, dataset_te)

    def test_lovasz_theta():
        """Test the Lovasz-theta kernel [decorator]."""
        try:
            gk = GraphKernel(kernel={"name": "lovasz_theta"},
                             verbose=verbose, normalize=normalize)
        except ImportError:
            return

        if verbose:
            print_kernel_decorator("Lovasz-theta", gk, dataset_tr, dataset_te)

    def test_svm_theta():
        """Test the SVM-theta kernel [decorator]."""
        gk = GraphKernel(kernel={"name": "svm_theta"},
                         verbose=verbose, normalize=normalize)

        if verbose:
            print_kernel_decorator("SVM-theta", gk, dataset_tr, dataset_te)

    def test_odd_sth():
        """Test the ODD-STh kernel [decorator]."""
        gk = GraphKernel(kernel={"name": "odd_sth"},
                         verbose=verbose, normalize=normalize)

        if verbose:
            print_kernel_decorator("ODD-STh", gk, dataset_tr, dataset_te)

    def test_propagation():
        """Test the Propagation kernel [decorator]."""
        gk = GraphKernel(kernel={"name": "propagation"},
                         verbose=verbose, normalize=normalize)

        if verbose:
            print_kernel_decorator("Propagation", gk, dataset_tr, dataset_te)

    def test_hadamard_code():
        """Test the Hadamard Code kernel [decorator]."""
        gk = GraphKernel(kernel=[{"name": "hadamard_code"},
                                 {"name": "subtree_wl"}],
                         verbose=verbose, normalize=normalize)

        if verbose:
            print_kernel_decorator("Hadamard-Code/Subtree-WL [Simple]",
                                   gk, dataset_tr, dataset_te)

    def test_edge_histogram():
        """Test the Edge Histogram kernel [decorator]."""
        gk = GraphKernel(kernel={"name": "edge_histogram"},
                         verbose=verbose, normalize=normalize)

        if verbose:
            print_kernel_decorator("Edge Histogram",
                                   gk, dataset_tr, dataset_te)

    def test_vertex_histogram():
        """Test the Vertex Histogram kernel [decorator]."""
        gk = GraphKernel(kernel={"name": "vertex_histogram"},
                         verbose=verbose, normalize=normalize)

        if verbose:
            print_kernel_decorator("Vertex Histogram",
                                   gk, dataset_tr, dataset_te)

    def test_multiscale_laplacian():
        """Test the Multiscale Laplacian kernel [decorator]."""
        gk = GraphKernel(kernel="ML", verbose=verbose, normalize=normalize)

        if verbose:
            print_kernel_decorator("Multiscale Laplacian",
                                   gk, dataset_attr_tr, dataset_attr_te)

    def test_graph_hopper():
        """Test the Graph Hopper kernel [decorator]."""
        gk = GraphKernel(kernel={"name": "graph_hopper"},
                         verbose=verbose, normalize=normalize)

        if verbose:
            print_kernel_decorator("Graph Hopper",
                                   gk, dataset_attr_tr, dataset_attr_te)

    def test_core_framework():
        """Test the Core Framework kernel [decorator]."""
        kernel = [{"name": "core_framework"}, {"name": "weisfeiler_lehman"},
                  {"name": "vertex_histogram"}]
        gk = GraphKernel(kernel=kernel, verbose=verbose, normalize=normalize)
        if verbose:
            print_kernel_decorator("Core Framework",
                                   gk, dataset_tr, dataset_te)

    def sec_to_time(sec):
        """Print time in a correct format."""
        dt = list()
        days = int(sec // 86400)
        if days > 0:
            sec -= 86400*days
            dt.append(str(days) + " d")

        hrs = int(sec // 3600)
        if hrs > 0:
            sec -= 3600*hrs
            dt.append(str(hrs) + " h")

        mins = int(sec // 60)
        if mins > 0:
            sec -= 60*mins
            dt.append(str(mins) + " m")

        if sec > 0:
            dt.append(str(round(sec, 2)) + " s")
        return " ".join(dt)

    def print_kernel_decorator(name, kernel, X, Y):
        """Print kernels in case of verbose execution."""
        if time_kernel:
            name += " [decorator]"
            print(str(name) + ":\n" + (len(str(name)) * "-") + "-\n")
            print("fit_transform\n-------------")

            # [time] fit_transform
            start = time()
            Kft = kernel.fit_transform(X)
            ft_time = time() - start

            print(Kft)
            print("[TIME] fit_transform:", sec_to_time(ft_time))
            print("\ntransform\n---------")

            start = time()
            Kt = kernel.transform(Y)
            t_time = time() - start
            print(Kt)
            print("[TIME] transform:", sec_to_time(t_time))
            print("[TIME] total:", sec_to_time(ft_time+t_time))
            print("--------------------------------------" +
                  "--------------------------------------\n")
        else:
            name += " [decorator]"
            print(str(name) + ":\n" + (len(str(name)) * "-") + "-\n")
            print("fit_transform\n-------------")
            print(kernel.fit_transform(X))

            print("\ntransform\n---------")
            print(kernel.transform(Y))
            print("--------------------------------------" +
                  "--------------------------------------\n")

if verbose and main:
    test_random_walk()
    test_shortest_path()
    test_weisfeiler_lehman()
    test_weisfeiler_lehman_optimal_assignment()
    test_pyramid_match()
    test_neighborhood_hash()
    test_graphlet_sampling()
    test_lovasz_theta()
    test_svm_theta()
    test_odd_sth()
    test_propagation()
    test_hadamard_code()
    test_neighborhood_pairwise_distance()
    test_multiscale_laplacian()
    test_subgraph_matching()
    test_vertex_histogram()
    test_edge_histogram()
    test_core_framework()
    test_graph_hopper()
