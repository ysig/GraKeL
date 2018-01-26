"""Tests for the GraphKernel class."""
import argparse

from grakel.dataset import fetch_dataset
from grakel.graph_kernels import GraphKernel

global verbose, main, development

if __name__ == '__main__':
    # Create an argument parser for the installer of pynauty
    parser = argparse.ArgumentParser(description='A test file for all kernels')

    parser.add_argument('--verbose', help='print kernels with their outputs on\
                        stdout', action="store_true")
    parser.add_argument('--problematic', help='allow execution of problematic\
                        test cases in development', action="store_true")
    parser.add_argument('--normalize', help='normalize the kernel output',
                        action="store_true")
    parser.add_argument('--ignore_warnings', help='ignore warnings produced by\
                        kernel executions', action="store_true")
    parser.add_argument(
        '--dataset',
        help='chose the datset you want the tests to be executed',
        type=str,
        default="MUTAG"
    )

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

    normalize = bool(args.normalize)

    dataset_name = args.dataset
else:
    import warnings
    warnings.filterwarnings('ignore', category=UserWarning)
    main, develop, verbose, problematic = True, False, False, False
    normalize = False
    dataset_name = "MUTAG"

global dataset_tr, dataset_te

dataset = fetch_dataset(dataset_name, with_classes=False, verbose=verbose).data
dataset_tr = dataset[:int(len(dataset)*0.8)]
dataset_te = dataset[int(len(dataset)*0.8):]


def test_subtree_wl():
    """Test the wl subtree kernel."""
    gk = GraphKernel(kernel={"name": "subtree_wl"}, verbose=verbose,
                     normalize=normalize)
    if verbose:
        print_kernel_decorator("Subtree WL", gk, dataset_tr, dataset_te)


def test_random_walk():
    """Test the simple random walk kernel."""
    gk = GraphKernel(kernel={"name": "random_walk"}, verbose=verbose,
                     normalize=normalize)
    if verbose:
        print_kernel_decorator("Random Walk", gk, dataset_tr, dataset_te)


def test_shortest_path():
    """Test Shortest Path kernel."""
    gk = GraphKernel(kernel={"name": "shortest_path"}, verbose=verbose,
                     normalize=normalize)
    if verbose:
        print_kernel_decorator("Shortest Path", gk, dataset_tr, dataset_te)


def test_graphlet_sampling():
    """Test the Graphlet Sampling Kernel."""
    gk = GraphKernel(kernel={"name": "graphlet_sampling", "n_samples": 200},
                     verbose=verbose, normalize=normalize)
    if verbose:
        print_kernel_decorator("Graphlet Sampling", gk, dataset_tr, dataset_te)


def test_weisfeiler_lehman():
    """Test the Weisfeiler Lehman kernel."""
    gk = GraphKernel(kernel=[{"name": "weisfeiler_lehman"},
                     {"name": "subtree_wl"}],
                     verbose=verbose, normalize=normalize)
    if verbose:
        print_kernel_decorator("WL/Subtree", gk, dataset_tr, dataset_te)


def test_pyramid_match():
    """Test Pyramid Match kernel."""
    gk = GraphKernel(kernel={"name": "pyramid_match"}, verbose=verbose,
                     normalize=normalize)
    if verbose:
        print_kernel_decorator("Pyramid Match", gk, dataset_tr, dataset_te)


def test_neighborhood_hash():
    """Test Neighborhood Hash kernel."""
    gk = GraphKernel(kernel={"name": "neighborhood_hash"}, verbose=verbose,
                     normalize=normalize)
    if verbose:
        print_kernel_decorator("Neighborhood Hash", gk, dataset_tr, dataset_te)


def print_kernel_decorator(name, kernel, X, Y):
    """Print kernels in case of verbose execution."""
    name += " [decorator]"
    print(str(name) + ":\n" + (len(str(name)) * "-") + "-\n")
    print("fit_transform\n-------------")
    print(kernel.fit_transform(X))
    print("\ntransform\n---------")
    print(kernel.transform(Y))
    print("--------------------------------------" +
          "--------------------------------------\n")


if verbose and main:
    test_subtree_wl()
    test_random_walk()
    test_shortest_path()
    test_weisfeiler_lehman()
    test_pyramid_match()
    test_neighborhood_hash()
    test_graphlet_sampling()

if verbose and develop:
    if problematic:
        pass
