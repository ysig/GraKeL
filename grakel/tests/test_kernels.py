"""Tests for the kernel sub-module."""
from grakel.dataset import load_dataset

from grakel.kernels import graphlet_sampling
from grakel.kernels import random_walk
from grakel.kernels import subtree_wl
from grakel.kernels import shortest_path
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

global MUTAG_tr, MUTAG_te

MUTAG = load_dataset("MUTAG", with_classes=False, verbose=verbose)
MUTAG_tr = MUTAG[:int(len(MUTAG)*0.8)]
MUTAG_te = MUTAG[int(len(MUTAG)*0.8):]


def test_subtree_wl():
    """Test the wl subtree kernel."""
    stwl_kernel = subtree_wl(verbose=verbose)
    if verbose:
        print_kernel("Subtree WL", stwl_kernel, MUTAG_tr, MUTAG_te)


def test_random_walk():
    """Test the simple random walk kernel."""
    rw_kernel = random_walk(verbose=verbose)
    if verbose:
        print_kernel("Random Walk", rw_kernel, MUTAG_tr, MUTAG_te)


def test_shortest_path():
    """Test Shortest Path kernel."""
    sp_kernel = shortest_path(verbose=verbose)
    if verbose:
        print_kernel("Shortest Path", sp_kernel, MUTAG_tr, MUTAG_te)


def test_graphlet_sampling():
    """Test the Graphlet Sampling Kernel."""
    gs_kernel = graphlet_sampling(verbose=verbose, n_samples=200)
    if verbose:
        print_kernel("Graphlet Sampling", gs_kernel, MUTAG_tr, MUTAG_te)


def test_weisfeiler_lehman():
    """Test the Weisfeiler Lehman kernel."""
    wl_st_kernel = weisfeiler_lehman(base_kernel=subtree_wl)
    if verbose:
        print_kernel("WL/Subtree", wl_st_kernel, MUTAG_tr, MUTAG_te)


def print_kernel(name, kernel, X, Y):
    """Print kernels in case of verbose execution."""
    print("\n" + str(name) + ":\n" + (len(str(name)) * "-") + "-")
    print("\nfit_transform\n-------------")
    print(kernel.fit_transform(X))
    print("\ntransform\n---------")
    print(kernel.transform(Y))
    print("---------------------------------------------------------------")


if verbose and main:
    test_subtree_wl()
    test_random_walk()
    test_shortest_path()
    test_graphlet_sampling()
    test_weisfeiler_lehman()
if verbose and develop:
    if problematic:
        pass
