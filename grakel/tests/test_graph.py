"""Tests for the Graph class."""
# Author: Ioannis Siglidis <y.siglidis@gmail.com>
# License: BSD 3 clause
import numpy as np
import numpy.testing as npt

from grakel.graph import Graph

global verbose

# Add extra arguments for allowing unit testing
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description='A test file for all `Graph` type objects')
    parser.add_argument(
        '--verbose',
        help='verbose outputs on stdout',
        action="store_true")
    parser.add_argument(
        '--ignore_warnings',
        help='ignore warnings produced by kernel executions',
        action="store_true")

    args = parser.parse_args()
    verbose = bool(args.verbose)

    if bool(args.ignore_warnings):
        import warnings
        warnings.filterwarnings('ignore', category=UserWarning)
else:
    import warnings
    warnings.filterwarnings('ignore', category=UserWarning)
    verbose = False


def test_graph_adjacency():
    """Testing Graph object consistency for an adjacency-type initialization object."""
    # Input
    X = np.array([[1, 1, 0, 3], [1, 0, 0, 2], [2, 3, 0, 1], [1, 0, 0, 0]])
    labels = {0: 'banana', 1: 'cherry', 2: 'banana', 3: 'cherry'}

    # try all formats
    g = dict()
    g["auto"] = Graph(X, labels, {}, "auto")
    g["dict"] = Graph(X, labels, {}, "dictionary")
    g["adjc"] = Graph(X, labels, {}, "adjacency")
    g["all"] = Graph(X, labels, {}, "all")

    # Desired output label group
    desired_output_label_group = {'cherry': [1, 3], 'banana': [0, 2]}

    for k in g.keys():
        gklg = g[k].get_label_group()
        if verbose:
            print(k)
            print(gklg, '\n')
        else:
            npt.assert_equal(desired_output_label_group, gklg)

    # Desired Shortest path matrix
    spm_do = [[0., 1., float("Inf"), 3.],
              [1., 0., float("Inf"), 2.],
              [2., 3., 0., 1.],
              [1., 2., float("Inf"), 0.]]

    for k in g.keys():
        spm, spl = g[k].build_shortest_path_matrix(algorithm_type="auto")
        if verbose:
            print(k)
            print(spm, '\n', spl, '\n')
        else:
            npt.assert_array_equal(spm, spm_do)
            npt.assert_equal(spl, labels)


def test_graph_edge_dictionary():
    """Testing Graph object consistency for an edge-dictionary-type initialization object."""
    # Input
    X = {'a': {'a': 1, 'b': 1, 'd': 3},
         'b': {'a': 1, 'd': 2},
         'c': {'a': 2, 'b': 3, 'd': 1},
         'd': {'a': 1}}

    labels = {'a': 'banana', 'b': 'cherry', 'c': 'banana', 'd': 'cherry'}

    # Test for all Graph formats
    g = dict()
    g["auto"] = Graph(X, labels, {}, "auto")
    g["dict"] = Graph(X, labels, {}, "dictionary")
    g["adjc"] = Graph(X, labels, {}, "adjacency")
    g["all"] = Graph(X, labels, {}, "all")

    # Desired output label group
    desired_output_label_group = {'cherry': set(['d', 'b']),
                                  'banana': set(['a', 'c'])}
    desired_output_label_group_idx = {'banana': set([0, 2]),
                                      'cherry': set([1, 3])}

    def proper_dict(x):
        return {key: set(x[key]) for key in x.keys()}

    for k in g.keys():
        gklg = g[k].get_label_group()
        if verbose:
            print(k)
            print(gklg, '\n')
        else:
            if (k == "adjc"):
                npt.assert_equal(
                    desired_output_label_group_idx,
                    proper_dict(gklg))
            else:
                npt.assert_equal(
                    desired_output_label_group,
                    proper_dict(gklg))

    # Desired Shortest path matrix
    spm_do = [[0., 1., float("Inf"), 3.],
              [1., 0., float("Inf"), 2.],
              [2., 3., 0., 1.],
              [1., 2., float("Inf"), 0.]]

    desired_labels = {0: 'banana', 1: 'cherry', 2: 'banana', 3: 'cherry'}

    for k in g.keys():
        spm, spl = g[k].build_shortest_path_matrix(algorithm_type="auto")
        if verbose:
            print(k)
            print(spm, '\n', spl, '\n')
        else:
            npt.assert_array_equal(spm, spm_do)
            npt.assert_equal(spl, desired_labels)


if __name__ == '__main__':
    test_graph_adjacency()
    test_graph_edge_dictionary()
