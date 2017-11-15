import numpy as np
from numpy.testing import assert_almost_equal

from kernels import *
X = np.array([[0,1,2,1,0,0.5,0],
              [1,0,0,0,1,2,0.5],
              [2,0,0,3,0,0,2],
              [1,0,3,0,0,0,0],
              [0,1,0,0,0,3,1],
              [0.5,2,0,0,3,0,0],
              [0,0.5,2,0,1,0,0]])
L = {0:'banana', 1:'cherry', 2:'banana', 3:'cherry', 4:'peach',5:'cherry',6:'lime'}

def test_dirac()
    print("Dirac:",dirac(X, X, L, L))
def test_random_walk_simple():
    print("Random Walk Simple:",random_walk(X, X, L, L, lamda=0.1, method_type="simple"))

def test_random_walk_sylvester():
    print("Random Walk Sylvester:",random_walk(X, X, L, L, lamda=0.1, method_type="sylvester"))

def test_shortest_path():
    print("Dijkstra:",shortest_path(X, X, L, L, "dijkstra"))
    print("Floyd Warshall:",(shortest_path(X, X, L, L, "floyd_warshall")))
    print("Auto:",shortest_path(X, X, L, L, "auto"))

def test_subtree_RG():
    print("Subtree [RG]:",subtree_RG(X,X,L,L,3))

def test_graphlets_sampling():
    print("Graphlets Sampling:",graphlets_sampling(X, X, 5, 0.05, 0.05,-1))

def test_weisfeiler_lehman():
    print("Weisfeiler_lehman:")

print("Dirac\n---------------------\n")
test_dirac()
print("Random Walk\n------------------\n")
test_random_walk_simple()
test_random_walk_sylvester()
print("\nShortest Path\n-----------------\n")
test_shortest_path()
print("\nSubtree\n------------------\n")
test_subtree_RG()
print("\nGraphlet\n------------------\n")
test_graphlets_sampling()
print("\nWeisfehler Lehman\n-----------------\n")

