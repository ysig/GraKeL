import numpy as np
from numpy.testing import assert_almost_equal

from grakel.kernels import *

X = np.array([[0,1,2,1,0,0.5,0],
              [1,0,0,0,1,2,0.5],
              [2,0,0,3,0,0,2],
              [1,0,3,0,0,0,0],
              [0,1,0,0,0,3,1],
              [0.5,2,0,0,3,0,0],
              [0,0.5,2,0,1,0,0]])
L = {0:'banana', 1:'cherry', 2:'banana', 3:'cherry', 4:'peach',5:'cherry',6:'lime'}

def test_dirac():
    print("Dirac:",dirac(X, X, L, L))

def test_random_walk_simple():
    print("Simple:",random_walk(X, X, lamda=0.1, method_type="simple"))

def test_random_walk_sylvester():
    print("Sylvester:",random_walk(X, X, lamda=0.1, method_type="sylvester"))

def test_shortest_path():
    print("Dijkstra:",shortest_path(X, X, L, L, "dijkstra"))
    print("Floyd Warshall:",(shortest_path(X, X, L, L, "floyd_warshall")))
    print("Auto:",shortest_path(X, X, L, L, "auto"))

def test_subtree_RG():
    print("Subtree [RG]:",subtree_RG(X,X,L,L,2))

def test_graphlets_sampling():
    print("Graphlets Sampling:",graphlets_sampling(X, X, 5, 0.05, 0.05,-1))

def test_weisfeiler_lehman():
    base_kernel = dict()
    base_kernel["dirac"] = lambda x, y: dirac_inner(x,y)
    base_kernel["shortest path"] = lambda x, y: shortest_path_inner(x,y)
    base_kernel["subtree"] = lambda x, y: subtree_RG_inner(x,y)
    for k in base_kernel.keys():
        print("Weisfeiler_lehman - "+str(k)+":",weisfeiler_lehman(X,X,L,L,base_kernel[k],5))

test_dirac()
test_random_walk_simple()
test_random_walk_sylvester()
test_shortest_path()
test_subtree_RG()
test_graphlets_sampling()
test_weisfeiler_lehman()
