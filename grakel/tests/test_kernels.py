import numpy as np
from numpy.testing import assert_almost_equal

from grakel.kernels import *

global X, L, Le, phi

X = np.array([[0,1,2,1,0,0.5,0],
              [1,0,0,0,1,2,0.5],
              [2,0,0,3,0,0,2],
              [1,0,3,0,0,0,0],
              [0,1,0,0,0,3,1],
              [0.5,2,0,0,3,0,0],
              [0,0.5,2,0,1,0,0]])
              
L = {0:'banana', 1:'cherry', 2:'banana', 3:'cherry', 4:'peach',5:'cherry',6:'lime'}

phi = {0:[0,1,2], 1:[0,0,1], 2:[1,1,0], 3:[3,0,1], 4:[0,4,0], 5:[1,1,1], 6:[0,1,0]}

Le = {(0,0):'b',(0,1):'c',(0,2):'d',(0,3):'b',(0,4):'c',(0,5):'a',(0,6):'d',
(1,0):'b',(1,1):'c',(1,2):'d',(1,3):'b',(1,4):'c',(1,5):'a',(1,6):'d',
(2,0):'d',(2,1):'c',(2,2):'y',(2,3):'b',(2,4):'r',(2,5):'q',(2,6):'a',
(3,0):'b',(3,1):'a',(3,2):'d',(3,3):'w',(3,4):'c',(3,5):'q',(3,6):'w',
(4,0):'a',(4,1):'c',(4,2):'u',(4,3):'q',(4,4):'t',(4,5):'a',(4,6):'r',
(5,0):'a',(5,1):'y',(5,2):'d',(5,3):'a',(5,4):'t',(5,5):'c',(5,6):'t',
(6,0):'w',(6,1):'c',(6,2):'d',(6,3):'w',(6,4):'r',(6,5):'c',(6,6):'w'}

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
    print("Subtree [RG]:",subtree_rg(X,X,L,L,2))

def test_graphlet_sampling():
    print("Graphlets Sampling:",graphlet_sampling(X, X, 5, 0.05, 0.05,-1))

def test_weisfeiler_lehman():
    base_kernel = dict()
    base_kernel["dirac"] = lambda x, y: dirac_inner(x,y)
    base_kernel["shortest path"] = lambda x, y: shortest_path_inner(x,y)
    base_kernel["subtree"] = lambda x, y: subtree_rg_inner(x,y)
    for k in base_kernel.keys():
        print("Weisfeiler_lehman - "+str(k)+":",weisfeiler_lehman(X,X,L,L,base_kernel[k],5))
        
def test_multiscale_laplacian():
    print("Multiscale Laplacian:",multiscale_laplacian(X,X,phi,phi))

def test_subgraph_matching():
    print("Subgraph Matching:",subgraph_matching(X,X,L,L,Le,Le))

"""
test_dirac()
test_random_walk_simple()
test_random_walk_sylvester()
test_shortest_path()
test_subtree_RG()
test_graphlet_sampling()
test_weisfeiler_lehman()
test_multiscale_laplacian()
"""
test_subgraph_matching()
