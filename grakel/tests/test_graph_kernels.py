import numpy as np
from numpy.testing import assert_almost_equal

from grakel.graph_kernels import GraphKernel

global X, L, k

X = np.array([[0,1,2,1,0,0.5,0],
              [1,0,0,0,1,2,0.5],
              [2,0,0,3,0,0,2],
              [1,0,3,0,0,0,0],
              [0,1,0,0,0,3,1],
              [0.5,2,0,0,3,0,0],
              [0,0.5,2,0,1,0,0]])

L = {0:'banana', 1:'cherry', 2:'banana', 3:'cherry', 4:'peach',5:'cherry',6:'lime'}

k = 5

def gk_test_dirac():
    XX = list(zip(k*[X],k*[L]))
    gk = GraphKernel(kernel={"name":"dirac"})
    gkf = gk.fit(XX)
    print("Dirac:",gkf.transform())

def gk_test_random_walk_simple():
    XX = list(zip(k*[X],k*[L]))
    gk = GraphKernel(kernel={"name":"random_walk", "lamda":0.1, "method_type":"simple"})
    gkf = gk.fit(XX)
    print("Simple:",gkf.transform())
    
def gk_test_random_walk_sylvester():
    XX = list(zip(k*[X],k*[L]))
    gk = GraphKernel(kernel={"name":"random_walk", "lamda":0.1, "method_type":"sylvester"})
    gkf = gk.fit(XX)
    print("Sylvester:",gkf.transform())

def gk_test_shortest_path():
    XX = list(zip(k*[X],k*[L]))
    gk = GraphKernel(kernel={"name":"shortest_path", "algorithm_type":"dijkstra"})
    gkf = gk.fit(XX)
    print("Dijkstra:", gkf.transform())

    gk = GraphKernel(kernel={"name":"shortest_path", "algorithm_type":"floyd_warshall"})
    gkf = gk.fit(XX)
    print("Floyd Warshall:", gkf.transform())
    
    
    gk = GraphKernel(kernel={"name":"shortest_path", "algorithm_type":"auto"})
    gkf = gk.fit(XX)
    print("Auto:", gkf.transform())

def gk_test_subtree_rg():
    XX = list(zip(k*[X],k*[L]))
    gk = GraphKernel(kernel={"name":"subtree_rg", "h":5})
    gkf = gk.fit(XX)
    print("Subtree [RG]:", gkf.transform())

def gk_test_graphlets_sampling():
    XX = list(zip(k*[X],k*[L]))
    gk = GraphKernel(kernel={"name":"graphlets_sampling", "k": 5, "delta" : 0.05, "epsilon" : 0.05, "a": -1})
    gkf = gk.fit(XX)
    print("Graphlets Sampling:", gkf.transform())

def gk_test_weisfeiler_lehman():
    XX = list(zip(k*[X],k*[L]))
    base_kernel = dict()
    base_kernel["dirac"] = {"name":"dirac"}
    base_kernel["shortest path"] = {"name":"shortest_path"}
    base_kernel["subtree"] = {"name":"subtree_rg"}
    for key in base_kernel.keys():
        gk = GraphKernel(kernel=[{"name":"weisfeiler_lehman","niter":5},base_kernel[key]])
        gkf = gk.fit(XX)
        print("Weisfeiler_lehman - "+str(key)+":",gkf.transform())

gk_test_dirac()
gk_test_random_walk_simple()
gk_test_random_walk_sylvester()
gk_test_shortest_path()
gk_test_subtree_rg()
gk_test_graphlets_sampling()
gk_test_weisfeiler_lehman()

