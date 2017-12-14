import itertools

import numpy as np
import numpy.testing as npt

from grakel.graph_kernels import GraphKernel

global verbose, main, development

main = True
verbose = False
development = False


global X, L, k

X = np.array([[0,1,2,1,0,0.5,0],
              [1,0,0,0,1,2,0.5],
              [2,0,0,3,0,0,2],
              [1,0,3,0,0,0,0],
              [0,1,0,0,0,3,1],
              [0.5,2,0,0,3,0,0],
              [0,0.5,2,0,1,0,0]])

# Route May 2001 map of Spirit Airlines
# Atlantic City, Chicago (O'Hare), Detroit, Fort Lauderdale, Fort Myers, Los Angeles, Melbourne, Myrtle Beach, Newark, New York (LaGuardia), Oakland, Orlando, Tampa, and Washington (Reagan National).
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
              
L = {0:'banana', 1:'cherry', 2:'banana', 3:'cherry', 4:'peach',5:'cherry',6:'lime'}

phi = {0:[0,1,2], 1:[0,0,1], 2:[1,1,0], 3:[3,0,1], 4:[0,4,0], 5:[1,1,1], 6:[0,1,0]}

Le = {(0,0):'b',(0,1):'c',(0,2):'d',(0,3):'b',(0,4):'c',(0,5):'a',(0,6):'d',
(1,0):'b',(1,1):'c',(1,2):'d',(1,3):'b',(1,4):'c',(1,5):'a',(1,6):'d',
(2,0):'d',(2,1):'c',(2,2):'y',(2,3):'b',(2,4):'r',(2,5):'q',(2,6):'a',
(3,0):'b',(3,1):'a',(3,2):'d',(3,3):'w',(3,4):'c',(3,5):'q',(3,6):'w',
(4,0):'a',(4,1):'c',(4,2):'u',(4,3):'q',(4,4):'t',(4,5):'a',(4,6):'r',
(5,0):'a',(5,1):'y',(5,2):'d',(5,3):'a',(5,4):'t',(5,5):'c',(5,6):'t',
(6,0):'w',(6,1):'c',(6,2):'d',(6,3):'w',(6,4):'r',(6,5):'c',(6,6):'w'}

k = 5

def gk_test_dirac():
    XX = list(zip(k*[X],k*[L]))
    gk = GraphKernel(kernel={"name":"dirac"})
    gkf = gk.fit(XX)
    if verbose:
        print("Dirac:", gkf.transform())
    else:
        XX_correct = np.full((k, k), 15)
        npt.assert_array_equal(XX_correct, gkf.transform())

def gk_test_random_walk_simple():
    XX = list(zip(k*[X],k*[L]))
    gk = GraphKernel(kernel={"name":"random_walk", "lamda":0.1, "method_type":"simple"})
    gkf = gk.fit(XX)
    if verbose:
        print("Simple:",gkf.transform())
    else:
        XX_correct = np.full((k, k), -30.912616526802676)
        npt.assert_array_equal(XX_correct, gkf.transform())
    
def gk_test_random_walk_sylvester():
    XX = list(zip(k*[X],k*[L]))
    gk = GraphKernel(kernel={"name":"random_walk", "lamda":0.1, "method_type":"sylvester"})
    gkf = gk.fit(XX)
    if verbose:
        print("Sylvester:",gkf.transform())
    else:
        XX_correct = np.full((k, k), -30.912616526802676)
        npt.assert_array_almost_equal(XX_correct, gkf.transform())
def gk_test_shortest_path():
    XX = list(zip(k*[X],k*[L]))
    gk = GraphKernel(kernel={"name":"shortest_path", "algorithm_type":"dijkstra"})
    gkf = gk.fit(XX)
    if verbose:
        print("Dijkstra:", gkf.transform())
    else:
        XX_correct = np.full((k, k), 66)
        npt.assert_array_equal(XX_correct, gkf.transform())
        
    gk = GraphKernel(kernel={"name":"shortest_path", "algorithm_type":"floyd_warshall"})
    gkf = gk.fit(XX)
    if verbose:
        print("Floyd Warshall:", gkf.transform())
    else:
        XX_correct = np.full((k, k), 66)
        npt.assert_array_equal(XX_correct, gkf.transform())

    
    gk = GraphKernel(kernel={"name":"shortest_path", "algorithm_type":"auto"})
    gkf = gk.fit(XX)
    if verbose:
        print("Auto:", gkf.transform())
    else:
        XX_correct = np.full((k, k), 66)
        npt.assert_array_equal(XX_correct, gkf.transform())

def gk_test_subtree_rg():
    XX = list(zip(k*[X],k*[L]))
    gk = GraphKernel(kernel={"name":"subtree_rg", "h":5})
    gkf = gk.fit(XX)
    if verbose:
        print("Subtree [RG]:", gkf.transform())
    else:
        bign = 15458150092069060998865743245478022133789841061117106540018232182730681287234085528132318819741260764317676931326456376488696020676136950620929717075828209329592328934916049        
        XX_correct = np.full((k, k), bign)
        npt.assert_array_equal(XX_correct, gkf.transform())

def gk_test_graphlets_sampling():
    XX = list(zip(k*[X],k*[L]))
    gk = GraphKernel(kernel={"name":"graphlets_sampling", "k": 5, "delta" : 0.05, "epsilon" : 0.05, "a": -1})
    gkf = gk.fit(XX)
    if verbose:
        print("Graphlets Sampling:", gkf.transform())
    else:
        XX_correct = np.full((k, k), 0.49863760218)
        npt.assert_array_almost_equal(XX_correct, gkf.transform())
        
def gk_test_weisfeiler_lehman():
    XX = list(zip(k*[X],k*[L]))
    base_kernel = dict()
    base_kernel["dirac"] = {"name":"dirac"}
    base_kernel["shortest path"] = {"name":"shortest_path"}
    base_kernel["subtree"] = {"name":"subtree_rg"}
    
    npt_v = dict()
    npt_v["dirac"] = 50
    npt_v["shortest path"] = 276
    npt_v["subtree"] = 15458150092069060998865743245478022133789841061117106540018232182730681287234085528132318819741260764317676931326456376488696020676136950620929717075828209329606781669103994

    for key in base_kernel.keys():
        gk = GraphKernel(kernel=[{"name":"weisfeiler_lehman","niter":5},base_kernel[key]])
        gkf = gk.fit(XX)
        if verbose:
            print("Weisfeiler_lehman - "+str(key)+":",gkf.transform())
        else:
            XX_correct = np.full((k, k), npt_v[key])
            npt.assert_array_equal(XX_correct, gkf.transform())
            
def gk_test_multiscale_laplacian():
    XX = k*[[X,phi]]
    gk = GraphKernel(kernel={"name":"multiscale_laplacian"})
    gkf = gk.fit(XX)
    
    if verbose:
        print("Multiscale Laplacian:", gkf.transform())

def gk_test_subgraph_matching():
    XX = list(k*[[X,L,Le]])
    gk = GraphKernel(kernel={"name":"subgraph_matching"})
    gkf = gk.fit(XX)
    
    if verbose:
        print("Subgraph Matching:", gkf.transform())

def gk_test_lovasz_theta():
    XX = list(k*[[D]])
    gk = GraphKernel(kernel={"name":"lovasz_theta"})
    gkf = gk.fit(XX)

    if verbose:
        print("Lovasz Theta:", gkf.transform())

def gk_test_svm_theta():
    XX = list(k*[[D]])    
    gk = GraphKernel(kernel={"name":"svm_theta"})
    gkf = gk.fit(XX)

    if verbose:
        print("SVM Theta:", gkf.transform())

def gk_test_neighborhood_pairwise_subgraph_distance():
    XX = list(k*[[X,L,Le]])
    gk = GraphKernel(kernel={"name":"neighborhood_pairwise_subgraph_distance"})
    gkf = gk.fit(XX)

    if verbose:
        print("NSPDK:", gkf.transform())

def gk_test_neighborhood_hash():
    XX = list(zip(k*[X],k*[L]))
    gk = GraphKernel(kernel={"name":"neighborhood_hash", "nh_type":"simple"})
    
    gkf = gk.fit(XX)
    if verbose:
        print("Neighborhood Hash - 'simple':", gkf.transform())
    
    gk = GraphKernel(kernel={"name":"neighborhood_hash", "nh_type":"count-sensitive"})
    gkf = gk.fit(XX)
    
    if verbose:
        print("Neighborhood Hash - 'count-sensitive':", gkf.transform())
   
def gk_test_odd_sth():
    l = {0:'A', 1:'B', 2:'C', 3:'D', 4:'E',5:'F',6:'G'}

    XX = list(k*[[X,l]])
    gk = GraphKernel(kernel={"name":"odd_sth"})
    gkf = gk.fit(XX)

    if verbose:
        print("ODD-STh:", gkf.transform())
 
def gk_test_wl_nh():
    XX = list(zip(k*[X],k*[L]))
    base_kernel = {"name":"neighborhood_hash", "nh_type":"count-sensitive"}
    
    gk = GraphKernel(kernel=[{"name":"weisfeiler_lehman","niter":5},base_kernel])
    gkf = gk.fit(XX)
    
    if verbose:
        print("Weisfeiler Lehman - Neighboorhood Hash:",gkf.transform())
        
def gk_test_propagation():
    XX = list(zip(k*[X],k*[L]))
    gk = GraphKernel(kernel={"name":"propagation"})
    gkf = gk.fit(XX)
    
    if verbose:
        print("Propagation:", gkf.transform())

def gk_test_pyramid_match():
    XX = list(zip(k*[X],k*[L]))
    gk = GraphKernel(kernel={"name":"pyramid_match"})
    gkf = gk.fit(XX)
    
    if verbose:
        print("pyramid match:", gkf.transform())

if verbose and main:
    gk_test_dirac()
    gk_test_random_walk_simple()
    gk_test_random_walk_sylvester()
    gk_test_shortest_path()
    gk_test_subtree_rg()
    gk_test_graphlets_sampling()
    gk_test_weisfeiler_lehman()
    gk_test_subgraph_matching()
    gk_test_lovasz_theta()
    gk_test_svm_theta()
    gk_test_neighborhood_pairwise_subgraph_distance()    
    gk_test_neighborhood_hash()
    gk_test_wl_nh()
    gk_test_odd_sth()
    gk_test_propagation()
    
if verbose and development:
    gk_test_pyramid_match()
#    gk_test_multiscale_laplacian()
