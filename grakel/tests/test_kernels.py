import numpy as np
import numpy.testing as npt

from grakel.kernels import *

global verbose, main, development

main = True
verbose = False 
development = True

global X, Y, L, Le, phi

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

LD = {0:'rock', 1:'post-punk', 2:'rock', 3:'indie', 4:'classical', 5:'pop', 6:'rock', 7:'punk', 8:'punk', 9:'indie', 10:'post-rock', 11:'post-punk', 12:'jazz', 13:'jazz', 14:'jazz', 15:'classical'}
              
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
    if verbose:
        print("Dirac:",dirac(X, X, L, L))
    else:
        npt.assert_equal(15, dirac(X, X, L, L))
def test_random_walk_simple():
    if verbose:
        print("Random Walk [Simple]:",random_walk(X, X, lamda=0.1, method_type="simple"))
    else:
        npt.assert_almost_equal(-30.912616526802676, random_walk(X, X, lamda=0.1, method_type="simple"), decimal=3)
def test_random_walk_sylvester():
    if verbose:
        print("Random Walk [Sylvester]:",random_walk(X, X, lamda=0.1, method_type="sylvester"))
    else:
        npt.assert_almost_equal(-30.912616526802676, random_walk(X, X, lamda=0.1, method_type="sylvester"), decimal=3)
def test_shortest_path():
    if verbose:
        print("Shortest Path [Dijkstra]:",shortest_path(X, X, L, L, "dijkstra"))
        print("Shortest Path [Floyd Warshall]:",(shortest_path(X, X, L, L, "floyd_warshall")))
        print("Shortest Path [Auto]:",shortest_path(X, X, L, L, "auto"))
    else:
        npt.assert_equal(66, shortest_path(X, X, L, L, "dijkstra")) 
        npt.assert_equal(66, shortest_path(X, X, L, L, "floyd_warshall")) 
        npt.assert_equal(66, shortest_path(X, X, L, L, "auto")) 

def test_subtree_RG():
    if verbose:
        print("Subtree [RG]:",subtree_rg(X,X,L,L,2))
    else:
        bign = 15458150092069060998865743245478022133789841061117106540018232182730681287234085528132318819741260764317676931326456376488696020676136950620929717075828209329592328934916049
        npt.assert_equal(bign, subtree_rg(X,X,L,L,2)) 

def test_graphlet_sampling():
    if verbose:
        print("Graphlets Sampling:",graphlet_sampling(X, X, 5, 0.05, 0.05,-1))
    else:
         npt.assert_almost_equal(0.4986,graphlet_sampling(X, X, 5, 0.05, 0.05,-1), decimal=3)
def test_weisfeiler_lehman():
    base_kernel = dict()
    base_kernel["dirac"] = lambda x, y: dirac_inner(x,y)
    base_kernel["shortest path"] = lambda x, y: shortest_path_inner(x,y)
    base_kernel["subtree"] = lambda x, y: subtree_rg_inner(x,y)
    
    if not verbose:
        npt_v = dict()
        npt_v["dirac"] = 50
        npt_v["shortest path"] = 276
        npt_v["subtree"] = 15458150092069060998865743245478022133789841061117106540018232182730681287234085528132318819741260764317676931326456376488696020676136950620929717075828209329606781669103994
    
    for k in base_kernel.keys():
        if verbose:
            print("Weisfeiler_lehman - "+str(k)+":",weisfeiler_lehman(X,X,L,L,base_kernel[k],5))
        else:
            npt.assert_equal(npt_v[k],weisfeiler_lehman(X,X,L,L,base_kernel[k],5))
            
def test_multiscale_laplacian():
    if verbose:
        print("Multiscale Laplacian:",multiscale_laplacian(X,X,phi,phi))

def test_subgraph_matching():
    if verbose:
        print("Subgraph Matching:",subgraph_matching(X,X,L,L,Le,Le))

def test_lovasz_theta():
    if verbose:
        print("Lovasz Theta:", lovasz_theta(D,D))

def test_svm_theta():
    if verbose:
        print("SVM Theta:", svm_theta(D,D))

def test_neighborhood_pairwise_subgraph_distance_kernel():
    if verbose:
        print("NSPDK:",neighborhood_pairwise_subgraph_distance(X, X, L, L, Le, Le))

def test_neighborhood_hash_kernel():
    if verbose:
        print("Neighborhood Hash - 'simple':",neighborhood_hash(X, X, L, L, nh_type='simple'))
        print("Neighborhood Hash - 'count-sensitive':",neighborhood_hash(X, X, L, L, nh_type='count-sensitive'))
        
def test_odd_sth():
    #x = np.array([[0, 1, 1, 0], [1, 0, 1, 1], [1, 1, 0, 1], [0, 1, 1, 0]])
    #l = {0: 's', 1:'e', 2:'b', 3:'d'}
    if verbose:
        print("ODD-STh:",odd_sth(X, X, L, L, h=None))

def test_propagation():
    if verbose:
        print("Propagation:", propagation(X,X,L,L))

def test_pyramid_match():
    if verbose:
        print("Pyramid Match:", pyramid_match(D,D,LD,LD))

if verbose and main:
    test_dirac()
    test_random_walk_simple()
    test_random_walk_sylvester()
    test_shortest_path()
    test_subtree_RG()
    test_graphlet_sampling()
    test_weisfeiler_lehman()
    test_subgraph_matching() 
    test_lovasz_theta()
    test_svm_theta()    
    test_neighborhood_pairwise_subgraph_distance_kernel()
    test_neighborhood_hash_kernel()
    test_odd_sth()
    test_propagation()
if verbose and development:
#    test_multiscale_laplacian()
    test_pyramid_match()
