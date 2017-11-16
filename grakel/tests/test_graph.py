import numpy as np
from numpy.testing import assert_almost_equal

from grakel.graph import graph

def test_graph_adjacency():
    X = np.array([[1,1,0,3],[1,0,0,2],[2,3,0,1],[1,0,0,0]])
    labels = {0:'banana', 1:'cherry', 2:'banana', 3:'cherry'}
    # try all formats
    g = dict()
    g["auto"] = graph(X,labels,"auto")
    g["dict"] = graph(X,labels,"dictionary")
    g["adjc"] = graph(X,labels,"adjacency")
    g["all"] = graph(X,labels,"all")
    for k in g.keys():
        print(k)
        print(g[k].get_label_group(),'\n')

    for k in g.keys():
        print(k)
        sp1, sp2 = g[k].build_shortest_path_matrix(algorithm_type="auto")
        print(sp1,'\n',sp2,'\n')

    #assert_almost_equal(estimator.predict(X), X[:, 0]**2)
def test_graph_edge_dictionary():
    X = {'a':{'a':1,'b':1,'d':3},'b':{'a':1,'d':2},'c':{'a':2,'b':3,'d':1},'d':{'a':1}}
    labels = {'a':'banana', 'b':'cherry', 'c':'banana', 'd':'cherry'}
    g = dict()
    g["auto"] = graph(X,labels,"auto")
    g["dict"] = graph(X,labels,"dictionary")
    g["adjc"] = graph(X,labels,"adjacency")
    g["all"] = graph(X,labels,"all")
    for k in g.keys():
        print(k)
        print(g[k].get_label_group(),'\n')

    for k in g.keys():
        print(k)
        sp1, sp2 = g[k].build_shortest_path_matrix(algorithm_type="auto")
        print(sp1,'\n',sp2,'\n')

