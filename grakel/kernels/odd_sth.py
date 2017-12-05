""" This file contains the ODD-Sth kernel
    as defined in [San Martino et al. 2009]
"""

import itertools

import numpy as np

from ..graph import graph

def odd_sth(X, Y, Lx, Ly, h=None):
    """ The ODD-Sth kernel as proposed
        in [San Martino et al. 2009]

        X,Y: Valid graph formats to be compared
    """
    Gx = graph(X,Lx)
    Gy = graph(Y,Ly)
    return float(odd_sth_matrix({0: Gx}, {0: Gy}, h=None)[0,0])

def odd_sth_matrix(Graphs_x, Graphs_y=None, h=None):
    """ The ODD-Sth kernel as proposed
        in [San Martino et al. 2009]

        Gx, Gy: Graph type objects
        h: depth of visit
    """
    h_x = len(Graphs_x.keys())
    
    if Graphs_y==None:
        h_y = h_x
        g_iter = Graphs_x.values()
        ng = h_x
        pairs = [(i,j) for i in range(0,h_x) for j in range(i, h_x)]
        offset = 0
    else:
        h_y = len(Graphs_y.keys())
        g_iter = itertools.chain(Graphs_x.values(), Graphs_y.values())
        ng = h_x+h_y
        pairs = list(itertools.product(range(h_x, ng), range(0,h_x)))
        offset = h_x
        
    big2Dag = None
    for g in g_iter:
        big2Dag = big_dag_append(make_big_dag(g, h=h), big2Dag, merge_features=False)
        
    C = dict()
    for v in big2Dag[0].keys():
        # number of identical subtrees
        # equal the D element
        C[v] = big2Dag[0][v][0]
    
    K = np.array(shape=(h_y, h_x))
    for (i,j) in pairs:
        for v in C.keys():
            #TODO write K calculation in terms of matrix multiplication
            K[i-offset,j] = big2Dag[0][v][1][i]*big2Dag[0][v][1][j]*C[v]
    
    if Graphs_y is None:
        K = np.tril(K) + np.tril(K, -1).T
    
    return K


def make_big_dag(g, h=None):
    """
        A function that creates a list
        of dags (edge_dictionaries)
        for the given graph, and the
        labels given whom a sorting has been made
    """
    big_dag = None
    for v in g.get_vertices(purpose='any'):
        dag_odd = bfs_odd(v, g, h=h)
        dag = tuple(hash_trees(bfs_odd_out))+([])+(dag_odd[1],dag_odd[3])
        big_dag = big_dag_append(dag, big_dag)
    
    _, _, D_ordering, _ = odd((big_dag[0], big_dag[2], big_dag[3]))
    big_dag = (big_dag[0], big_dag[1], sorted(list(big_dag[0].keys()), key = lambda x: D_ordering[x]), big_dag[2], big_dag[3])
        
    return big_dag
    
def bfs_odd(v,g, h=None):
    """
        create dags and apply topological sorting
    """
    vertices, edges = bfs(v, g, h=None)
    return odd(vertices, edges, g.get_labels(purpose='any'))
    
def bfs(v, g, h=None):
    """
        BFS exploration that returns a dag
    """
    if h is None:
        h = -1
    else:
        h = int(h)
        if h<=0:
            raise ValueError('depth of visits must be bigger than zero')
        
    q = [v]
    vertices = set()
    edges = dict()
    level = 0
    while len(q)>0:
        v = q.pop(0)
        if v not in edges:
            edges[v] = list()
        for n in g.neighbors(v, purpose='any'):
            edges[v].append(n)
            if n not in vertices:
                vertices.add(n)
        level+=1
        if level == h:
            break
            
    return vertices, edges

def odd(vertices, edges, labels):
    """
        topological sorting of a dag
    """
    # Kahn's algorithm for topological
    # sorting of dags
    # order graph:
    
    # Step-1: Compute in-degree (number of incoming edges) 
    # for each of the vertex present in the DAG and initialize
    # the count of visited nodes as 0.
    indegrees = dict()
    
    if type(vertices) is set:
        zero_indegrees = vertices
        visited_nodes = len(vertices)
    elif type(vertices) is dict:
        zero_indegrees = set(vertices.keys())
        visited_nodes = len(vertices.keys())
    else:
        raise ValueError('unsupported vertices type')
        
    for (k, e) in edges.items():
        for v in e:
            if v not in indegrees:
                indegrees[v] = 1
            else:
                indegreed[v] +=1
            zero_indegrees.discard(v)
    
    # Step-2: Pick all the vertices with in-degree as 
    # 0 and add them into a queue

    q = list(zero_indegrees)
    
    ordering = dict()
    while len(q)>0:
        # Step-3: Remove a vertex from the queue and then
        # Increment count of visited nodes by 1.
        # Decrease in-degree by 1 for all its neighboring nodes.
        # If in-degree of a neighboring nodes is reduced to zero,
        # then add it to the queue.
    
        q = sorted(q, key = lambda x: labels[x])
        e = q.pop(0)
        ordering[e] = visited_nodes
        for k in edges[e]:
            if k in indegrees[k]:
                if indegrees[k]==1:
                    indegrees.pop(k)
                    q.append(k)
        visited_nodes -= 1
    
    # apply the ordering                
    for k in edges.keys():
        edges[k] = sorted(edges[k], key = lambda x: ordering[x])
    
    return vertices, edges, ordering, labels

def hash_trees(tree):
    (vertex, edge, ordering, labels) = tree
    
    v_ordered = sorted(list(vertex), key = lambda x: ordering[x])
    vertices_hash_map = dict()
    vertices = dict()
    for v in v_ordered:            
        if v not in edge or len(edge[v])==0:
            if label[v] not in vertices_hash_map:
                vertices_hash_map[labels[v]] = list()
            vertices_hash_map[labels[v]] += v
            vertices[v] = (0,1,str(labels[v]))
        else:
            neighbors_ids = []
            d = 0
            for n in edge[v]:
                d += 1 + vertices[n][0]
                neighbors_ids.append(vertices[v][2])
            vertices[v]=(d, 1, str(labels[v])+'('+(','.join(n))+')')
            
    return (vertices, vertices_hash_map, v_ordered)
            
def big_dag_append(dag, big_dag = None, merge_features=True):
    """ Calculates the minimal DAG as 
        defined in [Aiolli, San Martino et al., (2006)]
        notated in [San Martino et al., (2009) as BigDAG.
        
        dags: a list of dags assumed to be in inverse
              topological order
    """
    if big_dag is None:
        D_labels = dict()
        D_hash_map = dict()
        D_vertex = dict()
        D_edges = dict()
        nodes_idx = 0
        nf = 1
    else:
        (D_vertex, D_hash_map, D_edges, _, D_labels) = big_dag
        nodes_idx = len(D_vertex.keys())
        if not merge_features:
            f = True
            for v in D_vertex.keys():
                D_vertex[v][1].append(0)
                if f:
                    nf = len(D_vertex[v][1])
                    f = False
            if f:
                nf = 1
                
    (vertices, vertices_hash_map, v_ordered, edges, labels) = dag
    for q in v_ordered:
        key = vertices[q][2]
        if key in D_hash_map:
            node = D_hash_map[key][0]
            # update frequency
            if merge_features:
                D_vertices[node][1] += vertices[q][1]
            else:
                D_vertices[node][1][-1] += vertices[q][1]
        else:
            D_labels[nodes_idx] = labels[q]
            D_edges[nodes_idx] = list()
            n = list()
            for c in edges[q]:
                ck = vertices[c][2]
                if ck in D_hash_map:
                    node = D_hash_map[key][0]
                    D_edges[nodes_idx].append(node)
                    n += D_vertex[nodes_idx][2]
                    d += D_vertex[nodes_idx][0]
            
            new_key = str(D_labels[nodes_idx])+'('+(','.join(n))+')'
            if new_key not in D_hash_map:
                D_hash_map[new_key] = list()
            D_hash_map[new_key].append(node_idx) 
            
            if merge_features:
                freq = 1
            else:
                freq = (nf-1)*[0]+[1]
                
            D_vertex[nodes_idx] = (d, freq, new_key)
            nodes_idx += 1
                
    return (D_vertex, D_hash_map, D_edges, D_labels)
    
