""" This file contains the ODD-Sth kernel
    as defined in [San Martino et al. 2009]
"""

import itertools

import numpy as np

from ..graph import graph

def odd_sth(X, Y):
    """ The ODD-Sth kernel as proposed
        in [San Martino et al. 2009]

        X,Y: Valid graph formats to be compared
    """
    Gx = graph(X)
    Gy = graph(Y)
    return svm_theta_inner(Gx, Gy)

def odd_sth_matrix(G_x, G_y=None, h=3):
    """ The ODD-Sth kernel as proposed
        in [San Martino et al. 2009]

        Gx, Gy: Graph type objects
    """
    # for all graphs
    # calculate ordered dags
    # calculate big dag (with frequencies)
    # from all big dags
    # calculate big dag (as hashes) 
    # of big dags
    # notated big2dag
    
    
    # 6.7
    # Sum u in V(Big2Dag)
    # Fu[i] * Fu[j] * C(u,u)
    # C: SubtreeKernel for trees
    # viswanathan
    # where u is subtree kernel
    return kernel

def create_sorted_dags(g):
    """
        A function that creates a list
        of dags (edge_dictionaries)
        for the given graph, and the
        labels given whom a sorting has been made
    """
    dags = list()
    for v in g.get_vertices():
        dags.append((bfs_odd(v,g))[0])
    return (dags, g.get_labesl())
    
def bfs_odd(v,g):
    """
        create dags and apply topological sorting
    """
    vertices, edges = bfs(v,g)
    vertices, edges, labels = odd(vertices, edges, g.get_labels())
    return vertices, edges, labels
    
def bfs(v,g):
    """
        BFS exploration that returns a dag
    """
    q = [v]
    vertices = set()
    edges = dict()
    level = 0
    while len(q)>0:
        v = q.pop(0)
        if v not in edges:
            edges[v] = list()
        for n in g.get_neighbors(v):
            edges[v].append(n)
            if n not in vertices:
                vertices.add(n)
        level+=1
        
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
    zero_indegrees = vertices.copy()
    for (i,k) in edges.items():
        if k not in outdegrees:
            indegrees[k] = 1
        else:
            indegreed[k] +=1
        zero_indegrees.dicard(k)
    
    # Step-2: Pick all the vertices with in-degree as 
    # 0 and add them into a queue

    q = list(zero_indegrees)
    visited_nodes = 0
    
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
    
    # apply the ordering                
    for k in edges.keys():
        edges[k] = sorted(edges[k], key = lambda x: ordering[x])
    
    return vertices, edges, labels
                
def calculate_big_dag(dags):
    
    
def st(u1,u2):
    """ The shortest-kernel for trees as proposed
        in [5.1, San Martino et al. 2009]
    """
    if (L(u1) != L(u2)):
        return 0
    else:
        # if u1, u2 preterminal nodes
            return lamda
        # else
        
