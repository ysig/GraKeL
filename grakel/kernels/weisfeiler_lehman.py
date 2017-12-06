""" This file contains the shortest graphlte
    as defined by [Shervashize et al., 2011]
"""
import itertools

import numpy as np

from ..graph import graph

def weisfeiler_lehman(X, Y, Lx, Ly, base_kernel, niter=5):
    """ Computes the Weisfeler Lehman as proposed
        at 2011 by shervashize et al.

        X,Y: Valid graph formats to be compared
        Lx, Ly: Valid labels for graphs
        base_kernel: A valid base kernel
    """
    Ga = graph(X,Lx)
    Gb = graph(Y,Ly)
    return weisfeiler_lehman_inner(Ga, Gb, base_kernel, niter)

def weisfeiler_lehman_inner(Ga, Gb, base_kernel, niter=5):
    """ Computes the Weisfeler Lehman as proposed
        at 2011 by shervashize et al.

        Ga,Gb: Graph type formats of the graphs
               to be compared
        Lx, Ly: Valid labels for graphs

        base_kernel: A valid kernel function of the form:
        graph_A, graph_B -> number
    """
    if niter<1:
        raise ValueError('n_iter must be an integer bigger than zero')
    
    Ga.desired_format("dictionary")
    Gb.desired_format("dictionary")
    
    La_original = Ga.get_labels(purpose="dictionary")
    Lb_original = Gb.get_labels(purpose="dictionary")
    
    Ga_edge_dictionary = Ga.edge_dictionary
    Gb_edge_dictionary = Gb.edge_dictionary

    WL_labels = dict()
    WL_labels_inverse = dict()
    # get all the distinct values of current labels
    distinct_values = sorted(list(set(La_original.values()).union(set(Lb_original.values()))))
    # assign a number to each label
    label_count = 0
    for dv in distinct_values:
        WL_labels[label_count] = dv
        WL_labels_inverse[dv] = label_count
        label_count +=1

    La = dict()
    Lb = dict()
    for k in La_original.keys():
        La[k] = WL_labels_inverse[La_original[k]]
    for k in Lb_original.keys():
        Lb[k] = WL_labels_inverse[Lb_original[k]]

    # add new labels
    Ga.relabel(La)
    Gb.relabel(Lb)

    kernel = base_kernel(Ga, Gb)
    for i in range(niter):
        label_set = set()
       
        # Find unique labels and sort
        # them for both graphs
        # Keep for each node the temporary
        La_temp = dict()
        for v in Ga_edge_dictionary.keys():
            nlist = list()
            for neighbor in Ga_edge_dictionary[v].keys():
                nlist.append(La[neighbor])
            credential = str(La[v])+","+str(sorted(nlist))
            La_temp[v] = credential
            label_set.add(credential)
            
        # Dictionary for second graph
        Lb_temp = dict()
        for v in Gb_edge_dictionary.keys():
            nlist = list()
            for neighbor in Gb_edge_dictionary[v].keys():
                nlist.append(Lb[neighbor])
            credential = str(Lb[v])+","+str(sorted(nlist))
            Lb_temp[v] = credential
            label_set.add(credential)

        label_list = sorted(list(label_set))
        for dv in label_list:
            WL_labels[label_count] = dv
            WL_labels_inverse[dv] = label_count
            label_count +=1

        # Recalculate labels
        La = dict()
        Lb = dict()
        for k in La_temp.keys():
            La[k] = WL_labels_inverse[La_temp[k]]
        for k in Lb_temp.keys():
            Lb[k] = WL_labels_inverse[Lb_temp[k]]        
            
        # relabel
        Ga.relabel(La)
        Gb.relabel(Lb)

        # calculate kernel 
        kernel += base_kernel(Ga, Gb)
        
    # Restore original labels  
    Ga.relabel(La_original)
    Gb.relabel(Lb_original)
    return kernel

def weisfeiler_lehman_matrix(Graphs_a, base_kernel_matrix, Graphs_b=None,  niter=5):
    """ Computes the Weisfeler Lehman as proposed
        at 2011 by shervashize et al.

        Graphs_{a,b}: dictionary of graph type objects with keys
                      from 0 to the number of values
        Lx, Ly: Valid labels for graphs

        base_kernel_matrix: A valid kernel function of the form:
            graph_dict_A, graph_dict_B -> n by m np.array
            where n = len(graph_dict_A.keys())
              and m = len(graph_dict_B.keys())
    """
    if niter<1:
        raise ValueError('n_iter must be an integer bigger than zero')
    
    target_is_self = (Graphs_b is None)

    ng_a = len(Graphs_a.keys())
      
    if Graphs_b==None:
        g_iter = [Graphs_a[i] for i in range(0, ng_a)]
        ng = ng_a
        ng_b = ng_a
    else:
        ng_b = len(Graphs_b.keys())
        g_iter = itertools.chain([Graphs_a[i] for i in range(0, ng_a)], [Graphs_b[i] for i in range(0, ng_b)])
        ng = ng_a + ng_b
    
    Gs = dict()    
    G_ed = dict()
    L_orig = dict()
    
    # get all the distinct values of current labels
    WL_labels = dict()
    WL_labels_inverse = dict()
    distinct_values = set()
    
    for (i, g) in enumerate(g_iter):
        g.desired_format("dictionary")
        Gs[i] = g
        G_ed[i] = g.edge_dictionary
        L_orig[i] = g.node_labels
        # calculate all the distinct values
        distinct_values |= set(L_orig[i].values())
    distinct_values = sorted(list(distinct_values))
    
    
    # assign a number to each label
    label_count = 0
    for dv in distinct_values:
        WL_labels[label_count] = dv
        WL_labels_inverse[dv] = label_count
        label_count +=1
        
    for i in range(ng):
        L = dict()
        for k in L_orig[i].keys():
            L[k] = WL_labels_inverse[L_orig[i][k]]

        # add new labels
        Gs[i].relabel(L)
    
    kernel = np.zeros(shape=(ng_b, ng_a))
    kernel = np.add(kernel,base_kernel_matrix(Graphs_a, Graphs_b))
        
    for i in range(niter):
        label_set = set()
        L_temp=dict()
        for j in range(ng):
            # Find unique labels and sort
            # them for both graphs
            # Keep for each node the temporary
            L_temp[j] = dict()
            for v in G_ed[i].keys():
                nlist = list()
                for neighbor in G_ed[j][v].keys():
                    nlist.append(Gs[j].node_labels[neighbor])
                credential = str(Gs[j].node_labels[v])+","+str(sorted(nlist))
                L_temp[j][v] = credential
                label_set.add(credential)
            
        label_list = sorted(list(label_set))
        for dv in label_list:
            WL_labels[label_count] = dv
            WL_labels_inverse[dv] = label_count
            label_count +=1

        # Recalculate labels
        for j in range(ng):
            L = dict()
            for k in L_temp[j].keys():
                L[k] = WL_labels_inverse[L_temp[j][k]]
            # relabel
            Gs[j].relabel(L)
            
        # calculate kernel 
        kernel = np.add(kernel, base_kernel_matrix(Graphs_a, Graphs_b))
        
    # Restore original labels  
    for i in range(ng):
        Gs[i].relabel(L_orig[i])
    
    if Graphs_b is None:
        kernel = np.triu(kernel) + np.triu(kernel, 1).T
    
    return kernel

