""" This file contains the shortest graphlte
    as defined by [Shervashize et al., 2011]
"""

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
    Ga.desired_format("dictionary")
    Gb.desired_format("dictionary")
    
    La_original = Ga.get_labels("dictionary")
    Lb_original = Gb.get_labels("dictionary")
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
    k = base_kernel(Ga, Gb)
    for i in range(niter):
        label_set = set()
       
        # Find unique labels and sort
        # them for both graphs
        # Keep for each node the temporary
        La_temp = dict()
        for v in Ga_edge_dictionary.keys():
            nlist = list()
            for neighbour in Ga_edge_dictionary[v].keys():
                nlist.append(La[neighbour])
            credential = str(La[v])+","+str(sorted(nlist))
            La_temp[v] = credential
            label_set.add(credential)
            
        # Dictionary for second graph
        Lb_temp = dict()
        for v in Gb_edge_dictionary.keys():
            nlist = list()
            for neighbour in Gb_edge_dictionary[v].keys():
                nlist.append(Lb[neighbour])
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
        k += base_kernel(Ga, Gb)
        
        
    # Restore original labels  
    Ga.relabel(La_original)
    Gb.relabel(Lb_original) 
    return k

