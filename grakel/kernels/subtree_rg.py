""" This file contains the subtree kernel
    as defined by [Shervashize et. al., 2011]
"""

import itertools

from ..graph import graph
from ..tools import nested_dict_get, nested_dict_add

def subtree_rg(X, Y, Lx, Ly, h=5):
    """ Computes the Ramon Gartner subtree kernel
        as proposed at 2011 by shervashize et. al.

        X, Y: Valid graph formats to be compared
        Lx, Ly: Valid labels for graphs
        base_kernel: A valid base kernel
    """
    Gx = graph(X,Lx)
    Gy = graph(Y,Ly)
    return subtree_rg_inner(Gx, Gy, h=5)

def subtree_rg_inner(Gx, Gy, h=5):
    """ Calculate The Ramon Gartner subtree kernel
        as proposed on [Ramon & Gartner, 2003]

        X_edge_dict, Y_edge_dict: The two graphs
            whose kernel is going to be computed
            as a dictionary of edges as follows:
                If (u,v) edge exists then edges["u"]["v"] has
                weight of the edge pointing from u to v

        Gx, Gy: Graph types objects

        height: The height of the subtree exploration         
    """
    Gx.desired_format("dictionary")
    Gy.desired_format("dictionary")
    kernel = 0
    dynamic_dictionary = dict()
    for u in Gx.edge_dictionary.keys():
        for v in Gy.edge_dictionary.keys():
            kernel += subtree_rg_core_dynamic(u,v,Gx,Gy,h,dynamic_dictionary)
    return kernel

def subtree_rg_core_dynamic(u, v, g_x, g_y, h, dynamic_dict, p_u=None, p_v=None):
    """ Calculate the inside of the summation
        of Ramon Gartner subtree kernel
        as proposed on [Ramon & Gartner, 2003]

        u,v: vertices that correspond to graphs g_x, g_y
        g_x, g_y: graph formats for the corresponding
        graphs 
        h: The height of the subtree exploration
        p_u, p_v: Needed for not to backtrack        
        
        This is an efficiency proposal using various
        data structures and dynamic programming
        assuming that for each node you take all 
        of its neighbors and not its preddecesor.
        A dictionary is being formed that holds
        level, u, v and pred of u and pred of v
    """
    
    if h==1:
        return int(g_x.label(u)==g_y.label(v))

    elif h>1:
        R = list()
        
        # Calculate R
        # First group all nodes with the same label
        # make a list of lists of all pairs for each label
        lbx = g_x.get_label_group()
        lby = g_y.get_label_group()
        Rset = []

        # Calculate neighbors and remove previous
        nx = g_x.neighbors(u)

        # Avoid going back
        if (p_u is not None):
            if p_u in nx:
                nx.remove(p_u)

        # do the same for y
        ny = g_y.neighbors(v)

        # Avoid going back
        if (p_v is not None):
            if p_v in ny:
                ny.remove(p_v)

        # What happens when one list has no neighbors?
        xy_and = len(nx)*len(ny)
        xy_or = len(nx)+len(ny)
        if(xy_or == 0):
            # both trees finish at the same point
            return int(g_x.label(u) == g_y.label(v))
        elif(xy_and == 0):
            # else trees are different so output zero (?)
            return 0
            
        # returns a list
        snx = set(nx)
        sny = set(ny)
        for kx in lbx:
            if kx in lby:
                # substract from the common label 
                # only the valid neighbors
                snxk = list(set(lbx[kx]).intersection(snx))
                snyk = list(set(lby[kx]).intersection(sny))
                if len(snxk) > 0 and len(snyk) > 0:
                    pair = [lbx[kx],lby[kx]]
                    Rset += list(itertools.product(*pair))

        # Designate all nodes with the same start
        # and all with the same finish
        right = dict()
        left = dict()

        # Dictionary to store for 
        # every pair the index
        Rset_dict = dict()
        Rset_inv_dict = dict()
        l = 0

        # assign values to indexes
        for (w,z) in Rset:
            if w not in right:
                right[w] = list()
            if z not in left:
                left[z] = list()
            right[w].append(l)
            left[z].append(l)
            Rset_dict[(w,z)] = l
            Rset_inv_dict[l] = (w,z)
            l+=1


        # Flatten both $right, $left to list of lists
        Kbins_flat_r = list()
        Kbins_flat_l = list()
        for k in right.keys():
            Kbins_flat_r.append(list(right[k]))
        for k in left.keys():
            Kbins_flat_l.append(list(left[k]))
            
        # calculate all possible combinations
        # product equals combinations
        # because the lists are disjoint
        setA = set(itertools.product(*Kbins_flat_r))
        setB = set(itertools.product(*Kbins_flat_l))

        # take the intersection of
        # found maximal valid R sets
        # all subsets of maximal sets are also valid
        Rmaximal = setA.intersection(setB)
        Valid_tuples = [Rset_inv_dict[i] for s in Rmaximal for i in s]
        
        # Calculate trees for the valid tuples
        # Rv assign values best on index
        Rv = dict()
        for (w,wp) in Valid_tuples:
            r = nested_dict_get(dynamic_dict, h-1, u, v, w, wp)
            if (r is None):
                kh = subtree_rg_core_dynamic(w, wp, g_x, g_y, h-1, dynamic_dict, u, v)
                nested_dict_add(dynamic_dict, kh, h-1, u, v, w, wp)
                Rv[Rset_dict[(w,wp)]] = kh
            else:
                Rv[Rset_dict[(w,wp)]] = r
        
        # Holds the values of all sets
        # inside R, that are not zero
        M_values = dict()
        for s in Rmaximal:
            # Keep only non zero elements
            # (all subsets containing them will be zero)
            non_zero_elements = set()
            for i in s:
                if(Rv[i] != 0):
                    non_zero_elements.add(i)

            if bool(non_zero_elements):
                # Can we calculate do calculation on subsets faster?
                # The answer is yes: dynamic programming
                plough_subsets(non_zero_elements, Rv, M_values)
        
        kernel = 0
        for v in M_values.values():
          kernel+=v
                
        return kernel
    else:
        # Raise Warning?
        return 0

def plough_subsets(initial_set, Rv, value):
    """ A function that finds all subset kernel
        values without repeating operations

        initial_set: The set you are operating (set)
        Rv: values for single elements (dictionary on integers)
        value: The dictionary of all sets and values (dictionary on frozensets)
    """
    # If only one element return its value
    if (len(initial_set) == 1):
        frozen_initial_set = frozenset(initial_set)
        p = initial_set.pop()
        value[frozen_initial_set] = Rv[p]
    
    elif (len(initial_set) > 1):
        # Explore all subsets for value dictionary to fill and calculate current
        flag = True
        frozen_initial_set = frozenset(initial_set)
        for s in initial_set:
            # copy set, delete element and if not calculated: plough
            temp_set = initial_set.copy()
            temp_set.discard(s)
            frozen_temp_set = frozenset(temp_set)
            if frozen_temp_set not in value:
                plough_subsets(temp_set, Rv, value)
            
            if(flag):
                value[frozen_initial_set] = Rv[s]*value[frozen_temp_set]
                flag = False
