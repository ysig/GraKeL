""" This file contains a set of functions implementing functions implementing
    various graph kernel metrics.

"""
import itertools

import numpy as np

from np.linalg import inv
from scipy.linalg import solve_sylvester

from graph import graph
from tools import inv_dict

def random_walk(X, Y,lamda=0.1, method_type="sylvester"):
    """ This is a function that implements the random walk kernel.

    X,Y: are two adjoint matrices representing relevant graphs to 
    be compared. Their format is a square array.

    """
    
    # calculate the product graph
    XY = np.kron(A,B)

    if(method_type == "simple"):
        # algorithm presented in [Kashima et al., 2003; Gärtner et al., 2003]
        # complexity of O(|V|^6)

        # XY is a square matrix
        s = XY.shape[0]
    
        k = np.dot(np.dot(np.ones(s),inv(I - lamda*XY)),np.ones(shape=(1,s))]))
    elif(method_type == "sylvester"):
        # algorithm presented in [Vishwanathan et al., 2006]
        # complexity of O(|V|^3)
        
        # Doing some initialization
        X_dimension = X.shape[0]
        Y_dimension = Y.shape[0]

        # For efficiency reasons multiply lambda
        # with the smallest e_{x,y} in dimension
        e_x = np.ones(X_dimension)
        e_y = np.ones(shape(1,Y_dimension))
        if(X_dimension < Y_dimension):
            e_x = np.dot(e_x,-lamda)
        else:
            e_y = np.dot(e_y,-lamda)

        # Prepare parameters for sylvester equation
        A = Y
        B = np.divide(inv(X), -lamda)
        C = np.dot(e_x ,np.dot(e_y,X))

        R = solve_sylvester(A, B, C)
        
        # calculate kernel
        k = - np.sum(np.sum(R,axis=0),axis=1)
    else:
        # raise exception
        # such method does not exist

    return k

def shortest_path(X, Y):
    """ This is a function that implements shortest path kernel
        as proposed by [Borgwardt & Kriegel, 2005]

        X, Y: to valid graph types i.e. edge dictionary or adjacency_mat
    """
    g_x = graph(X)
    g_y = graph(Y)

    S_x, L_x = g_x.build_shortest_path_matrix()
    S_y, L_y = g_y.build_shortest_path_matrix()

    return shortest_path_inner(S_x, S_y, L_x, L_y)

def shortest_path_inner(S_x, S_y, L_x={}, L_y={}):
    """ A function that calculates the inner summation part of the
        shortest path kernel as proposed by [Borgwardt & Kriegel, 2005]

        S_x, S_y: shortest path matrices for original graphs X, Y
        L_x, L_y: labels corresponding to index of S_x, S_y
                --note: this dictionary must be a surjection (one to many)

        Complexity: O(n^2*min(1,m)^2), where the number of vertices and m the mean maximum
        branching factor of the surjective labeling function of the Graph with minimum product
        over the two (not yet implemented!).
    """

    kernel = 0
    if(L_x == {} and L_y == {}):
    # Both non labeled case
        m = min(S_x.shape[0], S_y.shape[0])
        for i in range(0,m):
            for j in range(0,m):
                if(i!=j and (S_x[i,j] == S_y[i,j])):
                    kernel+=1
    elif(L_x != {} and L_y != {}):
        # Both labeled case
        L_y_inv = inv_dict(L_y)
        dim = S_x.shape[0]
        for i in range(0,dim):
            for j in range(0, dim):
                if(i!=j):
                    if((L_x[i] in L_y_inv) and (L_x[j] in L_y_inv)):
                        for k in L_y_inv[L_x[i]]:
                            for m in L_y_inv[L_x[j]]:
                                if(S_x[i,j] == S_y[k,m])
                                    kernel +=1
    else:
        # Case where only the one is labeled
        # Raise warning?
    return kernel

def inner_subtree_RG_dynamic(u, v, g_x, g_y, h, dynamic_dict, p_u=None, p_v=None)
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
		return (g_x.label(u)==g_y.label(v))
	elif h>1:
		R = list()
		
		
		# Calculate R
		# First group all nodes with the same label
		# make a list of lists of all pairs for each label
		lbx = g_x.label_group()
		lby = g_y.label_group()
		Rset = []
		# What happens when there are no neighbours?
		nx = g_x.neighbours(u)
		
		# Avoid going back
		if (p_u is not None):
			if p_u in nx:
				nx.remove(p_u)

		# do the same for y
		ny = g_y.neighbours(v)

		# Avoid going back
		if (p_v is not None):
			if p_v in ny:
				ny.remove(p_v)

		# What happens when one list has no neighbours?
		xy_and = len(nx)*len(ny)
		xy_or = len(nx)+len(ny)
		if(xy_or == 0):
			# both trees finish at the same point
			return (g_x.label(u)==g_y.label(v))
		elif(xy_and == 0):
			# else trees are different so output zero (?)
			return 0
			
		# returns a list
		snx = set(nx)
		sny = set(ny)
		for kx in lbx:
			if (kx is in lby):
				# substract from the common label 
				# only the valid neighbors
				snxk = list(set(lbx[k]).intersect(snx))
				snyk = list(set(lbx[y]).intersect(sny))
				if len(snxk) > 0 and len(snyk) > 0:
					pair = [lbx[kx],lby[ky]]
					Rset += list(itertools.product(*pair))

		# Designate all nodes with the same start
		# and all with the same finish
		right = dict()
		left = dict()
		l = 0
		# assign values to indexes
		for (w,z) in Rset:
			if w not in right:
				right[w] = list()
			if z not in left:
				left[z] = list()
			right[w].append(l)
			left[z].append(l)
			l+=1


		# Flatten both $right, $left to list of lists
		Kbins_flat_r = list()
		Kbins_flat_l = list()
		for k in right.keys():
			Kbins_flat_r.append(list(Kbins_flat_r[k]))
		for k in left.keys():
			Kbins_flat_l.append(list(Kbins_flat_l[k]))
		# calculate combinations
		# product equals combinations
		# because the lists are disjoint
		setA = set(itertools.product(*Kbins_flat_r))
		setB = set(itertools.product(*Kbins_flat_l))
		# take the intersect
		# found maximal valid R sets
		# allsubsets all also valid
		Rmaximal = setA.intersection(setB)

		# Rset: set of all pairs in R
		l = 0
		# Rv assign values best on index
		Rv = dict()
		for (w,wp) in Rset:
			r = nested_dict_get(dynamic_dict, k, u, v, w, wp)
			if (r is None):
				kh = inner_subtree_RG_dynamic(w, wp, g_x, g_y, h-1, dynamic_dict, u, v)
				nested_dict_add(dynamic_dict, kh, u, v, w, wp, h-1)
				Rv[l] = kh
			else:
				Rv[l] = r
			
		
		k = 0	
		for s in Rmaximal:
			# Can we calculate do calculation on subsets faster?
			# The answer is yes: dynamic programming
			non_zero_elements = set()
			for i in s:
				if(Rv[l] != 0):
					non_zero_elements.add()
			# Dictionary that keeps the subset values if they have been calculated
			# keys of type set
			value = dict()
			plough_subsets(ns.copy().discard(i), Rv, value):
			for v in value.values():
				k+=v

		return k
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
		value[frozenset(initial_set)] = Rv[initial_set.pop()]
	
	elif (len(initial_set) > 1):
		# If all elements explore all subsets
		# for value dictionary to fill and
		# take 
		flag = True
		for s in initial_set:
			temp_set = ns.copy().discard(s)
			frozen_temp_set = frozenset(temp_set)
			if q not in value:
				if(flag):
					value[frozenset(initial_set)] = Rv[s]*value[q]
				plough_for_subsets(temp_set, Rv, value)			
			else:
				if(flag):
					value[frozenset(initial_set)] = Rv[s]*value[q]
								flag = False
				
def subtree_RG(X_edge_dict, Y_edge_dict, L_x={}, L_y={}, h=5):
	""" Calculate The Ramon Gartner subtree kernel
		as proposed on [Ramon & Gartner, 2003]

		X_edge_dict, Y_edge_dict: The two graphs
			whose kernel is going to be computed
			as a dictionary of edges as follows:
				If (u,v) edge exists then edges["u"]["v"] has
                weight of the edge pointing from u to v

		X, Y: to valid graph types i.e. edge dictionary or adjacency_mat

		height: The height of the subtree exploration		 
	"""
	
	g_x = graph(X,L_x)
	g_y = graph(Y,L_y)

	kernel = 0
	dynamic_dictionary = dict()
	for u in X.edge_dictinary.keys():
		for v in Y.edge_dictionary.keys():
			inner_subtree_RG_dynamic(u,v,g_x,g_y,h,dynamic_dictionary)
