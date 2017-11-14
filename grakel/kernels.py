""" This file contains a set of functions implementing functions implementing
    various graph kernel metrics.

"""
import itertools
import math

import numpy as np

import pynauty

from np.linalg import inv
from scipy.linalg import solve_sylvester
from scipy.interpolate import interp1d

from graph import graph
from tools import inv_dict, matrix_to_dict

random.seed(352163)
np.random.seed(238537)

def random_walk(X, Y, Lx, Ly, lamda=0.1, method_type="sylvester"):
    """ This is a function that implements the random walk kernel.

        X, Y: to valid graph types i.e. edge dictionary or adjacency_mat
        algorithm_type: "dijkstra" or "floyd_warshall" or "auto" (for the best concerning current format)
    """
    g_x = graph(X)
    g_y = graph(Y)

    return random_walk_inner(X, Y,lamda, method_type)


def random_walk_inner(Gx, Gy,lamda=0.1, method_type="sylvester"):
    """ This is a function that implements the random walk kernel.

        Gx, Gy: are two graph type objects representing relevant graphs to 
        be compared. Their format is a s    uare array.
        lamda: the factor concerning summation
        method_type: "simple" [O(|V|^6)] or "sylvester" [O(|V|^3)]
    """
    Gx.desired_type("adjacency")
    Gy.desired_type("adjacency")
    X = Gx.adjacency_matrix
    Y = Gy.adjacency_matrix
    
    # calculate the product graph
    XY = np.kron(A,B)

    if(method_type == "simple"):
        # algorithm presented in [Kashima et al., 2003; Gärtner et al., 2003]
        # complexity of O(|V|^6)

        # XY is a s    uare matrix
        s = XY.shape[0]
    
        k = np.dot(np.dot(np.ones(s),inv(I - lamda*XY)),np.ones(shape=(1,s)))
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

        # Prepare parameters for sylvester e    uation
        A = Y
        B = np.divide(inv(X), -lamda)
        C = np.dot(e_x ,np.dot(e_y,X))

        R = solve_sylvester(A, B, C)
        
        # calculate kernel
        k = - np.sum(np.sum(R,axis=0),axis=1)
    else:
        pass
        # raise exception?
        # such method does not exist

    return k

def shortest_path(X, Y, Lx, Ly, algorithm_type="dijkstra"):
    """ This is a function that implements shortest path kernel
        as proposed by [Borgwardt & Kriegel, 2005]

        X, Y: to valid graph types i.e. edge dictionary or adjacency_mat
        algorithm_type: "dijkstra" or "floyd_warshall" or "auto" (for the best concerning current format)
    """
    g_x = graph(X)
    g_y = graph(Y)

    return shortest_path_inner(S_x, S_y, Lx, Ly, algorithm_type)

def shortest_path_inner(g_x, g_y, algorithm_type="dijkstra"):
    """ A function that calculates the inner summation part of the
        shortest path kernel as proposed by [Borgwardt & Kriegel, 2005]

        g_{x,y}: graph type objects
        algorithm_type: "dijkstra" or "floyd_warshall" or "auto" (for the best concerning current format)

        Complexity: O(n^2*min(1,m)^2), where the number of vertices and m the mean maximum
        branching factor of the surjective labeling function of the Graph with minimum product
        over the two (not yet implemented!).
    """
    S_x, Lx = g_x.build_shortest_path_matrix(algorithm_type)
    S_y, Ly = g_y.build_shortest_path_matrix(algorithm_type)

    kernel = 0
    Ly_inv = inv_dict(Ly)
    dim = S_x.shape[0]
    for i in range(0,dim):
        for j in range(0, dim):
            if(i!=j):
                if((Lx[i] in Ly_inv) and (Lx[j] in Ly_inv)):
                    for k in Ly_inv[Lx[i]]:
                        for m in Ly_inv[Lx[j]]:
                            if(S_x[i,j] == S_y[k,m]):
                                kernel +=1

    return kernel

def subtree_RG(X, Y, Lx, Ly, h=5):
    """ Computes the Weisfeler Lehman as proposed
        at 2011 by shervashize et. al.

        X, Y: Valid graph formats to be compared
        Lx, Ly: Valid labels for graphs
        base_kernel: A valid base kernel
    """
    Gx = graph(X,Lx)
    Gy = graph(Y,Ly)
    return subtree_RG_inner(Gx, Gy, h=5)

def subtree_RG_inner(Gx, Gy, h=5):
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
    Gx.desired_fromat("dictionary")
    Gy.desired_fromat("dictionary")
    kernel = 0
    dynamic_dictionary = dict()
    for u in X.edge_dictinary.keys():
        for v in Y.edge_dictionary.keys():
            kernel += core_subtree_RG_dynamic(u,v,Gx,Gy,h,dynamic_dictionary)
    return kernel

def core_subtree_RG_dynamic(u, v, g_x, g_y, h, dynamic_dict, p_u=None, p_v=None):
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
            if (kx in lby):
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
        # product e    uals combinations
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
            r = nested_dict_get(dynamic_dict, h-1, u, v, w, wp)
            if (r is None):
                kh = inner_subtree_RG_dynamic(w, wp, g_x, g_y, h-1, dynamic_dict, u, v)
                nested_dict_add(dynamic_dict, kh, h-1, u, v, w, wp)
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
            plough_subsets(ns.copy().discard(i), Rv, value)
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
        # Explore all subsets for value dictionary to fill and
        # and calculate current
        flag = True
        frozen_inital_set = frozenset(initial_set)
        for s in initial_set:
            temp_set = ns.copy().discard(s)
            frozen_temp_set = frozenset(temp_set)
            if frozen_temp_set not in value:
                plough_for_subsets(temp_set, Rv, value)            
            if(flag):
                value[frozenset_initial_set] = Rv[s]*value[frozen_temp_set]
                flag = False

def graphlets_sampling(X, Y, k = 5, delta=0.05, epsilon=0.05, a=-1):
    """ Applies the sampling random graph kernel as proposed
        by Shervashidze, Vishwanathan at 2009 (does not consider labels)
        
        X,Y: valid graph formats

    """
    Gx = graph(X,{})
    Gy = graph(Y,{})
    return graphlets_sampling_inner(Gx, Gy, k = 5, delta=0.05, epsilon=0.05, a=-1)
    
def graphlets_sampling_inner(Gx, Gy, k = 5, delta=0.05, epsilon=0.05, a=-1):
    """ Applies the sampling random graph kernel as proposed
        by Shervashidze, Vishwanathan at 2009 (does not consider labels)
        
        Gx, Gy: Graph type objects
        k: the dimension of the given graphlets
        delta : confidence level (typically 0.05 or 0.1)
        epsilon : precision level (typically 0.05 or 1)
        a : number of isomorphism classes of graphlets
    """
    Gx.desired_format("adjacency")
    Gy.desired_format("adjacency")
    X, Y = Gx.adjacency_matrix, Gy.adjacency_matrix

    fallback_map = {1: 1, 2: 2, 3: 4, 4: 8, 5: 19, 6: 53, 7: 209, 8: 1253, 9: 13599}
    if not k>=3:
        # Raise Warning
        exit(1)
    if(a==-1):
        if(k>9):
            # Raise warning for such size a is not known.
            # Use interpolations
            isomorphism_prediction = interp1d(list(fallback_map.keys()), list(fallback_map.values()), kind='cubic')
            a = isomorphism_prediction(k)
        else:
            a = fallback_map[k]
    
    nsamples = math.floor(2 * ( a* np.log10(2) + np.log10(1/delta) ) / (epsilon*epsilon))

    # Steps:
    # 1) Sample graphlets
    graphlets = dict()
    graphlet_set = set()
    # keeps a track of what 
    # to store graphlets as pynauty graphs
    to_edge_dict_binary = lambda x : matrix_to_dict(x, '==', 1, k, False)
    graph_bins = dict()
    nbins = 0
    for i in range(0,nsamples):
        # An upper triangular matrix is calculated
        # characteristics: at least one zero in every line
        # of the mat so in the upper corresponding
        # line and row.
        gr = np.empty(shape = (k,k))
        cert = list()
        f = True
        while f:
            for i in range(0,k-1):
                gr[i,i] = .0
                line = np.random.randint(2, size=k-1)
                while 1 not in line[:]:
                    line = np.random.randint(2, size=k-1)
                gr[0:(i-1),i] = line[0:i-1]
                gr[i,(i+1):k-1] = line[i:k-2]
                # Apply also for symmetric
                gr[(i+1):k,i] = line[i:k-2]
                gr[i,0:(i-1)] = line[0:i-1]
                cert += list(line[i:k-2])
            certificate = str(cert)
            f = certificate in graphlet_set
        graphlet_set.add(str(cert))
        graphlets[i] = pynauty.Graph(k, True, to_edge_dict_binary(gr))
        if i==0:
            graph_bins[0] = [0]
            nbins+=1
        else:
            newbin = True
            for j in range(nbins):
                if pynauty.isomorphic(graph_bins[j],graphlets[i]):
                    newbin = False
                    graph_bins[j].append(i)
                    break
            if newbin:
                graph_bins[nbins] = [i]
                nbins+=1

    # Produce Pij Matrix:
    # Based on the idea that 
    # if Pij = 1 and Pjk = 1 then Pik=1
    P = np.zeros(nsamples)
    for i in range(nbins):
        pair = list(graph_bins[i])
        for (i,j) in itertools.combinations(pair,2):
            P[i,j]=1

    # also Pij = Pji
    for i in range(0,nsamples):
        for j in range(i+1,nsamples):
            P[j,i]=P[j,i]

    # 2) Calculate fre    uencies for each graph
    # Check that if each matrix is a principal
    # minor of the adjoint and how many times
    i = 0

    # Fre    uency vector for x
    fx = np.zeros(nsamples)
    
    # Fre    uency vector for y
    fy = np.zeros(nsamples)

    # To transform the adjacency to edge dictionary
    # needed for nauty graph initialise a small lambda
    to_edge_dict_real = lambda x : matrix_to_dict(x, '>', .0, k, False)
    for c in itertools.combinations(list(range(k)),k):
        for s in range(0,nsamples):
            idxs = list(c)
            if(pynauty.isomorphic(pynauty.Graph(k, True, to_edge_dict_real(X[c,c])), graphlets[s])):
                fx[s]+=1
            if(pynauty.isomorphic(pynauty.Graph(k, True, to_edge_dict_real(Y[c,c])),graphlets[s])):
                fy[s]+=1

    # normalize fx
    sfx = np.sum(fx,axis=0)
    if(sfx != 0):
        fx = np.divide(fx,sfx)
    
    # normalize fy
    sfy = np.sum(fy,axis=0)
    if(sfy != 0):
        fy = np.divide(fy,sfy)

    # 3) Calculate the kernel
    kernel = np.dot(fx,np.dot(P,fy.T))

    return kernel

def weisfeiler_lehman(X, Y, Lx, Ly, base_kernel, niter):
    """ Computes the Weisfeler Lehman as proposed
        at 2011 by shervashize et. al.

        X,Y: Valid graph formats to be compared
        Lx, Ly: Valid labels for graphs
        base_kernel: A valid base kernel
    """
    Ga = graph(X,Lx)
    Gb = graph(Y,Ly)
    return weisfeiler_lehman_inner(Ga, Gb, niter)

def weisfeiler_lehman_inner(Ga, Gb, niter, base_kernel):
    """ Computes the Weisfeler Lehman as proposed
        at 2011 by shervashize et. al.

        Ga,Gb: Graph type formats of the graphs
               to be compared
        Lx, Ly: Valid labels for graphs

        base_kernel: A valid kernel function of the form:
        graph_A, graph_B -> number
    """
    Ga.desired_format("dictionary")
    Gb.desired_format("dictionary")
    
    La = Ga.get_labels("dictionary")
    Lb = Gb.get_labels("dictionary")
    Ga_edge_dictionary = Ga.edge_dictionary
    Gb_edge_dictionary = Gb.edge_dictionary

    WL_labels = OrderedDict()
    WL_labels_inverse = OrderedDict()
    distinct_values = sorted(list(set(La.values()).union(set(Lb.values()))))
    label_count = 0
    for dv in distinct_values:
        WL_labels[label_count] = dv
        WL_labels_inverse[dv] = label_count
        label_count +=1
    for k in La.keys():
        La[k] = WL_labels_inverse[La[k]]
    for k in Lb.keys():
        Lb[k] = WL_labels_inverse[Lb[k]]

    k = base_kernel(Ga_edge_dictionary, Gb_edge_dictionary, La, Lb)
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

        # relabel
        for k in La_temp.keys():
            La[k] = WL_labels_inverse[k]
        for k in Lb_temp.keys():
            Lb[k] = WL_labels_inverse[k]

        # calculate kernel
        k += base_kernel(Ga_edge_dictionary, Gb_edge_dictionary, La, Lb)
    return k

def dirac(input_a, input_b, Lx, Ly):
    """ The simple dirac kernel for labelled graphs
        k(X,Y) = \sum_{v \in V_{1}}\sum_{v' \in V_{2}}dirac(l(v),l(v'))

        edge_dict_{X,Y}: The edge dicitonaries of the two compared graphs
        L{x,y}: The corresponding label dictionaries.
    """
    Gx = graph(edge_dict_X, edge_dict_Y, Lx, Ly)
    Gy = graph(edge_dict_X, edge_dict_Y, Lx, Ly)
    return dirac_inner(Gx,Gy)

def dirac_inner(Gx,Gy):
    """ The simple dirac kernel for labelled graphs
        k(X,Y) = \sum_{v \in V_{1}}\sum_{v' \in V_{2}}dirac(l(v),l(v'))

        G_{x,y}: corresponding graph structures
    """
    Gx.desired_format("dictionary")
    Gy.desired_format("dictionary")
    
    # Calculate kernel
    linv_x = inv_dict(Gx.get_labels("dictionary"))
    linv_y = inv_dict(Gy.get_labels("dictionary"))
    k = 0
    for lx in linv_x:
        if ly in linv_y:
            k += len(linv_x[lx])*len(linv_y[ly])     
    return k

#def graphlets(X,Y):
    """ Calculates the full graphlet kernel
        as proposed by Shervashidze and Vishwanathan (2009)
        for graphlets of size 3,4,5
    """
