"""
    The main graph kernel object used for creating a graph kernel object
"""
import warnings

import numpy as np

from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.utils.validation import check_array

from .graph import graph
from .kernels import *

global supported_base_kernels, supported_general_kernels

supported_base_kernels = [
"dirac","random_walk",
"shortest_path","subtree_rg",
"graphlets_sampling","subgraph_matching",
"multiscale_laplacian",
"lovasz_theta", "svm_theta",
"neighborhood_hash", "neighborhood_pairwise_subgraph_distance",
"odd_sth"
]
    
supported_general_kernels = [
"weisfeiler_lehman"
]


class GraphKernel(BaseEstimator, TransformerMixin):
    """ A general class that describes all kernels

    Parameters
    ----------
    kernel : a list of dictionaries, or a single dicitonary
             that have the following structure:
                        "name": [str] - with the kernel name
                        "name_of_parameter_1": value
                        "name_of_parameter_2": value
                                ...
                        "name_of_parameter_k": value
             available "names" | "parametres" are:

             i)  base_kernels (the structure must always reach a base kernel):
                    - "dirac" | 

                    - "random_walk" | (o) "lamda": [float] < 1 , (o) method_type : [str] "sylvester"
                                                                                 "simple"
                    - "shortest_path" | (o) "algorithm_type": [str] "dijkstra"
                                                                "floyd_warshall"
                    - "subtree_rg" | (o) "h": [int]
                    
                    - "graphlets_sampling" | (o) "k": [int], (o) "delta" : [float], (o) "epsilon" : [float], (o) "a": [int]
                    
                    - "multiscale_laplacian" | (o) "L": [int], (o) "gamma"
                    
                    - "subgraph_matching" | (o) "kv": kernel function for nodes
                                                      (node_x, node_y, Lx, Ly) -> number
                                            (o) "ke": kernel function for edges
                                                      (edge_x, edge_y, Lx, Ly) -> number
                                            (o) "lw": a lambda weight function for cliques
                                                      set -> number
                    - "lovasz_theta" | (o) "n_samples": [int] > 1
                                       
                                       (o) "subsets_size_range": [touple] of two [int]
                                       
                                       (o) "metric": [function]
                                                     (number, number) -> number
                                       
                    - "svm_theta" | (o) "n_samples": [int] > 1
                                    
                                    (o) "subsets_size_range": [touple] with 2 [int] elements
                                    
                                    (o) "metric": [function]
                                                  (number, number) -> number
                                                  
                    - "neighborhood_hash" | (o) "nh_type": [str] "simple" or "count-sensitive"
                    
                                            (o) "R": [int] > 0
                                            
                                            (o) "bytes": [int] > 0
                                            
                    - "neighborhood_pairwise_subgraph_distance" | (o) "r": (int) positive integer
                    
                                                                  (o) "d": (int) positive integer
                    - "odd_sth" | (o) "h": (int) positive integer
                    
             ii) general_kernels (this kernel will consider the next kernel on list for nesting):
                    - "weisfeiler_lehman": "niter": [int]
                
                where (o): stands for optional parameters
                
    Attributes
    ----------
    X_graph: dictionary 
        Stores input graphs indexing from 0 to num__of_graphs-1
    num_of_graphs: int
        The number of graphs given in fit input.
    kernel: the full kernel applied between graph objects
        
        
    """

    num_of_graphs = 0
    X_graph = None

    def __init__(self, **kargs):
        if "kernel" in kargs:
            if (type(kargs["kernel"]) is dict):
                # allow single kernel dictionary inputs
                kernel, flag = self._make_kernel([kargs["kernel"]])
            elif (type(kargs["kernel"]) is list):
                kernel, flag = self._make_kernel(kargs["kernel"])
            else:
                raise ValueError('unsupported kernel format')
            if flag:
                self.kernel = lambda y: self.calculate_kernel_matrix(kernel, target_graph=y, kernel_type="matrix")
            else:
                self.kernel = lambda y: self.calculate_kernel_matrix(kernel, target_graph=y, kernel_type="pairwise")
        else:
            raise ValueError('kernel must be defined at the __init__ function of a graph kernel')

    def _make_kernel(self,kernel_list):
            """ Produces the desired kernel function.
                
                Parameters
                ----------
                kernel_list: list of kernel dictionaries 
                as defined at the documentation of class parameters
                
                Output
                ------
                kernel: function (of two arguments)
                        returns the kernel function
            """
            # If nesting type:
            kernel = kernel_list.pop(0)
            kernel_name = kernel.pop("name")
            if kernel_name in supported_base_kernels:
                if (len(kernel_list)!=0):
                    warnings.warn('rest kernel arguments are being ignored - reached base kernel')
                if kernel_name == "dirac":
                    return (lambda x, y: dirac_inner(x, y), False)
                elif kernel_name == "random_walk":
                    return (lambda x, y: random_walk_inner(x, y, **kernel), False)
                elif kernel_name == "shortest_path":
                    return (lambda x, y: shortest_path_inner(x, y, **kernel), False)
                elif kernel_name == "subtree_rg":
                    return (lambda x, y: subtree_rg_inner(x, y, **kernel), False)
                elif kernel_name == "graphlets_sampling":
                    nsamples, graphlets, P, graph_bins, nbins = sample_graphlets(**kernel)
                    print(str(nsamples)+" graphlets sampled")
                    if "k" in kernel:
                        return (lambda x, y: graphlet_sampling_core(x, y, nsamples, graphlets, P, graph_bins, nbins, k=kernel["k"]), False)
                    else:
                        return (lambda x, y: graphlet_sampling_core(x, y, nsamples, graphlets, P, graph_bins, nbins, k), False)
                elif kernel_name == "multiscale_laplacian":
                    return (lambda x, y: multiscale_laplacian_inner(x, y, **kernel), False)
                elif kernel_name == "subgraph_matching":
                    return (lambda x, y: subgraph_matching_inner(x, y, **kernel), False)
                elif kernel_name == "lovasz_theta":
                    return (lambda x, y: lovasz_theta_inner(x, y, **kernel), False)
                elif kernel_name == "svm_theta":
                    return (lambda x, y: svm_theta_inner(x, y, **kernel), False)
                elif kernel_name == "neighborhood_hash":
                    return (lambda x, y: neighborhood_hash_matrix(x, y, **kernel), True)
                elif kernel_name == "neighborhood_pairwise_subgraph_distance":
                    return (lambda x, y: neighborhood_pairwise_subgraph_distance_inner(x, y, **kernel), False)
                elif kernel_name == "odd_sth":
                    return (lambda x, y: odd_sth_matrix(x, y, **kernel), True)
            elif kernel_name in supported_general_kernels:
                if (len(kernel_list)==0):
                    raise ValueError(str(kernel_name)+' is not a base kernel')
                if kernel_name == "weisfeiler_lehman":
                    (bk, matrix_flag) = self._make_kernel(kernel_list)
                    if matrix_flag:
                        return (lambda x, y: weisfeiler_lehman_matrix(x, bk, y, **kernel), True)
                    else:
                        kernel["base_kernel"] = bk
                        return (lambda x, y: weisfeiler_lehman_inner(x,y,**kernel), False)
            else:
                raise ValueError('unsupported kernel: '+str(kernel_name))
                
    def fit(self, X, y=None):
        """A reference implementation of a fitting function for a transformer.

        Parameters
        ----------
        X : array-like or sparse matrix of shape = [n_samples, n_features]
            There must be two features: A valid graph structure (adjacency matrix
            or edge_dictionary, vertex_labels, node_labels) 
            The training input samples.
        y : None
            There is no need of a target in a transformer, yet the pipeline API
            requires this parameter.

        Returns
        -------
        self : object
            Returns self.
        """
        # Input validation
        # check_array(X)
        
        self.X_graph = dict()
        graph_idx = 0        
        for x in list(X):
            if len(x) == 1:
                self.X_graph[graph_idx] = graph(x[0],{},{},"all")
                graph_idx += 1
            elif len(x) == 2:
                self.X_graph[graph_idx] = graph(x[0],x[1],{},"all")
                graph_idx += 1
            elif len(x) == 3:
                self.X_graph[graph_idx] = graph(x[0],x[1],x[2],"all")
                graph_idx += 1 
            else:
                warnings.warn('input has an empty or unrecognised element')
        if not bool(self.X_graph):
            raise ValueError('unrecognized input')
        self.num_of_graphs = graph_idx

        # Return the transformer
        return self

    def transform(self, X=None):
        """ The transform function calculates the kernel matrix

        Parameters
        ----------
        X : array-like matrix of shape = [n_samples, n_features]
            There must be two features: A valid graph structure (adjacency matrix
            or edge_dictionary, edge_labels). If None the kernel matrix is calculated
            upon fit data.
            The training input samples.The input samples.

        Returns
        -------
        K : A numpy array of shape = [n_targets, n_input_graphs] corresponding
            to the kernel matrix, a calculation between all pairs of graphs 
            between target an features
        """
        # Check is fit had been called
        # check_is_fitted(self, ['X_graph_'])
        
        # Input validation
        # check_array(X)
        
        if(X is not None):
            target_graph = dict()
            num_of_targets = 0
            for x in X:
                if len(x) == 1:
                    target_graph[num_of_targets] = graph(x[0],{},{},"all")
                    num_of_targets += 1
                elif len(x) == 2:            
                    target_graph[num_of_targets] = graph(x[0],x[1],{},"all")
                    num_of_targets += 1
                elif len(x) == 3:
                    self.X_graph[graph_idx] = graph(x[0],x[1],x[2],"all")
                    num_of_targets += 1
                else:
                    warnings.warn('input has an empty or unrecognised element')
            if not bool(target_graph):
                raise ValueError('unrecognized input')
        else:
             target_graph = None           
        # Check that the input is of the same shape as the one passed
        # during fit.

        # Calculate kernel matrix

        return self.kernel(target_graph)
        
    def calculate_kernel_matrix(self, kernel, target_graph=None, kernel_type="pairwise"):
        """
        A function that calculates the kernel matrix given
        a target_graph and a kernel
        
        Parameters
        ----------
        target_graph : A dictionary from 0 to the number of graphs
                       of target "graph objects"
        kernel: A pairwise graph kernel (between "graph" type objects)
        kernel_type: [str] "matrix", "pairwise"
        
        Returns
        -------
        K : A numpy array of shape = [n_targets, n_input_graphs] corresponding
            to the kernel matrix, a calculation between all pairs of graphs 
            between target an features
        """
        if kernel_type not in ["matrix", "pairwise"]:
            raise ValueError('unsupported "kernel_type"')
        
        if kernel_type == "pairwise":
            if target_graph is None:
                target_graph = self.X_graph
                is_symmetric = True
                num_of_targets = self.num_of_graphs
            else:
                is_symmetric = False
                num_of_targets = len(target_graph.keys())
                
            K = np.zeros(shape = (num_of_targets,self.num_of_graphs))
            if is_symmetric:
                for i in range(0,num_of_targets):
                    for j in range(i,self.num_of_graphs):
                        K[i,j] = kernel(target_graph[i],self.X_graph[i])
                        if(i!=j):
                            K[j,i] = K[i,j]
            else:
                for i in range(0,num_of_targets):
                    for j in range(0,self.num_of_graphs):
                        K[i,j] = kernel(target_graph[i],self.X_graph[i])
            return K     
        else:
            return kernel(self.X_graph, target_graph)
