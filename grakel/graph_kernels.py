"""
    The main graph kernel object used for creating a graph kernel object
"""

import numpy as np

from sklearn.base import BaseEstimator, ClassifierMixin, TransformerMixin
from sklearn.utils.validation import check_X_y, check_array, check_is_fitted

from .graph import graph
from .kernels import dirac_inner, random_walk_inner, shortest_path_inner, subtree_rg_inner, sample_graphlets, graphlet_sampling_core, weisfeiler_lehman_inner

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
                self.kernel = self._make_kernel([kargs["kernel"]])
            elif (type(kargs["kernel"]) is list):
                self.kernel = self._make_kernel(kargs["kernel"])
        else:
            pass
            # Raise error Undefined kernel ?

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
            if kernel_name in ["dirac","random_walk","shortest_path","subtree_rg","graphlets_sampling"]:
                if (len(kernel_list)!=0):
                    pass
                    # Raise Warning: Arguments are being ignored reached base kernel?
                if kernel_name is "dirac":
                    return lambda x, y: dirac_inner(x,y)
                elif kernel_name is "random_walk":
                    return lambda x, y: random_walk_inner(x,y,**kernel)
                elif kernel_name is "shortest_path":
                    return lambda x, y: shortest_path_inner(x,y,**kernel)
                elif kernel_name is "subtree_rg":
                    return lambda x, y: subtree_rg_inner(x,y,**kernel)
                elif kernel_name is "graphlets_sampling":
                    nsamples, graphlets, P, graph_bins, nbins = sample_graphlets(**kernel)
                    print(str(nsamples)+" graphlets sampled")
                    if "k" in kernel:
                        return lambda x, y: graphlet_sampling_core(x,y,nsamples, graphlets, P, graph_bins, nbins, k=kernel["k"])
                    else:
                        return lambda x, y: graphlet_sampling_core(x,y,nsamples, graphlets, P, graph_bins, nbins, k)
            elif kernel_name in ["weisfeiler_lehman"]:
                if (len(kernel_list)==0):
                    pass
                    # Raise Error: $kernel_name is not a base kernel?
                if kernel_name is "weisfeiler_lehman":
                    bk = self._make_kernel(kernel_list)
                    kernel["base_kernel"] = bk
                    return lambda x, y: weisfeiler_lehman_inner(x,y,**kernel)

    def fit(self, X, y=None):
        """A reference implementation of a fitting function for a transformer.

        Parameters
        ----------
        X : array-like or sparse matrix of shape = [n_samples, n_features]
            There must be two features: A valid graph structure (adjacency matrix
            or edge_dictionary, edge_labels) 
            The training input samples.
        y : None
            There is no need of a target in a transformer, yet the pipeline API
            requires this parameter.

        Returns
        -------
        self : object
            Returns self.
        """

        self.X_graph = dict()
        graph_idx = 0
        for x in list(X):
            if len(x) == 1:
                self.X_graph[graph_idx] = graph(x[0],{},"all")
                graph_idx += 1
            elif len(x) == 2:
                self.X_graph[graph_idx] = graph(x[0],x[1],"all")
                graph_idx += 1
            else:
                pass
                # Raise warning input hole.?
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
        #check_is_fitted(self, ['input_shape_'])

        # Input validation
        # ~X = check_array(X)
        if(X is None):
            is_symmetric = True
            target_graph = self.X_graph
            num_of_targets = self.num_of_graphs
        else:
            is_symmetric = False
            target_graph = dict()
            num_of_targets = 0
            for x in X:
                if len(x) == 1:
                    target_graph[num_of_targets] = graph(x[0],{},"all")
                    num_of_targets += 1
                elif len(x) == 2:            
                    target_graph[num_of_targets] = graph(x[0],x[1],"all")
                    num_of_targets += 1
                else:
                    pass
                    # Raise warning passing null element?

        # Check that the input is of the same shape as the one passed
        # during fit.

        #if X.shape != self.input_shape_:
        #    raise ValueError('Shape of input is different from what was seen'
        #                     'in `fit`')
        
        # Calculate kernel matrix
        K = np.zeros(shape = (num_of_targets,self.num_of_graphs))
        if is_symmetric:
            for i in range(0,num_of_targets):
                for j in range(i,self.num_of_graphs):
                    K[i,j] = self.kernel(target_graph[i],self.X_graph[i])
                    if(i!=j):
                        K[j,i] = K[i,j]
        else:
            for i in range(0,num_of_targets):
                for j in range(0,self.num_of_graphs):
                    K[i,j] = self.kernel(target_graph[i],self.X_graph[i])
        return K
