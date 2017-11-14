"""
This is a module to be used as a reference for building other modules
"""
import numpy as np
from sklearn.base import BaseEstimator, ClassifierMixin, TransformerMixin
from sklearn.utils.validation import check_X_y, check_array, check_is_fitted
from graph import graph
from kernels import dirac_inner, random_walk_inner, shortest_path_inner, subtree_RG_inner, graphlets_sampling_inner, weisfeiler_lehman_inner
 
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
                    - "subtree_RG" | (o) "h": [int]
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
    def __init__(self, **kargs):
        if "kernel" in kargs:
            if (type(kargs["kernel"]) is dict):
                # allow single kernel dictionary inputs
                self.kernel = _make_kernel([kargs["kernel"]])
            elif (type(kargs["kernel"]) is list):
                self.kernel = _make_kernel(kargs["kernel"])
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
            if kernel_name in ["dirac","random_walk","shortest_path","subtree_RG","graphlets_sampling"]:
                if (len(kernel_list)!=0):
                    pass
                    # Raise Warning: Arguments are being ignored reached base kernel?
                if kernel_name is "dirac":
                    return lambda x, y: dirac_inner(x,y)
                elif kernel_name is "random_walk":
                    return lambda x, y: random_walk_inner(x,y,**kernel)
                elif kernel_name is "shortest_path":
                    return lambda x, y: shortest_path_inner(x,y,**kernel)
                elif kernel_name is "subtree_RG":
                    return lambda x, y: subtree_RG_inner(x,y,**kernel)
                elif kernel_name is "graphlets_sampling":
                    return lambda x, y: graphlets_sampling_inner(x,y,**kernel)
            elif kernel_name in ["weisfeiler_lehman"]:
                if (len(kernel_list)==0):
                    pass
                    # Raise Error: $kernel_name is not a base kernel?
                if kernel_name is "weisfeiler_lehman":
                    bk = _make_kernel(self,kernel_list)
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
        for x in X:
            self.X_graph[graph_idx] = graph(x[0],x[1],"all")
            graph_idx += 1
        self.num_of_graphs = graph_idx

        # Return the transformer
        return self

    def transform(self, X=None):
        """ The transform function calculates the kernel matrix

        Parameters
        ----------
        X : array-like or sparse matrix of shape = [n_samples, n_features]
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
        check_is_fitted(self, ['input_shape_'])

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
                target_graph[num_of_targets] = graph(x[0],x[1],"all")
                num_of_targets += 1

        # Check that the input is of the same shape as the one passed
        # during fit.

        #if X.shape != self.input_shape_:
        #    raise ValueError('Shape of input is different from what was seen'
        #                     'in `fit`')
        
        # Calculate kernel matrix
        K = np.empty(shape = (self.num_of_targets,self.num_of_graphs))
        if is_symmetric:
            for i in range(0,self.num_of_targets):
                for j in range(i,self.num_of_graphs):
                    K[i,j] = self.kernel(target_graph[i],self.X_graph[i])
                    if(i!=j):
                        K[j,i] = K[i,j]
        else:
            for i in range(0,self.num_of_targets):
                for j in range(0,self.num_of_graphs):
                    K[i,j] = self.kernel(target_graph[i],self.X_graph[i])
        return K
