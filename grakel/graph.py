""" A python file that implements a class and functions for graphs.

"""
import warnings
import collections
import operator
import itertools
import numbers

import numpy as np
import scipy

import cvxopt.base
import cvxopt.solvers

from sklearn.svm import OneClassSVM
from scipy.sparse import csr_matrix

from .tools import priority_dict, inv_dict, distribute_samples, nested_dict_add

np.random.seed(238537)
cvxopt.solvers.options['show_progress'] = False

class graph(object):
    """ A general graph class that supports adjacency, dictionary formats while beeing memory/computationaly sustainable.
    
    Parameters
    ----------
    initialization_object :  dict, list or array-like, square, default=None
        The initialisation object for the graph [**valid-graph-format**]. 
        If given a dictionary the input can be as follows:
            + 2-level nested dictionaries from edge symbols to weights
            + Dictionary of symbols to list of symbols (unweighted)
            + Dictionary of tuples to weights (weighted)
            + Iterable of tuples of len 2 (unweighted)
            + Iterable of tuples of len 3 (vertex, vertex, weight)
        
        If given a array the input can be as follows:
            + array-like lists of lists
            + np.array
            + sparse matrix (scipy.sparse.csr.csr_matrix)
    
    node_labels : dict, default=None
        A label dictionary corresponding to all vertices of the graph:
            + for adjacency matrix labels should be given to numbers starting 
              from 0 and ending in N-1, where the matrix has size N by N
            + for dictionary labels should correspond to all keys
    
    edge_labels : dict, default=None
        A labels dictionary corresponding to all edges of the graph keys: tuples, value: label
    
    graph_format : str, valid_values={"dictionary", "adjacency", "all", "auto"}, default=None
        Defines the internal representation of the graph object which can be a dictionary as a matrix, or both:
          + for dictionary: "dictionary"
          + for adjacency_matrix: "adjacency"
          + for both: "all"
          + for the current_format (if existent): "auto"
          
    Attributes
    ----------
    n : int
        The one dimension of the (square) adjacency matrix.
        Signifies the number of vertices.
        
    adjacency_matrix: np.array
        The adjacency_matrix corresponding to the graph.
        
    index_node_labels : dict
        Label dictionary for indeces, of adjacency matrix.
        Keys are valid numbers from 0 to n-1.
        
    index_edge_labels : dict
        Label dictionary for edges, of adjacency matrix.
        Keys are tuple pairs of symbols with valid numbers from 0 to n-1, that have a positive adjacency matrix value.
        
    vertices: set
        The set of vertices corresponding to the edge_dictionary representation.
        
    edge_dictionary: dict
        A 2-level nested dictionary from edge symbols to weights.
        
    node_labels : dict
        Label dictionary for nodes.
        Keys are vertex symbols inside vertices.
        
    edge_labels : dict
        Label dictionary for edges.
        Keys are tuple pairs of symbols inside vertices and edges.
        
    edsamic : dict
        *Edge-Dictionary-Symbols-Adjacency-Matrix-Index-Correspondance*.
        A dictionary which translates dictionary symbols to adjacency matrix indexes, when storing both formats.
        
    shortest_path_mat : np.array, square
        Holds the shortest path matrix.
    
    label_group : dict
        A 2-level nested dict that after the first level of a pair tuple for purpose: "adjacency", "dictionary"
        and "vertex","edge" specification holds the inverse map of labels.

    laplacian_graph : np.array
        Holds the graph laplacian.
    
    lovasz_theta : float
        Holds the lovasz theta of the current graph.
        
    svm_theta : float
        Holds the svm theta of the current graph.
        
    metric_subgraphs_dict : dict
        Stores a calculation for subgraphs for each tuple consistent of:
        ("lovasz"/"svm", the number of samples, minimum sample size, maximum sample size)
        to avoid recalculation.
        
    _format: str, valid_values={"adjacency", "dictionary", "all"}
        Private attribute that keeps the current format.
    """
    
    # Adjacency Preset
    n = -1
    adjacency_matrix = np.empty(0)
    index_node_labels = dict()
    index_edge_labels = dict()
    
    # Dictionary Preset
    vertices = set()
    edge_dictionary = dict()
    node_labels = dict()
    edge_labels = dict()
    edsamic = dict()

    # Other Preset
    shortest_path_mat = None
    label_group = dict()
    laplacian_graph = None
    lovasz_theta =  None
    svm_theta = None
    metric_subgraphs_dict = dict()
    
    def __init__(self, initialization_object=None, node_labels=None, edge_labels=None, graph_format="auto"):
        if graph_format in ["adjacency","dictionary","auto","all"]:
            self._format = graph_format
            if (initialization_object is not None):
                self.build_graph(initialization_object, node_labels, edge_labels)
            elif graph_format is "auto":
                raise ValueError('no initialization object - format must not be auto')
        else:
            raise ValueError('Invalid graph format.\nValid graph formats are "all", "dictionary", "adjacency", "auto"')
            
    def build_graph(self, g, node_labels=None, edge_labels=None):
        """ Builds a graph structure given a supported graph representation.
        
        Parameters
        ----------
        g : a valid graph format
            Similar to intialization object (of ``__init__``).
            
        node_labels: dict, default=None
            Node labels dictionary relevant to g format.
            
        edge_labels: dict, default=None
            Edge labels dictionary relevant to g format.
            
        Returns
        -------
        None.
        """

        self.shortest_path_mat = None
        self.label_group = None
        case = 0
        if g is not None:
            if is_adjacency(g):
                # Input is considered an adjacency matrix
                case = 1
                
                # Assign labels for nodes
                self.index_node_labels = node_labels

                # Assign labels for edges
                self.index_edge_labels = edge_labels
                
                if(self._format is "auto"):
                    self._format = "adjacency"
            elif is_edge_dictionary(g):
                # Input is considered as an edge dictionary
                case = 2
                
                # Assign labels for nodes
                self.node_labels = node_labels
                
                # Assign labels for edges
                self.edge_labels = edge_labels
                
                if(self._format is "auto"):
                    self._format = "dictionary"
            else:
                raise ValueError('Unsupported input type. For more information check the documentation.')

        
        # If graph is of one type prune the other
        if self._format is "adjacency":
            edge_dictionary = None

        elif self._format is "dictionary":
            adjacency_matrix = None
        
        # Import the given input properly
        if (case==1):
            self._import_adjacency(g)
        elif (case==2):
            self._import_dictionary(g)

    def change_format(self, graph_format):
        """ Changes the format of the graph from an existing to an other.

        Paramters
        ---------
        graph_format : str, valid_values={"dictionary", "adjacency", "all"}
            Defines the internal representation of the graph object which can be a dictionary as a matrix, or both:
              + for dictionary: "dictionary"
              + for adjacency_matrix: "adjacency"
              + for both: "all"
        
        Returns
        -------
        None.
        """
        if graph_format not in ["all", "dictionary", "adjacency"]:
            raise ValueError('Invalid graph format for function change format. Valid formats are "all", "dictionary", "adjacency"')
        else:
            if (graph_format is not self._format):
                past_format = self._format
                self._format = graph_format
                if past_format is "adjacency":
                    self._import_adjacency()
                elif past_format is "dictionary":
                    self._import_dictionary()
                else:
                    if self._format is "dictionary":
                        self.n = -1
                        self.adjacency_matrix = None
                        self.edsamic = None
                        self.index_node_labels = None
                        self.index_edge_labels = None
                    else:
                        self.vertices = None
                        self.node_labels = None
                        self.edge_labels = None
                        self.edge_dictionary = None

    def desired_format(self, graph_format, warn=False):
        """ Changes the format to include the desired.

        Parameters
        ----------
        graph_format : str, valid_values={"dictionary", "adjacency", "all"}
            Defines the internal representation of the graph object which can be a dictionary as a matrix, or both:
              + for dictionary: "dictionary"
              + for adjacency_matrix: "adjacency"
              + for both: "all"
        
        warn : bool, default=False
            Warn the user if the format of the graph is being changed.
        
        Returns
        -------
        None.
        """
        if graph_format is "all":
            self.change_format(graph_format)
        elif graph_format is "dictionary":
            if self._format not in ["all","dictionary"]:
                if warn:
                    warnings.warn('changing format from "dictionary" to "all"')
                self.change_format("all")
        elif graph_format is "adjacency":
            if self._format not in ["all","adjacency"]:
                warnings.warn('changing format from "adjacency" to "all"')
                self.change_format("all")
    
    def construct_labels(self, label_type="vertex", purpose="adjacency"):
        """ A method that constructs labels (if user does not provide).
            
        Parameters
        ----------
        label_type: str, valid_values={"vertex", "edge"}, default="vertex"
            What kind of labels are going to be constructed
            
        purpose: str, valid_values={"adjacency", "dictionary"}, default="adjacency"
            Defines if the labels correspond to dictionary or adjacency matrix
            
        Returns
        -------
        None.
        """
        if (purpose == "adjacency"):
            if (label_type == "vertex"):
                nodes = list(range(0,self.n))
                self.index_node_labels = dict(zip(nodes,nodes))
            elif (label_type=="edge"):
                nodes = list(range(0,self.n))
                self.index_edge_labels = {(i,j):self.adjacency_matrix[i,j] for i in nodes for j in nodes}
            else:
                warnings.warn('unsupported label type')
        elif (purpose == "dictionary"):
            if (label_type == "vertex"):
                nodes = sorted(list(self.vertices))
                self.node_labels = dict(zip(nodes,nodes))
            elif (label_type=="edge"):
                self.index_edge_labels = {(i,j):(0 if (i not in self.edge_dictionary) or (j not in self.edge_dictionary[i]) else self.edge_dictionary[i][j]) for i in nodes for j in nodes}
            else:
                warnings.warn('unsupported label type')
                                
    def convert_labels(self, target_format="dictionary", purpose="all", init=False):
        """ A method that converts labels to a desired format. 

        Parameters
        ----------
        target_format: str, valid_values={"adjacency", "dictionary"}, default="dictionary"
            Defines what the target format for conversion will be.
            
        purpose: str, valid_values={"adjacency", "dictionary", "all"}, default="all"
            Defines if the labels will be converted for dictionary, adjacency matrix or both.
       
        init: bool, default=False
            An override parameter for format checks, usefull for initialisation.
        
        Returns
        -------
        None.
        """

        if (target_format == "adjacency"):
            if (self._format != "adjacency" or init):
                if bool(self.index_node_labels) and purpose in ['all', 'vertex']:
                    warnings.warn('overriding existing node labels for indexes')             
                if bool(self.index_edge_labels) and purpose in ['all', 'edges']:
                    warnings.warn('overriding existing edge labels for indexes')
                cond_labels_nodes = purpose in ['all', 'vertex'] and bool(self.node_labels)
                cond_labels_edges = purpose in ['all', 'edges'] and bool(self.edge_labels)
                if cond_labels_nodes or cond_labels_edges:
                    lov_sorted = sorted(list(self.vertices))
                    lv = len(lov_sorted)
                    if cond_labels_nodes:
                        self.index_node_labels = {i: self.node_labels[k] for (i,k) in zip(list(range(0,self.n)),lov_sorted)}
                    if cond_labels_edges:
                        self.index_edge_labels = {(i,j): self.edge_labels[(k,q)] for (i,k) in zip(list(range(0,lv)),lov_sorted) for (j,q) in zip(list(range(0,lv)),lov_sorted)}
                else:
                    if not init:
                        warnings.warn('no labels to convert from, for the given purpose')
            else:
                warnings.warn('labels already defined for that format - nothing to convert')
        elif (target_format == "dictionary"):
            if (self._format != "dictionary" or init):
                if bool(self.node_labels) and purpose in ['all', 'vertex']:
                    warnings.warn('overriding existing node labels, for dictionary symbols')
                    
                if bool(self.edge_labels) and purpose in ['all', 'edges']:
                    warnings.warn('overriding existing edge labels, for dictionary symbols')
                
                cond_labels_nodes = purpose in ['all', 'vertex'] and bool(self.index_node_labels)
                cond_labels_edges = purpose in ['all', 'edges'] and bool(self.index_edge_labels)
                if cond_labels_nodes or cond_labels_edges:
                    if cond_labels_nodes:
                        self.node_labels = self.index_node_labels
                    if cond_labels_nodes:
                        self.edge_labels = self.index_edge_labels
                else:
                    if not init:
                        warnings.warn('no labels to convert from, for the given purpose')
            else:
                warnings.warn('labels already defined for the given format')
    
    def label(self, obj, label_type="vertex", purpose="dictionary"):
        """ Returns the label of a vertex.
        
        Parameters
        ----------
        obj : hashable
            The candidate labeled object for the corresponding purpose and label type
            
        label_type : str, valid_values={"vertex", "edge"}, default="vertex"
            Defines if the lookup of the label will be done for vertices or edges.
        
        purpose : str, valid_values={"dictionary", "adjacency", "any"}, default="dictionary"
            Defines if the lookup will be done on the existing ("any" - if "all" default is adjacency)
            to the "dictionary" or to the "adjacency" format of the graph.
            
        Returns
        -------
        label : str, list, ..., *valid-label-type*
            Returns the label of the current object on the defined lookup.
            
        """
        if purpose == "any":
            if self._format == "all":
                purpose = "adjacency"
            else:
                purpose = self._format
                
        if purpose not in ["dictionary","adjacency"]:
            raise ValueError('unrecognized purpose')
        
        if (label_type == "vertex"):
            vertex = obj
            labels = self.get_labels(purpose=purpose, label_type=label_type)
            
            if labels is not None:
                if (vertex in labels):
                    return labels[vertex]
                else:
                    raise ValueError('no label assigned to this vertex symbol')
            else:
                warnings.warn('labels are not set - default value returned')
                return vertex
        elif (label_type == "edge"):
            edge = obj
            labels = self.get_labels(purpose=purpose, label_type=label_type)
            
            if not bool(labels):
                if (edge in labels):
                    return labels[edge]
                else:
                    raise ValueError('no label assigned to this edge symbol')
            else:
                warnings.warn('labels are not set - default value returned')
                if (edge[0] not in self.edge_dictionary) or (edge[1] not in self.edge_dictionary[edge[0]]):
                    return .0
                else:
                    return self.edge_dictionary[edge[0]][edge[1]]
    
    def relabel(self, new_labels, purpose="dictionary", label_type="vertex"):
        """ A method relabels the graph object, supporting the current format.
        
        Parameters
        ----------    
        new_labels : dict
            The new labels corresponding to the label type and purpose.
        
        purpose : str, valid_values={"dictionary", "adjacency"}, default="dictionary"
            Defines if the new labels are given for "adjacency" or "dictionary".
        
        label_type : str, valid_values={"vertex", "edge"}, default="vertex"
            Defines if the new labels are for vertices or edges.
               
        Returns
        -------
        None.
        """
        
        # checks
        if purpose not in ["dictionary","adjacency"]:
            raise ValueError('purpose can either be dictionary or adjacency')
            
        if label_type not in ["vertex", "edge"]:
            raise ValueError('label_type can either be vertex or edge')
        
        if new_labels is None or not bool(new_labels):
            warnings.warn('user must provide new labels, input is None')
        elif label_type=="vertex":
            if purpose is "dictionary":
                if self._format in ["dictionary","all"]:
                    self.node_labels = new_labels
                else:
                    raise ValueError('graph is in a format that supports labels for dictionary symbols')
                
                # Transform to a valid format
                if self._format is "all":
                    self.convert_labels("adjacency","vertex")
                if bool(self.label_group) and ((label_type,purpose) in self.label_group):
                    self.label_group[(label_type,purpose)] = dict()  
            elif purpose is "adjacency":
                if self._format in ["all","adjacency"]:
                    self.index_node_labels = new_labels
                else:
                    raise ValueError('graph is in a format that supports labels for indexes')
                
                if self._format is "all":
                    self.convert_labels("dictionary","vertex")
                if bool(self.label_group) and (label_type,purpose) in self.label_group:
                    self.label_group[(label_type,purpose)] = dict()
            else:
                warnings.warn('no labels to convert from, for the given purpose')
        elif label_type=="edge":
            if purpose is "dictionary":
                if self._format in ["dictionary","all"]:
                    self.edge_labels = new_labels
                else:
                    raise ValueError('graph is in a format that supports labels for dictionary symbols')
                
                # Transform to a valid format
                if self._format is "all":
                    self.convert_labels("adjacency","edge")
                if (label_type,purpose) in self.label_group:
                    self.label_group[(label_type,purpose)] = dict()  
 
            elif purpose is "adjacency":
                if self._format in ["all","adjacency"]:
                    self.index_edge_labels = new_labels
                else:
                    raise ValueError('graph is in a format that supports labels for indexes')

                if self._format is "all":
                    self.convert_labels("dictionary","edge")
                if (label_type,purpose) in self.label_group:
                    self.label_group[(label_type,purpose)] = dict()  

            else:
                warnings.warn('no labels to convert from, for the given purpose')       
        else:
            warnings.warn('unrecognized label type')

        
    def build_shortest_path_matrix(self, algorithm_type="auto", clean=False, labels="vertex"):
        """ A method that builds and returns the shortest path matrix between all nodes.

        Parameters
        ----------
        algorithm_type : str, valid_values={"auto", "adjacency", "dictionary"}, default="auto"
            Defines which shortest-path algorithm will be used for building the shortest path matrix.
             + "dijkstra" : choses the dijkstra algorithm (Matrix computation complexity: :math:`O(|V|(|E|+|V|) \log{|V|})`
             + "floyd_warshall" : choses the floyd-warshall algorithm (Matrix computation complexity: :math:`O(|V|^3)`)
             + "auto" : choses the best possible algorithm for the current format
        
        clean : bool, default=False
            Construct the shortest path matrix from scratch or output existing if exists
            
        labels : str, valid_values={"vertex", "edge", "all"}, default="vertex" 
            Returns labels corresponding for the indexes of the shortest path matrix
            for vertices, for edge (only for the valid ones on the original graph) or for both ("all")
            
        Returns
        -------
        shortest_path_matrix : np.array, shape=(:math:`|V|`,:math:`|V|`)
            The produced shortest path matrix.
            
        vertex_labels : dict
            The vertex labels, outputed only if *labels* parameter is either "vertex" or "all".
            
        edge_labels : dict
            The edge labels, outputed only if *labels* parameter is either "edge" or "all".
        """
        if labels not in ['vertex', 'edge', 'all']:
            raise ValueError('only labels for vertices and edges exist')
        
        if clean:
            self.shortest_path_mat = None
            
        if self.shortest_path_mat is not None:
            return self.shortest_path_mat, self.index_node_labels

        # Assign the desired algorithm
        if algorithm_type is "auto":
            if self._format in ["all","dictionary"]:
                algorithm_type = "dijkstra"
            elif self._format is "adjacency":
                algorithm_type = "floyd_warshall"
        
        if algorithm_type is "dijkstra":
            # all format required
            self.desired_format("all", warn=True)
            if (bool(self.edsamic)):
                indexes = self.edsamic
            else:
                indexes = dict(zip(list(range(self.n)),list(range(self.n))))
            
            # calculate shortest path matrix
            shortest_path_mat = np.full([self.n, self.n],float("Inf"))
            for k in indexes.keys():
                dict_fd, _ = dijkstra(self.edge_dictionary,k)
                for s in dict_fd.keys():
                    shortest_path_mat[indexes[k],indexes[s]] = dict_fd[s]

        elif algorithm_type is "floyd_warshall":
            self.desired_format("adjacency", warn=True)
            shortest_path_mat = floyd_warshall(self.adjacency_matrix)
        
        self.shortest_path_mat = shortest_path_mat
        if labels is "all":
            return shortest_path_mat, self.get_labels(), self.get_labels("edge")
        if labels is "edge":
            return shortest_path_mat, self.get_labels("edge")
        if labels is "vertex":
            return shortest_path_mat, self.get_labels()
            
    def get_labels(self, label_type="vertex", purpose="adjacency"):
        """ Return labels corresponding to the purpose.
        
        Parameters
        ----------
        label_type : str, valid_values={"vertex", "edge"}, default="vertex"
            Defines if the the labels will correspond to vertices or edges.
        
        purpose : str, valid_values={"dictionary", "adjacency", "any"}, default="dictionary"
            Defines if the labels will correspond to "dictionary" (symbols), to "adjacency"
            indexes, or to "any" valid format (if "all" the result is for "adjacency")
        
        Returns
        -------
        labels : dict,
            Returns the labels for the given type and purpose.
        """
        if (purpose == "adjacency"):
            case = True
        elif (purpose == "dictionary"):
            case = False
        elif (purpose == "any"):
            if self._format in ['all', 'adjacency']:
                case = True
            else:
                case = False
        else:
            raise ValueError('unsupported label purpose')
            
        if case:
            self.desired_format("adjacency", warn=True)
            if label_type == "vertex":
                if not bool(self.index_node_labels):
                    self.construct_labels(label_type, purpose)
                return self.index_node_labels
            elif label_type == "edge":
                if not bool(self.index_edge_labels):
                    self.construct_labels(label_type, purpose)
                return self.index_edge_labels
            else:
                raise ValueError('label type can only be "vertex" or "edge"')
        else:
            self.desired_format("dictionary", warn=True)
            if label_type == "vertex":
                if not bool(self.node_labels):
                    self.construct_labels(label_type, purpose)
                return self.node_labels
            elif label_type == "edge":
                if not bool(self.edge_labels):
                    self.construct_labels(label_type, purpose)
                return self.edge_labels
            else:
                raise ValueError('label type can only be "vertex" or "edge"')
    
    def get_label_group(self, label_type="vertex", purpose="dictionary"):
        """ A function that calculates the inverse dictionary for vertex labels (once).
            
        Parameters
        ----------
        label_type : str, valid_values={"vertex", "edge"}, default="vertex"
            Defines if the the labels-group will correspond to vertices or edges.
        
        purpose : str, valid_values={"dictionary", "adjacency", "any"}, default="dictionary"
            Defines if the labels-group will correspond to "dictionary" (symbols), to "adjacency" (indexes)
            or to "any" valid format (if "all" the result is for "adjacency")
        
        Returns
        -------
        label_group : dict
            Returns the inverse label group.
        """
        if not bool(self.label_group):
           self.label_group = dict()
        if not (label_type,purpose) in self.label_group:
           self.label_group[(label_type,purpose)] = dict()
        if not bool(self.label_group[(label_type,purpose)]):
            # calculate the label_group
            self.label_group[(label_type,purpose)] = inv_dict(self.get_labels(label_type,purpose))
        return self.label_group[(label_type,purpose)]
                
    def neighbors(self, vertex, purpose="any", with_weights = False):
        """ Find all neighbors of a vertex.
           
        Parameters
        ----------
        vertex : hashable
            The vertex, which neighbors we are searching for.
        
        purpose : str, valid_values={"adjacency", "dictionary", "any"}, default="any"
            Defines if the vertex is given for the "dictionary" format of the graph (symbol)
            to the "adjacency" (index) or to "any" existing format (if "all" the expected type is for "adjacency").
        
        with_weights : bool, default=False
            Defines if the neighbours will be outputed with weights.
                
        Returns
        -------
        neighbors : list or dict
            The neighbors of the given vertex. 
                + with_weights=False: list of neighbor vertices
                + with_weights=True: dictionary between neighbor vertices and edge labels 
        """
        if purpose in ['adjacency', 'dictionary', "any"]:
            if purpose == 'dictionary':
                self.desired_format('dictionary')
                case = True
            if purpose == 'adjacency':
                self.desired_format('adjacency')
                case = False
            if purpose == "any":
                if self._format in ['all','adjacency']:
                    case = False
                else:
                    case = True
        else:
            raise ValueError('purpose is either "adjacency", "dictionary" or "any"')
                
        if case:
            if vertex in self.edge_dictionary:
                if not with_weights:
                    return list(self.edge_dictionary[vertex].keys())
                else:
                    return self.edge_dictionary[vertex]
            else:
                raise ValueError('vertex not inside edge dictionary')
        else:
            idx = int(vertex)
            if 0 <= idx < self.n:
                out_idx = np.where(self.adjacency_matrix[idx,:] > 0)
                ns = out_idx[0].tolist()
                if not with_weights:
                    return ns
                else:
                    return dict(zip(ns,self.adjacency[idx,out_idx].tolist())[0])
            else:
                raise ValueError("item with index "+str(idx)+" does not exist")

    def _make_edsamic(self, vertices):
        """ A method to produce the *edge-symbols-adjacency-matrix-correspondance-dictionary*.
        
        Parameters
        ----------
        vertices : list-able, sort-able
            The vertices of the graph.
        
        all_out : bool, default=True
            A variable that defines, if everything valuable will be outputed.
            
        Returns
        -------
        edsamic : dict
            The produced edsamic
        n : int
            The number of vertieces
        lov_sorted : list
            Sorted list of vertices.
            
        """
        if bool(vertices):
            # vertices to list
            lov = list(vertices)

            # length is adjacency matrix size
            n = len(lov)

            # sort lov indexes for labeling
            lov_sorted = sorted(lov)
            
            # edge_labels_adjacency_matrix_correspondance_dictionary            
            edsamic = dict(zip(lov_sorted,list(range(0,n))))
            
            return edsamic, n, lov_sorted
        else:
            warnings.warn('wrong format: edge dictionary must exist')
            return None


    def _import_adjacency(self, adjacency_matrix=None, init=True):
        """ A function that creates a graph object representation given its adjacency matrix.
        
        Parameters
        ----------
        adjacency_matrix : array-like, default=None
            If given a array the input can be as follows:
            + array-like lists of lists
            + np.array
            + sparse matrix (scipy.sparse.csr.csr_matrix)
            If None, imports the existing array, inside self.adjacency_matrix
        
        init : bool, default=True
            A parameter used to defined initialization.
            
        Returns
        -------
        None.
        """
        if adjacency_matrix is not None:
            # calculate graph size
            is_adj, adjacency_matrix = is_adjacency(adjacency_matrix, True)
            
            if not is_adj:
                raise ValueError('unsupported format type for adjacency matrix')
                
            n = adjacency_matrix.shape[0]
            
            if n != adjacency_matrix.shape[1]:
                raise ValueError('input matrix must be squared')
            
            # import_adjacency
            if (self._format is "all" or self._format is "adjacency"):
                self.adjacency_matrix = adjacency_matrix
            self.n = n
        else:
            n = self.n
            adjacency_matrix = self.adjacency_matrix

        # construct a dictionary out of the adjacency
        if self._format in ["all", "dictionary"]:
            vertices = set(list(range(0, n)))
            
            self.vertices = vertices
            self.edge_dictionary = {i: dict() for i in range(0,n)}
            
            idx_i, idx_j = np.where(adjacency_matrix > 0)
            for (i,j) in zip(idx_i,idx_j):
                self.edge_dictionary[i][j] = adjacency_matrix[i,j]

            # Add labels
            self.convert_labels(target_format="dictionary", purpose="all", init=init)
            
            # Prune not interesting format
            if self._format is "dictionary":
                self.n = -1
                self.adjacency_matrix = None
                self.edsamic = None
                self.index_node_labels = None
                self.index_edge_labels = None
            else:
                self.edsamic = dict(zip(range(0,n),range(0,n)))
            
    def _import_dictionary(self, edge_dictionary=None, init=True):
        """ A function that creates a graph object representation given its edge dictionary.

        Parameters
        ----------
        edge_dictionary :  dict-like, default=None
            The edge_dictionary the input can be as follows:
                + 2-level nested dictionaries from edge symbols to weights
                + Dictionary of symbols to list of symbols (unweighted)
                + Dictionary of tuples to weights (weighted)
                + Iterable of tuples of len 2 (unweighted)
                + Iterable of tuples of len 3 (vertex, vertex, weight)
           If None, imports from the existing edge_dictionary
            
        init : bool, default=True
            A parameter used to defined initialization.
            
        Returns
        -------
        None.
        """
        if edge_dictionary is not None:
            # find vertices, refine dictionary
            vertices = set()
            edge_dictionary_refined = dict()
            is_edge_dict, vertices, edge_dictionary_refined = is_edge_dictionary(edge_dictionary, True)
            if not is_edge_dict:
                raise ValueError('unsupported edge_dictionary format')
            
            # Save dictionary, vertices @self if needed        
            if (self._format == "all" or self._format == "dictionary"):
                self.edge_dictionary = edge_dictionary_refined
            self.vertices = vertices
        else:
            vertices = self.vertices
            edge_dictionary_refined = self.edge_dictionary

        # Create and store the adjacency matrix
        if self._format in ["adjacency","all"]:
            edsamic, self.n, lov_sorted = self._make_edsamic(vertices)
            # Initialize adjacency_matrix
            adjacency_matrix = np.zeros(shape = (self.n,self.n))
            
            # Produce and save adjacency matrix
            for k in edge_dictionary_refined.keys():
                for l in edge_dictionary_refined[k].keys():
                    adjacency_matrix[edsamic[k],edsamic[l]] = edge_dictionary_refined[k][l]
                
            self.adjacency_matrix = adjacency_matrix
            
            if self._format == "all":
                self.edsamic = edsamic
            
            # Add labels
            self.convert_labels(target_format="adjacency", purpose="all", init=init)
            
            # Prune for a certain format
            if self._format == "adjacency":
                self.edge_dictionary = None
                self.vertices = None
                self.node_labels = None
                self.edge_labels = None


    def laplacian(self, save=True):
        """ Calculates the laplacian of the given graph.
        
        Parameters
        ----------
        save: bool, default=True
            Optional parameter to store the matrix.
            
        Returns
        -------
        laplacian: array-like
            Returns the graph laplacian
        """ 
        if self.laplacian_graph is not None:
            laplacian_graph = self.laplacian_graph
        else:
            self.desired_format("adjacency", warn=True)
            laplacian_graph = laplacian(self.adjacency_matrix)
            
            if save:
                self.laplacian_graph = laplacian_graph
        return laplacian_graph
    
    def calculate_lovasz_theta(self, from_scratch=False):
        """ Calculates the lovasz theta for the given graph.
            
        Parameters
        ----------
        from_scratch: bool, default=False
            Defines if the the lovasz_theta will be calculated from scratch.
        
        Returns
        -------
        lovasz_theta: float
            Returns the lovasz theta number.
        """
        if not from_scratch:
            if self.lovasz_theta is not None:
                return self.lovasz_theta
        
        # Calculate dependent on the format
        # sparse matrix and required graph 
        # parameters for defining the convex 
        # optimization problem
        if self._format == "dictionary":
            nv = len(vertices)
            if nv == 1:
                return 1.0
                
            else:
                enum_vertices = {i: v for (i,v) in enumerate(sorted(list(vertices)))}
                
                # Calculate the sparse matrix needed to feed
                # into the convex optimization problem
                q = {(v,n) for v in list(vertices) for n in self.edge_dictionary[v].keys() if v!=n}
                
                x_list = list()
                e_list = list()
                ei = 0
                for i in range(0, nv):
                    na = enum_vertices[i]
                    for j in range(i+1,nv):
                        nb = enum_vertices[j]
                        if (na,nb) in q:
                            e_list.append(int(ei)), x_list.append(int(i*nv+j))
                            if (nb,na) in q:
                                e_list.append(int(ei)), x_list.append(int(i*nv+j))
                            ei+=1
                ne = ei-1
                
        elif self._format in ["adjacency","all"]:
            if self.n == 1:
                return 1.0
            else:
                tf_adjm = self.adjacency_matrix>0
                i_list, j_list = np.nonzero(np.triu(tf_adjm.astype(int), k=1))
                
                nv = self.n
                ne = len(i_list)
                
                x_list = list()
                e_list = list()
                for (e,(i,j)) in enumerate(zip(i_list,j_list)):
                    e_list.append(int(e)), x_list.append(int(i*nv+j))
                    if tf_adjm[i,j]:
                        e_list.append(int(e)), x_list.append(int(i*nv+j))
        
        # Add on the last row, diagonal elements
        e_list = e_list+(nv*[ne])
        x_list = x_list+[int(i*nv + i) for i in range(0,nv)]
        
        # initialise g sparse (to values -1, based on two list that 
        # define index and one that defines shape
        g_sparse = cvxopt.spmatrix(-1, x_list, e_list, (nv*nv,ne+1))
        
        # Initialise optimization parameters
        h = cvxopt.matrix(-1.0, (nv,nv))
        c = cvxopt.matrix([0.0]*ne + [1.0])
                
        # Solve the convex optimization problem
        sol = cvxopt.solvers.sdp(c, Gs=[g_sparse], hs=[h])
        self.lovasz_theta = sol['x'][ne]
        return self.lovasz_theta
    
    def calculate_svm_theta(self, from_scratch=False):
        """ Calculates the svm theta for the given graph.
            
        Parameters
        ----------
        from_scratch: bool, default=False
            Defines if the the svm_theta number will be calculated from scratch.
        
        Returns
        -------
        svm_theta: float
            Returns the svm theta number.
        """
        if not from_scratch:
            if self.svm_theta is not None:
                return self.svm_theta
        
        # Calculate dependent on the format
        # sparse matrix and required graph 
        # parameters for defining the convex 
        # optimization problem
        if self._format == "dictionary":
            nv = len(vertices)
            lov_sorted = sorted(list(vertices))
            vertices = {v: i for (i,v) in enumerate(lov_sorted)}
            
            K = np.zeros(shape=(nv,nv))
            idx = zip(*[(vertices[v],vertices[n]) for v in lov_sorted for n in self.neighbors(v)])
            idx1 = list(idx[0])
            idx2 = list(idx[1])
            
            K[idx1,idx2] = 1
            np.fill_diagonal(K, 1)
                        
        elif self._format in ["adjacency","all"]:
            K = self.adjacency_matrix>0
            np.fill_diagonal(K, False)
            K.astype(int)
        
        svm = OneClassSVM(kernel="precomputed")
        svm.fit(K)
                
        self.svm_theta = np.sum(np.abs(svm.dual_coef_[0]),axis=0)
        return self.svm_theta        
    
    def calculate_subgraph_samples_metric_dictionary(self, metric_type, n_samples=50, subsets_size_range=(2,8), from_scratch=False, save=True):
        """ A function useful for calculating a graph metric as the lovasz theta kernel or the svm-theta kernel, on a set of randomly sampled subgraphs, producing a dictionary of subgraph levels and sets.
            
        Parameters
        ----------
        n_samples : int, default=50
            The number of samples that will be sampled.
            
        subsets_size_range : tuple, len=2, default=(2,8)
            A tuple having the min and the max subset size.
        
        from_scratch : bool
            Defines if the the metric will be calculated from scratch.

        save : bool
            Determine if the metric output will be stored.
            
        Returns
        -------
        level_values : dict
            Returns a dictionary with levels (subsite size) and the lovasz value of all the sampled subgraphs.
        """
        if metric_type == "svm":
            metric = lambda G: G.calculate_svm_theta()
        elif metric_type == "lovasz":
            metric = lambda G: G.calculate_lovasz_theta()            
        else:
            raise ValueError('unsupported metric type')
            
        min_ss, max_ss = subsets_size_range[0], subsets_size_range[1]
        # If precalculated and the user wants it, output
        if not from_scratch:
            if bool(self.metric_subgraphs_dict):
                if (metric_type, n_samples, min_ss, max_ss) in self.metric_subgraphs_dict:
                    return self.metric_subgraphs_dict[(metric_type, n_samples, min_ss, max_ss)]
        
        
        self.desired_format("adjacency", warn=True)
        
        ## Calculate subsets
        samples_on_subsets = distribute_samples(self.n, subsets_size_range, n_samples)
        
        # Calculate level dictionary with lovasz values
        level_values=dict()
        for level in samples_on_subsets.keys():
            subsets = set()
            for k in range(0,samples_on_subsets[level]):
                while True:
                    indexes = np.random.choice(self.n, level, replace=False)
                    if tuple(indexes) not in subsets:
                        subsets.add(tuple(indexes))
                        break
                # calculate the lovasz value for that level
                if level not in level_values:
                    level_values[level] = list()
                level_values[level].append(metric(graph(self.adjacency_matrix[:,indexes][indexes,:])))
                
        if save:
            self.metric_subgraphs_dict[(metric_type, n_samples, min_ss, max_ss)] = level_values
        
        return level_values
    
    
    def get_vertices(self, purpose="adjacency"):
        """ A method that returns an iterable of vertices.
        
        Parameters
        ----------
        purpose : str, valid_values={"adjacency", "dictionary", "any"}, default="adjacency"
            Defines if the vertices will correspond given for the "dictionary" format of the graph (symbol)
            to the "adjacency" (index) or to "any" existing format (if "all" the expected type is for "adjacency")
        
        Returns
        -------
        vertices : iterable
            Returns an iterable on vertices
        """
        if purpose not in ["adjacency", "dictionary", "any"]:
            raise ValueError('purpose is either "adjacency" of "dictionary"')

        if purpose == "any":
            if self._format in ['all','adjacency']:
                purpose = "adjacency"
            else:
                purpose = "dictionary"
                
        if purpose == "adjacency":
            self.desired_format("adjacency", warn=True)
            return range(0,self.n)
        if purpose == "dictionary":
            self.desired_format("dictionary", warn=True)
            return self.vertices
        
    def get_edges(self, purpose="adjacency", with_weights=False):
        """ A method that returns an iterable of edges as tuples.
        
        Parameters
        ----------
        purpose : str, valid_values={"adjacency", "dictionary"}, default="adjacency"
            Defines if the edges is given for the "dictionary" format of the graph (symbol) to the "adjacency" (index).
        
        Returns
        -------
        vertices : list
            Returns a list of tuples for edges.
        """
        if purpose not in ["adjacency", "dictionary"]:
            raise ValueError('purpose is either "adjacency" of "dictionary"')
        
        if purpose == "adjacency":
            self.desired_format("adjacency", warn=True)
            idx_i, idx_j = np.where(self.adjacency_matrix > 0)
            edges = zip(idx_i, idx_j)
            if with_weights:
                return list(zip(edges,self.adjacency_matrix[idx_i, idx_j]))
            else:
                return list(edges)
        if purpose == "dictionary":
            self.desired_format("dictionary", warn=True)
            if with_weights:
                return [(i,j) for i in self.edge_dictionary.keys() for j in self.edge_dictionary[i].keys()]
            else:
                return [((i,j), self.edge_dictionary[i][j]) for i in self.edge_dictionary.keys() for j in self.edge_dictionary[i].keys()]
    
    def get_adjacency_matrix(self):
        """ A format agnostic method that returns the adjacency matrix.
        
        Parameters
        ----------
        None.
        
        Returns
        -------
        adjacency_matrix : np.array
            Returns the adjacency matrix of the current graph.
        """
        if self._format is "dictionary":
            A = np.zeros(shape=(len(self.vertices),len(self.vertices)))
            v_map = {v:i  for (v,i) in enumerate(sorted(list(self.vertices)))}
            for (k,v) in self.edge_dictionary.items():
                i = v_map[k]
                for (kv,w) in v.items():
                    A[i,v_map[kv]] = w
            return A
        else:
            return self.adjacency_matrix
            
    def get_edge_dictionary(self):
        """ A format agnostic method that returns the edge dictionary
        
        Parameters
        ----------
        None.
        
        Returns
        -------
        edge_dictionary : dict
            Returns the edge_dictionary of the current graph.
        """
        if self._format is "dictionary":
            idx_i, idx_j = np.where(self.adjacency_matrix > 0)
            edge_dictionary = {i:dict() for i in range(0, self.n)}
            for (i,j) in zip(idx_i,idx_j):
                edge_dictionary[i][j] = self.adjacency_matrix[i,j]
            return edge_dictionary
        else:
            return self.edge_dictionary
            
    
    def nv(self):
        """ A method that returns the number of vertices for any existing format.
        
        Parameters
        ----------
        None.
        
        Returns
        -------
        num_of_vertices : int
            Returns the number of vertices.
        """
        if self._format in ['all','adjacency']:
            return self.n
        else:
            return len(self.vertices)
            
    def produce_neighborhoods(self, r=3, purpose="adjacency", with_distances=False, d=-1):
        """ Calculates neighborhoods for each node of a Graph up to a depth.
            
        Parameters
        ----------
        r : int
            The neighborhood depth (radius).
            
        purpose : str, valid_values={"adjacency", "dictionary"}, default="adjacency"
            Defines if the neighborhood symbols will correspond given for the "dictionary" 
            format of the graph (symbol) to the "adjacency" (index).
            
        with_distances: bool, default=False
            Defines if we need to calculate BFS distances for each pair.
            
        d : int, default=-1
            Maximum distance considered. If -1 is provided the distance is as max as the radius r.
            
        Returns
        -------
        N : dict
            A level, vertex nested dictionary of lists, corresponding to the neighbours of level
            :math:`l` for a certain vertex :math:`v`.
        
        D : dict
            For each level, set of tuples of nodes connected in that level.
            Appears only if with_distances is *True*.
            
        Dist_pair : dict
            A dictionary of all pairs and their distances.
        """
        N = dict()

        if purpose not in ["adjacency", "dictionary"]:
            raise ValueError('purpose is either "adjacency" or "dictionary"')
            
        if r < 0:
            raise ValueError('r must be positive or equal to zero')

        if with_distances and d<0:
            d = r
            warnings.warn('negative d as input - d set to r')
            
        # initialization
        N[0] = dict(zip(self.get_vertices(purpose),self.get_vertices(purpose)))
        
        if with_distances:
            D = dict()
            D[0] = set(zip(self.get_vertices(purpose),self.get_vertices(purpose)))
            Dist_pair = {(v,v): 0 for v in self.get_vertices(purpose)}
            
        if r > 0:
            N[1] = dict()
            if with_distances and d>=1:
                D[1] = set()
            for i in self.get_vertices(purpose):
                ns = list(self.neighbors(i,purpose))
                N[1][i] = sorted([i]+ns)
                if with_distances and d>=1:
                    dset = {(i,n) for n in ns}
                    Dist_pair.update(zip(dset,len(dset)*[1]))
                    D[1] |= dset 
            # calculate neighborhoods
            # by a recursive formula
            # for all levels from 1 to r
            for level in range(1,max(r,d)):
                N[level+1] = dict()
                if with_distances and level<=d-1:
                    D[level+1] = set()
                for i in self.get_vertices(purpose=purpose):
                    neighbours = set()
                    for w in N[level][i]:
                        neighbours |= set(N[level][w])
                    N[level+1][i] = sorted(list(neighbours))
                    if with_distances and level<=d-1:
                        dset = {(i,j) for j in (neighbours - set(N[level][i]))}
                        Dist_pair.update(zip(dset,len(dset)*[level+1]))
                        D[level+1] |= dset
            if with_distances:
                for level in range(r+1,d):
                    N.pop(level,None)
        
        if with_distances:     
            return N, D, Dist_pair
        else:
            return N
        
    def get_subgraph(self, vertices):
        """A method that creates a graph object subgraph in the same format as the original graph.
    
        Parameters
        ----------
        vertices: iterable
            An iterbale vertices extracted from the original graph.
        
        Returns
        -------
        subgraph : graph
            The induced subgraph, from the input vertices.
        """
        if vertices is None:
            raise ValueError('vertices must not be null')
        if type(vertices) in [str,int]:
            vertices = {vertices}
        else:
            vertices = set(vertices).copy()

        subgraph = graph(graph_format=self._format)        
        if self._format == 'adjacency':
            for v in vertices:
                if v<0 or v>=self.n:
                    raise ValueError('vertices are not valid for the original graph format')
            recipe = [tuple(['enum_vertices','default']), tuple(['add_adjacency'])]
            
        elif self._format == 'dictionary':
            for v in vertices:
                if v not in self.edge_dictionary:
                    raise ValueError('vertices are not valid for the original graph format')
            
            recipe = [tuple(['get_correct', 'default']), tuple(['add_adjacency'])]
            get_correct = lambda i: i
        else:
            fv = False
            fa = False
            for v in vertices:
                if not fv and v not in self.edge_dictionary:
                    fv = True        
                if not fa and (v<0 or v>=self.n):
                    fa = True
                if fa and fv:
                    raise ValueError('vertices are not valid for the original graph format')
    
            if not fa and fv:
                recipe = [tuple(['enum_vertices','default']), tuple(['get_correct', 'edsamic']), tuple(['add_adjacency']), tuple(['add_edge_dict'])]
            elif not fa:
                recipe = [tuple(['enum_vertices','default']), tuple(['get_correct', 'edsamic']), tuple(['add_adjacency']), tuple(['idxs_to_nodes']), tuple(['add_edge_dict'])]
            elif not fv:
                recipe = [tuple(['enum_vertices','edsamic']), tuple(['get_correct', 'edsamic']), tuple(['add_adjacency']), tuple(['add_edge_dict'])]
                
        for ingredient in recipe:
            if ingredient[0] == 'idxs_to_nodes':
                vertices = {self.edsamic[i] for i in vertices}
                
            elif ingredient[0] == 'enum_vertices':
                if ingredient[1] == 'edsamic':
                    inverse_edsamic = inv_dict(self.edsamic)
                    lov_sorted = sorted([int(inverse_edsamic[v][0]) for v in vertices])
                    new_indexes, inv_new_indexes = {}, {}
                    for (i,l) in enumerate(lov_sorted):
                        new_indexes[i], inv_new_indexes[l] = l, i
                    subgraph.edsamic = {self.edsamic[inv_new_indexes[nik]]: nik for nik in new_indexes.keys()}
                else:
                    lov_sorted = sorted(list(vertices))
                    new_indexes = {l:i for (i,l) in enumerate(lov_sorted)}
                
            elif ingredient[0] == 'add_adjacency':
                subgraph.adjacency_matrix = self.adjacency_matrix[lov_sorted,:][:,lov_sorted]
                subgraph.n = len(new_indexes.keys())
                if bool(self.index_node_labels):
                    subgraph.index_node_labels = {new_indexes[k]: self.index_node_labels[k] for k in self.index_node_labels.keys() if k in vertices}
                if bool(self.index_edge_labels):
                    subgraph.index_edge_labels = {(new_indexes[i],new_indexes[j]): self.index_edge_labels[(i,j)] for (i,j) in self.index_edge_labels if (i in vertices and j in vertices)}
                    
            elif ingredient[0] == 'add_edge_dict':
                subgraph.edge_dictionary = {get_correct(i):{get_correct(j): self.edge_dictionary[i][j] for j in self.edge_dictionary[i] if j in vertices} for i in self.edge_dictionary if i in vertices}
                subgraph.vertices = vertices
            
                if bool(self.node_labels):
                    subgraph.node_labels = {get_correct(k): self.node_labels[k] for k in self.node_labels.keys() if k in vertices}
                if bool(self.edge_labels):
                    subgraph.index_edge_labels = {(get_correct(i),get_correct(j)): self.edge_labels[(i,j)] for (i,j) in self.edge_labels if (i in vertices and j in vertices)}
            elif ingredient[0] == 'get_correct':
                if ingredient[1] == 'edsamic':
                    get_correct = lambda i: self.edsamic[new_indexes[i]]
                else:
                    get_correct = lambda i: i
            
        return subgraph


def is_adjacency(g, transform=False):
    """ A function that defines if input is in a valid adjacency matrix format.
    
    Parameters
    ----------
    g : Object
        The input object.
    transform : bool, default=False
        Defines if the input will be transformed to the internal adjacency matrix support format.
    
    Returns
    -------
    is_adjacency : bool
        A variable that determines if the input is a valid adjacency matrix.
    g_transformed : np.array
        Holds the transformed object to an np.array.
        This output appears **only** if transform parameter is True.
    """
    
    if type(g) in [np.array, np.ndarray] and len(g.shape)==2:
        if transform:
            return True, g
        else:
            return True
    elif type(g) is scipy.sparse.csr.csr_matrix:
        if transform:
            return True, g.todense()
        else:
            return True
    elif type(g) is list and all(isinstance(l, list) and all(isinstance(i, numbers.Number) for l in g for i in l)):
        if transform:
            return True, np.array(g)
        else:
            return True
    else:
        if transform:
            return False, None
        else:
            return False
                    
def is_edge_dictionary(g, transform=True):
    """ A function that defines if input is in a valid edge dictionary format.
    
    Parameters
    ----------
    g : Object
        The input object.
    transform : bool, default=False
        Defines if the input will be transformed to the internal edge dictionary support format.
    
    Returns
    -------
    is_edge_dictionary : bool
        A variable that determines if the input is a valid edge dictionary.
    g_transformed_vertices : set
        Holds the transformed object vertices of the edge_dictionary.
    g_transformed_edge_dict : dict
        Holds the transformed object as a 2-level edge dict.
        This output appears **only** if transform parameter is True.
    """
    if type(g) is dict:
        if all(type(k) is tuple and len(k)==2 and isinstance(n, numbers.Number) for (k,n) in g.items()):
            if transform:
                vertices_key = set()
                vertices_val = set()
                edge_dict = dict()
                for (k,v) in g.items():
                    vertices_key.add(k[0])
                    vertices_val.add(k[1])
                    nested_dict_add(edge_dict, v, k[0], k[1])
                
                # intialise empty edges
                v_dif = vertices_val - vertices_key
                if len(v_dif) > 0:
                    for v in v_dif:
                        edge_dict[v] = dict()
                return True, vertices_key | vertices_val, edge_dict
            else:
                return True
        if all(isinstance(d, list) for d in g.values()):
            if transform:
                vertices_key = set()
                vertices_val = set()
                edge_dict = dict()
                for (k,v) in g.items():
                    vertices_key.add(k)
                    vertices_val |= set(v)
                    for kp in v:
                        nested_dict_add(edge_dict, 1, k, kp)
                        
                # intialise empty edges
                v_dif = vertices_val - vertices_key
                if len(v_dif) > 0:
                    for v in v_dif:
                        edge_dict[v] = dict()
                return True, vertices_key | vertices_val, edge_dict
            else:
                return True
        if all(isinstance(d, dict) and all(isinstance(n, numbers.Number) for n in d.values()) for d in g.values()):
            if transform:
                vertices_key = set(g.keys())
                vertices_val = {kp for k in g.keys() for kp in g[k].keys()}
                v_dif = vertices_val - vertices_key
                
                # intialise empty edges
                if len(v_dif) > 0:
                    for v in v_dif:
                        g[v] = dict()
                return True, vertices_key | vertices_val, g
            else:
                return True
    if isinstance(g, collections.Iterable):
        if all(type(t) is tuple and len(t)==2 for t in g):
            if transform:
                vertices_key = set()
                vertices_val = set()
                edge_dict = dict()
                for (v, u) in g:
                    vertices_key.add(v)
                    vertices_val.add(u)
                    nested_dict_add(g_transformed, 1, v, u)
                v_dif = vertices_val - vertices_key
                
                # intialise empty edges
                if len(v_dif) > 0:
                    for v in v_dif:
                        edge_dict[v] = dict()
                return True, vertices_key | vertices_val, edge_dict
            else:
                return True
        elif all(type(t) is tuple and len(t)==3 for t in g):
            if transform:
                vertices_key = set()
                vertices_val = set()
                edge_dict = dict()
                for (v, u, w) in g:
                    vertices_key.add(v)
                    vertices_val.add(u)
                    nested_dict_add(g_transformed, w, v, u)
                v_dif = vertices_val - vertices_key
                
                # intialise empty edges
                if len(v_dif) > 0:
                    for v in v_dif:
                        edge_dict[v] = dict()
                return True, vertices_key | vertices_val, edge_dict
            else:
                return True
    if transform:
        return False, None
    else:
        return False
        
def laplacian(A):
    """ Calculates the laplacian given the adjacency matrix.
    
    Parameters
    ----------
    A : np.array, square
        The adjacency matrix of a graph.
    
    Returns
    -------
    laplacian_graph : np.array
        The produced laplacian graph.    
    """        
    laplacian_graph = np.zeros(shape=(n, n))
    for i in range(0, n):
        for j in range(0, n):
            if i!=j and A[i,j]!=.0:
                laplacian_graph[i,j] = -A[i,j]
                laplacian_graph[i,i]+= A[i,j]
    return laplacian_graph

def dijkstra(edge_dictionary, start_vertex, end_vertex=None):
    """ Implementation of the dijkstra algorithm[fn:1].
    
    Parameters
    ----------
    edge_dictionary: dict
        A 2-level nested dictionary of symbols, with value corresponding to the weight.

    start_vertex: hashable
        The start vertex symbol (should exists as a key inside edge_dictionary).
    
    end_vertex: hashable
        The end vertex symbol (should exists as a key inside edge_dictionary).
        
    Returns:
    --------
    dict_fd : dict
        The dictionary of final distances.
    dict_pred : dict
        The dictionary of predecessors.
        
    [fn:1] The majority of this function code came from: `Link text <http://code.activestate.com/recipes/119466-dijkstras-algorithm-for-shortest-paths/>`_
    """

    dict_fd = {}    # dictionary of final distances
    dict_pred = {}    # dictionary of predecessors
    queue = priority_dict()   # est.dist. of non-final vert.
    queue[start_vertex] = 0
    
    for v in queue:
        dict_fd[v] = queue[v]
        if v == end_vertex:
            break
        
        for w in edge_dictionary[v]:
            vwLength = dict_fd[v] + edge_dictionary[v][w]
            if w in dict_fd:
                if vwLength < dict_fd[w]:
                    raise(ValueError, "Dijkstra: found better path to already-final vertex")
            elif w not in queue or vwLength < queue[w]:
                queue[w] = vwLength
                dict_pred[w] = v
    
    return dict_fd, dict_pred

def floyd_warshall(adjacency_matrix):
    """ Floyd Warshall calculates the matrix of shortest paths between all pairs.

    Parameters
    ----------
    adjacency_matrix : np.array, square
        The adjacency matrix of the graph, on which the distances are being calculated.
    
    Returns
    -------
    dist : np.array
        The shortest path matrix as produced by floyd warshall
    """
    n = adjacency_matrix.shape[0]
    
    dist = np.empty((n, n))
    # Initialization
    for i in range(0,n):
        for j in range(0,n):
            if (i==j):
                dist[i,j] = 0.
            elif (adjacency_matrix[i,j] == 0):
                dist[i,j] = float("Inf")
            else:
                dist[i,j] = adjacency_matrix[i,j]

    # Calculation
    for k in range(0,n):
        for i in range(0,n):
            for j in range(0,n):
                if (dist[i,j] > dist[i,k] + dist[k,j]):
                    dist[i,j] = dist[i,k] + dist[k,j]

    return dist
