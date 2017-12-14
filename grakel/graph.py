""" A python file that implements classes, functions for graphs

"""
import warnings
import collections
import operator
import itertools

import numpy as np

import cvxopt.base
import cvxopt.solvers

from .tools import priority_dict, inv_dict, distribute_samples
from sklearn.svm import OneClassSVM

np.random.seed(238537)
cvxopt.solvers.options['show_progress'] = False

class graph(object):
    """ A general graph class that supports adjacency, dictionary formats while beeing memory/computationaly sustainable
    
    
    Parameters
    ----------
    initialization_object :  dict or array-like, square 
        The initialisation object for the graph. 
        If given a dictionary the input can be as follows:
            + 2-level nested dictionaries from edge symbols to weights
            + Dictionary of symbols to list of symbols (unweighted)
            + Dictionary of tuples to weights (weighted) #TODO
            + Set or List of tuples (unweighted) #TODO
            
        If given a array the input can be as follows:
            + array-like lists of lists #TODO
            + np.array
            + sparse matrix #TODO
        
    node_labels : dict
        A label dictionary corresponding to all vertices of the graph:
            + for adjacency matrix labels should be given to numbers starting 
              from 0 and ending in N-1, where the matrix has size N by N
            + for dictionary labels should correspond to all keys
    
    edge_labels : dict
        A labels dictionary corresponding to all edges of the graph keys: touples, value: label
    
    graph_format : str, valid_values={"dictionary", "adjacency", "all", "auto"}
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
        """ builds a graph structure given a supported graph representation
        
            g: - for an adjacency matrix input an ndarray
               - for edge_dictionary input a dictionary
        """

        self.shortest_path_mat = None
        self.label_group = None
        case = 0
        if g is not None:
            if type(g) is np.array or type(g) is np.ndarray:
                # Input is considered an adjacency matrix
                case = 1
                
                # Assign labels for nodes
                self.index_node_labels = node_labels

                # Assign labels for edges
                self.index_edge_labels = edge_labels
                
                if(self._format is "auto"):
                    self._format = "adjacency"
            elif type(g) is dict:
                # Input is considered as an edge dictionary
                case = 2
                
                # Assign labels for nodes
                self.node_labels = node_labels
                
                # Assign labels for edges
                self.edge_labels = edge_labels
                
                if(self._format is "auto"):
                    self._format = "dictionary"
            else:
                raise ValueError('Unsupported input type.\nValid input formats are numpy arrays and edge dictionaries.\nFor more information check the documentation')

        
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
        """ Changes the format of the graph from an existing to an other

            graph_format: Is the internal represantation of the graph object be a dictionary as a matrix, or both
                - for dictionary: "dictionary"
                - for adjacency_matrix: "adjacency"
                - for both: "all"
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
        """ Changes the format to include the desired

            graph_format: Is the internal represantation of the graph object be a dictionary as a matrix, or both
                - for dictionary: "dictionary"
                - for adjacency_matrix: "adjacency"
                - for both: "all"
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
        """ A method that constructs labels (if user does not provide)
            
            label_type: "vertex", "edge"
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
                                
    def convert_labels(self, target_format = "dictionary", purpose="all", init=False):
        """ A method that converts labels to a desired format 

            target_format: "dictionary" or "adjacency"
            purpose: "all", "edge" or "vertex"
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
        """ Returns the label of a vertex
            
            obj: a valid vertex
            label_type="vertex"
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
        """ A method to relabel the graph
            
            new_labels: the new_labels
            label_type: "vertex" or "edge"
            purpose: "dicitonary", "adjacency"
        """
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
        """ A method that builds and returns the shortest path matrix between all nodes

            algorithm_type: "auto" for "dictionary" or "all" format - Dijkstra
                                   for "adjacency" - Floyd Warshall
                            "dijkstra" - Dijkstra
                            "floyd_warshall" - Floyd Warshall
            labels: "vertex", "edge", "all"
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
            shortest_path_mat = floyd_warshall(self.adjacency_matrix,self.n)
        
        self.shortest_path_mat = shortest_path_mat
        if labels is "all":
            return shortest_path_mat, self.get_labels(), self.get_labels("edge")
        if labels is "edge":
            return shortest_path_mat, self.get_labels("edge")
        if labels is "vertex":
            return shortest_path_mat, self.get_labels()
            
    def get_labels(self, label_type="vertex", purpose="adjacency"):
        """ Return labels corresponding to the purpose
            
            purpose: if "adjacency" for indexes
                     if "dictionary" for nodes
            label_type: if "vertex" labels for vertices
                        if "edge" labels for edges
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
        """ A function that calculates the inverse dictionary
            for vertex labels (once)
            
            label_type: "vertex" or "edge"
            purpose: "dictionary" or "adjacency"
        """
        if not bool(self.label_group):
           self.label_group = dict()
        if not (label_type,purpose) in self.label_group:
           self.label_group[(label_type,purpose)] = dict()
        if not bool(self.label_group[(label_type,purpose)]):
            # calculate the label_group
            self.label_group[(label_type,purpose)] = inv_dict(self.get_labels(label_type,purpose))
        return self.label_group[(label_type,purpose)]
                
    def neighbors(self, vertex, purpose='any', with_weights = False):
        """ Find all neighbors of a vertex
            
            vertex: a valid vertex inside the matrix
                    that is not a sink

            with_weights: Determines the type of the output
                    if False: list of neighbor vertices
                    if True: dictionary between neighbor
                             vertices and edge labels 
            purpose: vertex is for 'adjacency' or 'dictionary'
                     'any': if format is 'all' default 'dictionary'
                            else the given format
        """
        if purpose in ['adjacency', 'dictionary', 'any']:
            if purpose == 'dictionary':
                self.desired_format('dictionary')
                case = True
            if purpose == 'adjacency':
                self.desired_format('adjacency')
                case = False
            if purpose == 'any':
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

    def _make_edsamic(self, vertices, all_out = True):
        """ A method to produce the edge symbols
            adjacency matrix correspondance dictionary
            
            all_out: output everything valuable
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
            if all_out:
                return edsamic, n, lov_sorted
            else:
                return edsamic, n
        else:
            warnings.warn('wrong format: edge dictionary must exist')
            return None


    def _import_adjacency(self, adjacency_matrix=None, init=True):
        """ A function that creates a graph object
            representation given its adjacency matrix

            adjacency_matrix: A square numpy ndarray
            save_matrix: an override variable that allows
                         the matrix to be stored internally
        """
        if adjacency_matrix is not None:
            # calculate graph size
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
        """ A function that creates a graph object
            representation given its edge dictionary

            edge_dictionary: input a dictionary as follows:
                                If (u,v) edge exists then edges["u"]["v"] has
                                weight of the edge pointing from u to v

            save_dictionary: an override variable that allows the matrix to be stored internally

        """
        if edge_dictionary is not None:
            # find vertices, refine dictionary
            vertices = set()
            edge_dictionary_refined = dict()
            for u in edge_dictionary.keys():
                vertices.add(u)
                # TODO support list of tuples
                if type(edge_dictionary[u]) is list:
                    vertices.union(set(edge_dictionary[u]))
                    edge_dictionary_refined[u] = dict(zip(edge_dictionary[u],len(edge_dictionary[u])*[1.0]))
                elif type(edge_dictionary[u]) is dict:
                    vertices.union(set(edge_dictionary[u].keys()))
                    edge_dictionary_refined[u] = edge_dictionary[u]
                else:
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
        """ Calculates the laplacian of the given graph
            save: optional parameter to store the matrix
        """ 
        if self.laplacian_graph is not None:
            laplacian_graph = self.laplacian_graph
        else:
            self.desired_format("adjacency", warn=True)
            laplacian_graph = laplacian(self.adjacency_matrix, self.n)
            
            if save:
                self.laplacian_graph = laplacian_graph
        return laplacian_graph
    
    def calculate_lovasz_theta(self, from_scratch=False):
        """ Calculates the lovasz theta for the given graph.
            
            from_scratch: boolean flag that chooses with False to
                          recalculate the lovasz number if it has
                          already been calculated.
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
            
            from_scratch: boolean flag that chooses with False to
                          recalculate the svm theta number if it has
                          already been calculated.
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
        """ A function useful for calculating a graph metric
            as the lovasz theta kernel or the svm-theta kernel,
            on a set of randomly sampled subgraphs, producing
            a dictionary of subgraph levels and sets
                        
            n_samples: the number of samples that will be sampled
            subsets_size_range: A touple having the min and the max subset size
            from_scratch: boolean flag that chooses with False to
                          recalculate the lovasz number if it has
                          already been calculated.
            save: boolean flag that determines if the output will be saved
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
            purpose: "adjacency", "dictionary", "any"
        """
        if purpose not in ["adjacency", "dictionary", "any"]:
            raise ValueError('purpose is either "adjacency" of "dictionary"')

        if purpose == 'any':
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
        """ A method that returns an iterable of edges as touples.
            purpose: "adjacency", "dictionary"
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
        """ A format agnostic method that returns the adjacency matrix
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
        """ A method that returns the number of vertices
            for any existing format
        """
        if self._format in ['all','adjacency']:
            return self.n
        else:
            return len(self.vertices)
            
    def produce_neighborhoods(self, r=3, purpose="adjacency", with_distances=False, d=-1):
        """ Calculates neighborhoods for each node
            of a Graph up to a depth c.
            
            G: a graph type object
            r: neighborhood depth
            purpose: node symbols for "adjacency", "dictionary"
            with_distances: a flag that defines if we need to calculate
                            BFS distances for each pair.
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
        """A method that creates a graph object subgraph
           in the same format as the original graph
    
           vertices: an iterable of vertices convertable to a set
                     vertices a subset of verices of the original graph
                     
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
        
def laplacian(A, n=-1):
    """ Calculates the laplacian given the adjacency matrix
    
        A: a numpy array of a square matrix corresponding
           to the adjacency matrix of a graph
        n: x dimension of the matrix [int] - optional
    """        

    if (n<=0):
        n = A.shape[0]
        
    laplacian_graph = np.zeros(shape=(n, n))
    for i in range(0, n):
        for j in range(0, n):
            if i!=j and A[i,j]!=.0:
                laplacian_graph[i,j] = -A[i,j]
                laplacian_graph[i,i]+= A[i,j]
    return laplacian_graph

def dijkstra(edge_dictionary, start_vertex, end_vertex=None):
    """ Implementation of the dijkstra algorithm
    
        Find shortest paths from the start vertex to all
        vertices nearer than or equal to the end_vertex.

        The majority of this function code came from:
        http://code.activestate.com/recipes/119466-dijkstras-algorithm-for-shortest-paths/ 

        edge_dictionary: input a dictionary as follows:
                             If (u,v) edge exists then edges["u"]["v"] has
                             weight of the edge pointing from u to v

        start_vertex: A valid vertex symbol that exists inside edge_dictionary as a key
        
        end_vertex: A valid vertex symbol that exists inside edge_dictionary as a key
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

def floyd_warshall(adjacency_matrix, n=-1):
    """ Floyd Warshall calculates the matrix of shortest paths between all pairs

        adjacency_matrix : a square nd array
        n : length of the one dimension of the ndarray [optional]
    """
    
    if n < 0:
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
