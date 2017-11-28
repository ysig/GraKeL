""" A python file that implements classes, functions for graphs

"""
import warnings
import collections
import operator

import numpy as np

import cvxopt.base
import cvxopt.solvers

from .tools import priority_dict, inv_dict, distribute_samples
from sklearn.svm import OneClassSVM

np.random.seed(238537)
cvxopt.solvers.options['show_progress'] = False

class graph(object):
    """ Definition of a graph

    G = (V,E,W)
    
    V: a set of vertices
    E: a set of edges between vertices of V
    W: a set of weights between each edge in E

    """
   
    # n is the number of vertices
    n = -1

    # A represantation of the graph as a dictionary
    # between vertices.
    # If (u,v) edge exists then edges["u"]["v"] has
    # weight of the edge pointing from u to v
    
    edge_dictionary = dict()
    # A set of vertices corresponding to the edge dictionary
    # for fast indexing
    vertices = set()
    
    # adjacency matrix corresponding to current graph
    adjacency_matrix = np.empty(0)

    # a labels dictionary
    node_labels = dict()
    edge_labels = dict()

    # edge_dictionary_symbols_adjacency_matrix_index_correspondance
    # A dictionary which links adjacency matrix numbering 
    # with edge labels of input or existing dictionary
    edsamic = dict()

    # A direct label holder for adjoint idxs when having dictionary with symbols
    index_node_labels = dict()    
    index_edge_labels = dict()
    
    # a variable holding the shortest path matrix
    shortest_path_mat = None

    # label_group: an inverse map between the labels and the nodes
    label_group = dict()
    
    # keeps the laplacian graph
    laplacian_graph = None
    
    # an edge labels dictionary
    edge_labels = dict()

    # lovasz theta
    lovasz_theta =  None
    
    # svm_theta
    svm_theta = None
    
    # metrics subgraphs dictionary
    metric_subgraphs_dict = dict()
    
    def __init__(self, initialization_object=None, node_labels=None, edge_labels=None, graph_format="auto"):
        
        """ Creates a new graph object

            initialization_object: An input object given to initialise the given graph
                - for an adjacency matrix input a square numpy array
                - for edge_dictionary input a dictionary as follows:
                    If (u,v) edge exists then edges["u"]["v"] has
                    weight of the edge pointing from u to v
                    If no weights for an edge (u,v) then edges[u]
                    must be a list and v must be inside                   

            node_labels: A label dictionary corresponding to all vertices of the graph
                - for adjacency matrix labels should be given to numbers starting
                  from 0 and ending in N-1, where the matrix has size N by N
                - for dictionary labels should correspond to all keys
            edge_labels: A labels dictionary corresponding to all edges of the graph
                         keys: touples, value: label
            graph_format: Is the internal represantation of the graph object be a dictionary as a matrix, or both
                - for dictionary: "dictionary"
                - for adjacency_matrix: "adjacency"
                - for both: "all"
                - for the current_format (if existent): "auto"
        """


        if graph_format in ["adjacency","dictionary","auto","all"]:
            self._format = graph_format
            if (initialization_object is not None):
                self.build_graph(initialization_object, node_labels, edge_labels)
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
        """ Changes the format of the graph
            from an existing to an other

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

    def desired_format(self,graph_format):
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
                self.change_format("all")
        elif graph_format is "adjacency":
            if self._format not in ["all","adjacency"]:
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
                    warnings.warn('no labels to convert from, for the given purpose')
            else:
                warnings.warn('labels already defined for that format - nothing to convert')
        elif (target_format == "dictionary"):
            if (self._format != "dictionary" or init):
                if bool(self.node_labels) and purpose in ['all', 'vertex']:
                    warnings.warn('overriding existing node labels, for dictionary symbols')
                    
                if bool(self.edge_labels) and purpose in ['all', 'edges']:
                    warnings.warn('overriding existing edge labels, for dictionary symbols')
                if purpose in ["all", "vertex"] and bool(self.index_node_labels):
                    self.node_labels = self.index_node_labels
                if purpose in ["all", "edges"] and bool(self.index_edge_labels):
                    self.edge_labels = self.index_edge_labels
                else:
                    warnings.warn('no labels to convert from, for the given purpose')
            else:
                warnings.warn('labels already defined for the given format')
    
    def label(self, obj, label_type="vertex", purpose="dictionary"):
        """ Returns the label of a vertex
            
            obj: a valid vertex
            label_type="vertex"
        """
        if (label_type == "vertex"):
            vertex = obj
            if purpose == "dictionary":
                labels = self.node_labels
            elif purpose == "adjacency":
                labels = self.index_node_labels
            else:
                raise ValueError('unrecognized purpose')
                
            if not bool(labels):
                if (vertex in labels):
                    return labels[vertex]
                else:
                    raise ValueError('no label assigned to this vertex symbol')
            else:
                warnings.warn('labels are not set - default value returned')
                return vertex
        elif (label_type == "edge"):
            edge = obj
            if purpose == "dictionary":
                labels = self.edge_labels
            elif purpose == "adjacency":
                labels = self.index_edge_labels
            else:
                raise ValueError('unrecognized purpose')
                
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

        
    def build_shortest_path_matrix(self, algorithm_type="auto", clean=False):
        """ A method that builds and returns the shortest path matrix between all nodes

            algorithm_type: "auto" for "dictionary" or "all" format - Dijkstra
                                   for "adjacency" - Floyd Warshall
                            "dijkstra" - Dijkstra
                            "floyd_warshall" - Floyd Warshall
        """
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
            self.desired_format("all")
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

            shortest_path_labels = self.get_labels()
        elif algorithm_type is "floyd_warshall":
            self.desired_format("adjacency")
            shortest_path_mat = floyd_warshall(self.adjacency_matrix,self.n)
            shortest_path_labels = self.get_labels()

        self.shortest_path_mat = shortest_path_mat
        return shortest_path_mat, shortest_path_labels
    
    def get_labels(self, label_type="vertex", purpose="adjacency"):
        """ Return labels corresponding to the purpose
            
            purpose: if "adjacency" for indexes
                     if "dictionary" for nodes
            label_type: if "vertex" labels for vertices
                        if "edge" labels for edges
        """
        if (purpose == "adjacency"):
            self.desired_format("adjacency")
            if label_type == "vertex":
                if not bool(self.index_node_labels):
                    self.construct_labels(label_type, purpose)
                return self.index_node_labels
            elif label_type == "edge":
                if not bool(self.index_edge_labels):
                    self.construct_labels(label_type, purpose)
                return self.index_edge_labels
        elif (purpose == "dictionary"):
            self.desired_format("dictionary")
            if label_type == "vertex":
                if not bool(self.node_labels):
                    self.construct_labels(label_type, purpose)
                return self.node_labels
            elif label_type == "edge":
                if not bool(self.edge_labels):
                    self.construct_labels(label_type, purpose)
                return self.edge_labels
        else:
            raise ValueError('unsupported label purpose')
    
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
                
    def neighbours(self, vertex, with_weights = False):
        """ Find all neighbours of a vertex
            
            vertex: a valid vertex inside the matrix
                    that is not a sink

            with_weights: Determines the type of the output
                    if False: list of neighbor vertices
                    if True: dictionary between neighbour
                             vertices and edge labels 
            
        """
        if self._format in ['dictionary','all']:
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
                if not with_weights:
                    out = list()
                    for i in range(0,self.n):
                        if i == idx:
                            continue
                        if(self.adjacency_matrix[idx,i] != .0):
                            out.append(i)
                    return out
                else:
                    out = dict()
                    for i in range(idx,self.n):
                        if i == idx:
                            continue
                        element = self.adjacency_matrix[idx,i]
                        if(element != .0):
                            out[i] = element
                    return edge_dictionary[vertex]
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


    def _import_adjacency(self, adjacency_matrix=None):
        """ A function that creates a graph object
            representation given its adjacency matrix

            adjacency_matrix: A square numpy ndarray
            save_matrix: an override variable that allows the matrix to be stored internally
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
            edge_dictionary = dict()
            for i in range(0,n):
                for j in range(0,n):
                    if (adjacency_matrix[i,j] != 0):
                        if i not in edge_dictionary:
                            edge_dictionary[i] = dict()
                        edge_dictionary[i][j] = adjacency_matrix[i,j]
            self.vertices = vertices
            self.edge_dictionary = edge_dictionary

            # Add labels
            self.convert_labels(target_format="dictionary", purpose="all", init=True)
            
            # Prune not interesting format
            if self._format is "dictionary":
                self.n = -1
                self.adjacency_matrix = None
                self.edsamic = None
                self.index_node_labels = None
                self.index_edge_labels = None

            
    def _import_dictionary(self, edge_dictionary=None):
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
                if type(edge_dictionary[u]) is list:
                    vertices.union(set(edge_dictionary[u]))
                    edge_dictionary_refined[u] = dict(zip(edge_dictionary[u],len(edge_dictionary[u])*[1.0]))
                elif type(edge_dictionary[u]) is dict:
                    vertices.union(set(edge_dictionary[u].keys()))
                    edge_dictionary_refined[u] = edge_dictionary[u]
                else:
                    raise ValueError('unsupported edge_dictionary format')
            
            # Save dictionary, vertices @self if needed        
            if (self._format is "all" or self._format is "dictionary"):
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
            
            if self._format is "all":
                self.edsamic = edsamic
            
            # Add labels
            self.convert_labels(target_format="adjacency", purpose="all", init=True)
            
            # Prune for a certain format
            if self._format is "adjacency":
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
            self.desired_format("adjacency")
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
        if self._format is "dictionary":
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
        if self._format is "dictionary":
            nv = len(vertices)
            lov_sorted = sorted(list(vertices))
            vertices = {v: i for (i,v) in enumerate(lov_sorted)}
            
            K = np.zeros(shape=(nv,nv))
            idx = zip(*[(vertices[v],vertices[n]) for v in lov_sorted for n in self.neighbours(v)])
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
        if metric_type is "svm":
            metric = lambda G: G.calculate_svm_theta()
        elif metric_type is "lovasz":
            metric = lambda G: G.calculate_lovasz_theta()            
        else:
            raise ValueError('unsupported metric type')
            
        min_ss, max_ss = subsets_size_range[0], subsets_size_range[1]
        # If precalculated and the user wants it, output
        if not from_scratch:
            if bool(self.metric_subgraphs_dict):
                if (metric_type, n_samples, min_ss, max_ss) in self.metric_subgraphs_dict:
                    return self.metric_subgraphs_dict[(metric_type, n_samples, min_ss, max_ss)]
        
        
        self.desired_format("adjacency")
        
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
