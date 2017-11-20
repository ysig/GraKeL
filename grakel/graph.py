""" A python file that implements classes, functions for graphs

"""
import collections
import operator
import numpy as np

from .tools import priority_dict, inv_dict

class graph(object):
    """ Definition of a graph

    G = (V,E,W)
    
    V: a set of vertices
    E: a set of edges between vertices of V
    W: a set of weights between each edge in E

    """
   
    # n is the number of vertices
    n = 0

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
    labels = dict()

    # edge_labels_adjacency_matrix_correspondance_dictionary
    # A dictionary which links adjacency matrix numbering 
    # with edge labels of input or existing dictionary
    elamcd = dict()

    # A direct label holder for adjoint idxs when having dictionary with symbols
    index_labels = dict()    

    # a variable holding the shortest path matrix
    shortest_path_mat = None

    # label_group: an inverse map between the labels and the nodes
    label_group = dict()

    def __init__(self, initialization_object=None, labels=None, graph_format="auto"):
        
        """ Creates a new graph object

            initialization_object: An input object given to initialise the given graph
                - for an adjacency matrix input a square numpy array
                - for edge_dictionary input a dictionary as follows:
                    If (u,v) edge exists then edges["u"]["v"] has
                    weight of the edge pointing from u to v
                    If no weights for an edge (u,v) then edges[u]
                    must be a list and v must be inside                   

            labels: A label dictionary corresponding to all vertices of the graph
                - for adjacency matrix labels should be given to numbers starting
                  from 0 and ending in N-1, where the matrix has size N by N
                - for dictionary labels should correspond to all keys

            graph_format: Is the internal represantation of the graph object be a dictionary as a matrix, or both
                - for dictionary: "dictionary"
                - for adjacency_matrix: "adjacency"
                - for both: "all"
                - for the current_format (if existent): "auto"
        """


        if graph_format in ["adjacency","dictionary","auto","all"]:
            self._format = graph_format
            self.build_graph(initialization_object,labels)
        else:
            pass
            # Raise an exception if graph format is not valid?

    def build_graph(self, g, labels=None):
        """ builds a graph structure given a supported graph representation
        
            g: - for an adjacency matrix input an ndarray
               - for edge_dictionary input a dictionary
        """

        self.labels = labels
        self.shortest_path_mat = None
        self.label_group = None

        case = 0
        if g is not None:
            if type(g) is np.array or type(g) is np.ndarray:
                # Input is considered an adjacency matrix
                case = 1
                if(self._format is "auto"):
                    self._format = "adjacency"
            elif type(g) is dict:
                # Input is considered as a edge dictionary
                case = 2
                if(self._format is "auto"):
                    self._format = "dictionary"
            else:
                pass
                # Raise exception: "Unsupported input type"

        # If graph is of one type prune the other        
        if self._format is "adjacency":
            edge_dictionary = None

        elif self._format is "dictionary":
            adjacency_matrix = None

        if (case==1):
            self._import_adjacency(g)
        elif (case==2):
            self._import_dictionary(g)

        if self.labels is None:
            if self._format in ["dictionary","all"]:
                nodes = sorted(list(self.vertices))
                self.labels = dict(zip(nodes,nodes))
            else:
                nodes = list(range(0,self.n))
            self.labels = dict(zip(nodes,nodes)) 

    def change_format(self, graph_format):
        """ Changes the format of the graph
            from an existing to an other

            graph_format: Is the internal represantation of the graph object be a dictionary as a matrix, or both
                - for dictionary: "dictionary"
                - for adjacency_matrix: "adjacency"
                - for both: "all"
        """
        if graph_format not in ["all", "dictionary", "adjacency"]:
            pass
            # Raise exception ?
        else:
            if (graph_format is not self._format):
                past_format = self._format
                self._format = graph_format
                if past_format is "adjacency":
                    self._import_adjacency()
                    if self._format is "dictionary":
                        self.n = 0
                        self.adjacency_matrix = None
                        self.elamcd = None
                elif past_format is "dictionary":
                    self._import_dictionary()
                    if self._format is "adjacency":
                        self.edge_dictionary = None
                        self.vertices = None
                else:
                    if self._format is "dictionary":
                        self.n = 0
                        self.adjacency_matrix = None
                        self.elamcd = None
                    else:
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
            if self._format not in ["all","edge_dictionary"]:
                self.change_format("all")
        elif graph_format is "adjacency":
            if self._format not in ["all","adjacency"]:
                self.change_format("all")

    def get_label_group(self):
        """ A function that calculates the inverse dictionary
            for labels (once)
        """
        if not bool(self.label_group):
            label_group = dict()
            # calculate the label_group
            if not bool(self.labels):
                if (self._format == "adjacency"):
                    for i in range(0,self.n):
                        label_group[i] = [i]
                else:
                    for k in vertices.set():
                        label_group[k] = [k]
                self.label_group = label_group
            else:
                self.label_group = inv_dict(self.labels)
        return self.label_group

    def label(self,vertex):
        """ Returns the label of a vertex
            vertex: a valid vertex
        """
        if not bool(self.labels):
            if (vertex in self.labels):
                return self.labels[vertex]
            else:
                pass
                # Raise exception: No label assigned to thise vertex symbol?
        else:
            return vertex


    def neighbours(self, vertex, with_weights = False):
        """ Find all neighbors of a vertex
            
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
                # Raise warning: vertex not inside matrix ?
                return None
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
                pass
                # Raise exception: "item with index ",idx," does not exist"?

    def relabel(self,new_labels,labels_for="dictionary"):
        """ Adds new labels dictionary """
        if new_labels is None or not bool(new_labels):
            # Raise warning: User must provide new labels?
            pass
        elif labels_for is "dictionary":
            if self._format in "dictionary":
                self.labels = new_labels
                if self.shortest_path_mat is not None:
                    lov = list(self.vertices)
                    shortest_path_mat_dim = len(lov)
                    lov_sorted = sorted(lov)

                    # shortest path labels
                    sp_labels =  dict()
                    if bool(self.splabels):
                        tmp = [self.labels[k] for k in lov_sorted]
                        sp_labels = dict(zip(list(range(0,shortest_path_mat_dim)),tmp))
            elif self._format in "all":
                self.labels = new_labels
                if bool(self.index_labels):
                    for k in elamcd.keys():
                        self.index_labels[elamcd[k]] = new_labels[k]
                    if bool(self.sp_labels):
                        self.sp_labels = self.index_labels
            else:
                pass
                # Raise exception?
        elif labels_for is "adjacency":
            if self._format is "all":
                if bool(self.index_labels):
                    self.index_labels = new_labels
                    # move from index labels to dictionary key labels
                    ielamcd = inv_dict(self.elamcd)
                    for k in ielamcd.keys():
                        self.labels[ielamcd[k]] = new_labels[k]
                else:
                    self.labels = new_labels
                if bool(self.sp_labels):
                    self.sp_labels = self.labels
            elif self._format is "adjacency":
                self.labels = new_labels
                if bool(self.sp_labels):
                    self.sp_labels = self.labels
            else:
                pass
                # Raise exception?
        else:
            pass
            # Raise exception?
                
                
    def build_shortest_path_matrix(self, algorithm_type="auto", clean=False):
        """ A method that builds and returns the shortest path matrix between all nodes

            algorithm_type: "auto" for "dictionary" or "all" format - Dijkstra
                                   for "adjacency" - Floyd Warshall
                            "dijkstra" - Dijkstra
                            "floyd_warshall" - Floyd Warshall
        """
        if clean:
            self.shortest_path_mat = None
            self.shortest_path_labels = None
        
        if self.shortest_path_mat is not None:
            return self.shortest_path_mat, self.shortest_path_labels

        # Assign the desired algorithm
        if algorithm_type is "auto":
            if self._format in ["all","dictionary"]:
                algorithm_type = "dijkstra"
            elif self._format is "adjacency":
                algorithm_type = "floyd_warshall"
        
        if algorithm_type is "dijkstra":
            # Prepare for algorithm implementation
            self.desired_format("dictionary")
            if self._format is "dictionary":
                
                # vertices to list
                lov = list(self.vertices)
                shortest_path_mat_dim = len(lov)
                lov_sorted = sorted(lov)

                # make indexes - same with elamcd           
                indexes = dict(zip(lov_sorted,list(range(0,shortest_path_mat_dim))))

                # shortest path labels
                sp_labels =  dict()
                if bool(self.labels):
                    tmp = [self.labels[k] for k in lov_sorted]
                    sp_labels = dict(zip(list(range(0,shortest_path_mat_dim)),tmp))
                else:
                    sp_labels = dict(zip(list(range(self.n)),list(range(self.n))))

            elif self._format is "all":
                if (bool(self.elamcd)):
                    indexes = self.elamcd
                else:
                    indexes = dict(zip(list(range(self.n)),list(range(self.n))))
                sp_labels = self.get_labels()
                shortest_path_mat_dim = len(indexes)
            
            # calculate shortest path matrix
            shortest_path_mat = np.full([shortest_path_mat_dim, shortest_path_mat_dim],float("Inf"))
            for k in indexes.keys():
                dict_fd, _ = dijkstra(self.edge_dictionary,k)
                for s in dict_fd.keys():
                    shortest_path_mat[indexes[k],indexes[s]] = dict_fd[s]

            shortest_path_labels = sp_labels
        elif algorithm_type is "floyd_warshall":
            self.desired_format("adjacency")
            shortest_path_mat = floyd_warshall(self.adjacency_matrix,self.n)
            shortest_path_labels = self.get_labels()

        self.shortest_path_mat = shortest_path_mat
        self.shortest_path_labels = shortest_path_labels
        return shortest_path_mat, shortest_path_labels
    
    def get_labels(self, purpose="adjacency"):
        """ Return labels corresponding to the purpose
            
            purpose: if "adjacency" for indexes
                     if "dictionary" for nodes
        """
        if (purpose == "adjacency"):
            if bool(self.index_labels):
                return self.index_labels
            else:
                return self.labels
        else:
            return self.labels

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
                pass
                # Raise an exception if matrix is not squared?
            # import_adjacency
            if (self._format is "all" or self._format is "adjacency"):
                self.adjacency_matrix = adjacency_matrix
            self.n = n
        else:
            n = self.n
            adjacency_matrix = self.adjacency_matrix

        # construct a dictionary out of the adjacency
        if self._format is "all" or self._format is "dictionary":
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
                    pass
                    # Raise exception unsupported edge_dictionary format ?
            
            # Save dictionary, vertices @self if needed        
            if (self._format is "all" or self._format is "dictionary"):
                self.vertices = vertices
                self.edge_dictionary = edge_dictionary_refined
        else:
            vertices = self.vertices
            edge_dictionary_refined = self.edge_dictionary

        # Create and store the adjacency matrix
        if self._format in ["adjacency","all"]:
            # vertices to list
            lov = list(vertices)

            # length is adjacency matrix size
            self.n = len(lov)

            # sort lov indexes for labeling
            lov_sorted = sorted(lov)

            # edge_labels_adjacency_matrix_correspondance_dictionary            
            elamcd = dict(zip(lov_sorted,list(range(0,self.n))))

            # index labels
            index_labels = dict()
            if bool(self.labels):
                tmp = [self.labels[k] for k in lov_sorted]
                index_labels = dict(zip(list(range(0,self.n)),tmp))

            # Initialize adjacency_matrix
            adjacency_matrix = np.zeros(shape = (self.n,self.n))
            
            # Produce and save adjacency matrix
            for k in edge_dictionary_refined.keys():
                for l in edge_dictionary_refined[k].keys():
                    adjacency_matrix[elamcd[k],elamcd[l]] = edge_dictionary_refined[k][l]
                
            self.adjacency_matrix = adjacency_matrix
        
            if self._format is "adjacency":
                self.labels = index_labels
            elif self._format is "all":
                self.index_labels = index_labels
                self.elamcd = elamcd


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
 
# Exceptions
# ----------
#
# Always use class exceptions instead of the old-style string exceptions.
class MessageError(Exception):
    """ Base class for errors in the email package.

    Exceptions are classes too! Hence exception names are CamelCase.

    """
    
    pass
