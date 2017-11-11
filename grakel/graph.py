""" A python file that implements classes, functions for graphs

"""

import numpy as np

from tools import priority_dict, inv_dict

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
    edges_dictionary = dict()
    # A set of vertices corresponding to the edge dictionary
    # for fast indexing
    vertices = set()
    
    # adjacency matrix corresponding to current graph
    adjacency_matrix = np.ndarray.empty()

    # a labels dictionary
    labels = dict()

    # edge_labels_adjacency_matrix_correspondance_dictionary
    # A dictionary which links adjacency matrix numbering 
    # with edge labels of input or existing dictionary
    elamcd = dict()

    # a variable holding the shortest path matrix
    shortest_path_mat = None

    # label_group: an inverse map between the labels and the nodes
    label_group = dict()

    def __init__(self, initialization_object=None, labels=None, graph_format="dictionary"):
        
        """ Creates a new graph object

            initialization_object: An input object given to initialise the given graph
                - for an adjacency matrix input a square numpy ndarray
                - for edge_dictionary input a dictionary as follows:
                    If (u,v) edge exists then edges["u"]["v"] has
                    weight of the edge pointing from u to v

            labels: A label dictionary corresponding to all vertices of the graph
                - for adjacency matrix labels should be given to numbers starting
                  from 1 and ending in N, where the matrix has size N by N
                - for dictionary labels should correspond to all keys

            graph_format: Is the internal represantation of the graph object be a dictionary as a matrix, or both
                - for dictionary: "dictionary"
                - for adjacency_matrix: "adjacency"
                - for both: "all"
                
        """


        if self.format in ["adjacency","dictionary","all"]:
            self.format = graph_format
            self.build_graph(initialization_object,labels)
        else:
            # Raise an exception if graph format is not valid?
        
    def build_graph(self, g, labels=None):
        """ builds a graph structure given a supported graph representation
        
            g: - for an adjacency matrix input an ndarray
               - for edge_dictionary input a dictionary
        """

        self.labels = labels
        self.shortest_path_mat = None
        self.label_group = None

        # If graph is of one type prune thee other
        if self.format is "adjacency":
            edge_dictionary = None

        elif self.format is "dictionary":
            adjacency_matrix = None

        if g is not None:
            if type(g) is np.ndarray:
                # Input is considered an adjacency matrix
                                
                self._import_adjacency(g)

            elif type(g) is dictionary:
                # Input is considered as a edge dictionary

                self._import_dictionary(g)

            else:
                # Raise exception: "Unsupported input type"

    def change_format(self, graph_format):
        """ Changes the format of the graph
            from an existing to an other

            graph_format: Is the internal represantation of the graph object be a dictionary as a matrix, or both
                - for dictionary: "dictionary"
                - for adjacency_matrix: "adjacency"
                - for both: "all"
        """
        if graph_format not in ["all", "dictionary", "adjacency"]
            # Raise exception ?
        else:
            if (graph_format is not self.format):
                if self.format is "adjacency":
                    self._import_adjacency(self.adjacency,False)
                    if graph_format is "dictionary"
                        self.n = 0
                        self.adjacency_matrix = None
                        self.elamcd = None
                elif self.format is "dictionary":
                    self._import_dictionary(self.edge_dictionary,False)
                    if graph_format is "adjacency":
                        self.edge_dictionary = None
                        self.vertices = None
                else:
                    if graph_format is "dictionary":
                        self.n = 0
                        self.adjacency_matrix = None
                        self.elamcd = None
                    else:
                        self.edge_dictionary = None

    def get_label_group(self):
        """ A function that calculates the inverse dictionary
            for labels (once)
        """
        if self.label_group is None:
            label_group = dict()
            # calculate the label_group
            if (self.labels == {} or self.labels is None):
                if (self.format == "adjacency"):
                    for i in range(0,n):
                        label_group[i] = [i]
                else
                    for k in vertices.set():
                        label_group[k] = [k]
                self.label_group = label_group
            else:
                self.label_group = inv_dict(self.labels)
        else:
            return self.label_group

    def label(self,vertex):
        """ Returns the label of a vertex
            vertex: a valid vertex
        """
        if(self.labels == {} or self.labels is None):
            return vertex
        else:
            if (vertex in self.labels):
                return self.labels[vertex]
            else:
                # Raise exception


    def neighbours(self, vertex, with_weights = False):
        """ Find all neighbors of a vertex
            
            vertex: a valid vertex inside the matrix
                    that is not a sink

            with_weights: Determines the type of the output
                    if False: list of neighbor vertices
                    if True: dictionary between neighbour
                             vertices and edge labels 
            
        """
        if self.format in ['dictionary','all']:
            if vertex in edges_dictionary:
                if not with_weights:
                    return edge_dictionary[vertex].keys()
                else:
                    return edge_dictionary[vertex]
            else:
                # Raise exception: vertex not inside matrix
                return None
        else:
            idx = int(vertex)
            if 1 < idx < self.n:
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
                    for i in range(0,self.n):
                        if i == idx:
                            continue
                        element = self.adjacency_matrix[idx,i]
                        if(element != .0):
                            out[i] = element
                    return edge_dictionary[vertex]
            else:
                # Raise exception:
                # "item with index ",idx," does not exist"

    def build_shortest_path_matrix(self):
        """ A method that builds and returns the shortest path matrix between all nodes """

        # If dijkstra algorithm is going to be applied
        # labels must be numbered and added to a dictionary
        if self.format is "dictionary":
            temp_labels = dict()
            edge_labels_count = 0
            for k in edge_dictionary.keys():
                if k not in temp_labels:
                    temp_labels[k] = edge_labels_count
                    edge_labels_count += 1
                for l in edge_dictionary[k].keys():
                    if l not in temp_labels: 
                        temp_labels[k] = edge_labels_count
                        edge_labels_count += 1
             shortert_path_mat_dim = edge_labels_count - 1

        elif self.format is "all":
            temp_labels = self.elamcd
            shortert_path_mat_dim = len(self.elamcd.keys())

        # Calculate the shortest path matrix
        if self.format is in ["dictionary","all"]:
           shortest_path_mat = np.empty([short_path_mat_dim, short_path_mat_dim])
           for k in temp_labels.keys():
               (dict_fd,_) = dijkstra(self.edge_dictionary,k):
               for s in dict_fd.keys():
                   shortest_path_mat[temp_labels[k],temp_labels[s]] = dict_fd[s]

            # Create assign labels to shortest path idxs
            shortest_path_labels = dict()
            if (self.labels != []):
                for k in temp_labels.keys():
                    shortest_path_labels[temp_labels[k]] = labels[temp_labels[k]]
            else:
                shortest_path_labels = temp_labels
        else:
            shortest_path_mat = floyd_warshall(self.adjacency_matrix,self.n)
            shortest_path_labels = labels

        self.shortest_path_mat = shortest_path_mat
        return shortest_path_mat, shortest_path_labels

    def _import_adjacency(self,adjacency_matrix,save_matrix=True):
        """ A function that creates a graph object
            representation given its adjacency matrix

            adjacency_matrix: A square numpy ndarray
            save_matrix: an override variable that allows the matrix to be stored internally
        """

        # calculate graph size
        self.n = g.shape[0]

        # Raise an exception if matrix is not squared?

        # import_adjacency
        if (save_matrix is True) and (self.format is "all" or self.format is "adjacency"):
            self.adjacency_matrix = adjacency_matrix

        # construct a dictionary out of the adjacency
        if self.format is "all" or self.format is "dictionary":
            vertices = set(list(range(0, self.n)))
            edge_dictionary = dict()
            for i in range(0,self.n):
                for j in range(0,self.n):
                    if (adjacency_matrix[i,j] != 0):
                        if i not in edge_dictionary:
                            edge_dictionary[i] = OrderedDict()
                        edge_dictionary[i][j] = adjacency_matrix[i,j]
            self.vertices = vertices
            self.edge_dictionary = edge_dictionary
                    
    def _import_dictionary(self,edge_dictionary,save_dictionary=True):
        """ A function that creates a graph object
            representation given its edge dictionary

            edge_dictionary: input a dictionary as follows:
                                If (u,v) edge exists then edges["u"]["v"] has
                                weight of the edge pointing from u to v

            save_dictionary: an override variable that allows the matrix to be stored internally

        """

        # Save the dictionary internally if that is the case
        if (save_dictionary is True) and (self.format is "all" or self.format is "dictionary"):
            vertices = dict()
            edge_dictionary = dict()
            for u in edge_dictionary.keys():
                edge_dictionary[u] = OrderedDict(edge_dictionary[u])
                vertices.add(u)
                for v in edge_dictionary[u].keys():
                    vertices.add(v)
            self.vertices = vertices
            self.edge_dictionary = edge_dictionary

        # Create and store the adjacency matrix
        if self.format is "all" or self.format is "adjacency":  
            # First count edge labels
            edge_labels_count = 0
            # edge_labels_adjacency_matrix_correspondance_dictionary
            elamcd = dict()
            for k in edge_dictionary.keys():
                if k not in elamcd:
                    elamcd[k] = edge_labels_count
                    edge_labels_count += 1
                for l in edge_dictionary[k].keys():
                    if l not in elamcd: 
                        elamcd[k] = edge_labels_count
                        edge_labels_count += 1

            # Initialize adjacency_matrix
            self.n = edge_labels_count-1
            adjacency_matrix = np.zeros(shape = (self.n,self.n))
            
            # Produce and save adjacency matrix
            for k in edge_dictionary.keys():
                for l in edge_dictionary[k].keys():
                    adjacency_matrix[elamcd[k],elamcd[l]] = edge_dictionary[k][l]
            self.adjacency_matrix = adjacency_matrix


def dijkstra(edge_dictionary,start_vertex,end_vertex=None):
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
                    raise ValueError, "Dijkstra: found better path to already-final vertex"
            elif w not in queue or vwLength < queue[w]:
                queue[w] = vwLength
                dict_pred[w] = v
    
    return (dict_fd,dict_pred)

def floyd_warshall(adjacency_matrix,n = -1):
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


# Whitespace in Expressions and Statements
def example():
    """ How to use whitespace in expressions and statements. """

    # Around brackets, parentheses and braces:-
    #
    # BAD
    spam( ham[ 1 ], { eggs: 2 } )

    # GOOD
    spam(ham[1], {eggs : 2})

    # Around a comma, semicolon, or colon, as in (although multi-statement
    # lines are STRONGLY discouraged!):-
    #
    # BAD
    if x == 4 : print x , y ; x , y = y , x

    # GOOD
    if x == 4: print x, y; x, y = y, x

    # Function/method calls.
    #
    # BAD
    spam (1)

    # GOOD
    spam(1)

    # Indexing/slicing.
    #
    # BAD
    dict ['key'] = list [index]

    # GOOD
    dict['key'] = list[index]

    # Groups of assignment statements.
    #
    # BAD
    x             = 1
    y             = 2
    long_variable = 3

    # GOOD
    x = 1
    y = 2
    long_variable = 3

    # These binary operators should aways surrounded with a single space.
    #
    # assignment (=), comparisons (==, <, >, !=, <>, <=, >=, in, not in, is,
    # is not), Booleans (and, or, not).

    # Use your better judgment for the insertion of spaces around arithmetic
    # operators.  Always be consistent about whitespace on either side of a
    # binary operator.
    #
    # Some examples:
    i = i+1
    submitted = submitted + 1
    x = x*2 - 1
    hypot2 = x*x + y*y
    c = (a+b) * (a-b)
    c = (a + b) * (a - b)

    return


# Don't use spaces around the '=' sign when used to indicate a keyword
# argument or a default parameter value.  For instance:
def complex(real, imag=0.0):
    return magic(r=real, i=imag)

#### EOF ######################################################################
