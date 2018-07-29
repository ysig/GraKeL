"""Main .pyx file for the _isomorphism subpackage of grakel.kernels"""
# Author: Ioannis Siglidis <y.siglidis@gmail.com>
# This file is a modification and extension of the [GNU LPGL] licensed 
# PyBliss which can be found at: http://www.tcs.hut.fi/Software/bliss/
# PyBliss and Bliss are copyright of their respective owners.
# License: BSD 3 clause"
import cython
import sys
import types
import collections

from grakel.kernels._isomorphism import intpybliss
from functools import total_ordering
from six import iteritems

def _report(perm, args):
    [reporter_func, reporter_args, bliss_map_inv] = args;
    if reporter_func:
        p = {}
        for (index,image) in enumerate(perm):
            p[bliss_map_inv[index]] = bliss_map_inv[image]
        reporter_func(p, reporter_args)


class Graph:
    @total_ordering
    class _Vertex:
        def __init__(self, name, color):
            assert type(color) is int
            assert 0 <= color and color < pow(2,31)
            assert isinstance(name, collections.Hashable)
            assert type(name) is int
            self.edges = set()
            self.name = name
            self.color = color
        
        def __hash__(self):
            return self.name.__hash__()

        def __eq__(self, other):
            return str(self.name) == str(other.name)
    
        def __lt__(self, other):
            return str(self.name) < str(other.name)

    """
    The class for undirected graphs.
    """
    def __init__(self, *args):
        """
        Create a new empty graph.
        """
        self._vertices = {}

        if len(args) > 0:
            assert len(args) == 2
            n_nodes, edges = args
            # Initialise vertices
            for v in range(n_nodes):
                self.add_vertex(v)

            # Initialise edges
            for u, v in edges:
                self.add_edge(u, v)

    def __hash__(self):
        return hash(str(self))

    def __str__(self):
        vnames = sorted(self._vertices.keys())
        s = "<"
        for name in vnames:
            vertex = self._vertices[name]
            nnames = sorted(neighbour.name for neighbour in vertex.edges)
            s += "("+str(name)+" "+str(vertex.color)+" "+str(nnames)+")"
        s += ">"
        return s

    def __eq__(self, other):
        """
        Check whether the graph is equal to the graph 'other'.
        That is, whether the graphs have the same vertices colored with
        the same colors and the same edges.
        """
        return self.is_equal(other)

    def __ne__(self, other):
        """
        Check whether the graph is not equal to the graph 'other'.
        """
        return not self.is_equal(other)

    def __le__(self, other):
        """
        Ordering of graphs is not defined, will raise an Error.
        """
        raise RuntimeError("Operation not implemented")

    def __lt__(self, other):
        """
        Ordering of graphs is not defined, will raise an Error.
        """
        raise RuntimeError("Operation not implemented")

    def __ge__(self, other):
        """
        Ordering of graphs is not defined, will raise an Error.
        """
        raise RuntimeError("Operation not implemented")

    def __gt__(self, other):
        """
        Ordering of graphs is not defined, will raise an Error.
        """
        raise RuntimeError("Operation not implemented")

    def __cmp__(self, other):
        """
        Ordering of graphs is not defined, will raise an Error.
        """
        raise RuntimeError("Operation not implemented")

    def nof_vertices(self):
        """
        Get the number of vertices in the graph.
        """
        return len(self._vertices)

    def get_vertices(self):
        """
        Return the (list of) vertices in the graph.
        """
        return self._vertices.keys()

    def add_vertex(self, v, color = 0):
        """
        Add a new vertex named 'v' with color 'c' in the graph.
        If the vertex already exists, do nothing.
        The color must be an integer between 0 and 2^31-1.
        The vertex name v should be hashable and immutable,
        preferably an integer or a string not containing spaces.
        Returns False if the vertex already existed, True otherwise.
        """
        assert type(color) is int
        assert 0 <= color and color < pow(2,31)
        if v in self._vertices:
            return False
        self._vertices[v] = self._Vertex(v, color)
        return True
        
    def del_vertex(self, v):
        """
        Delete the vertex named 'v' from the graph.
        If the vertex does not exist, do nothing.
        The vertex name v should be hashable and immutable,
        preferably an integer or a string not containing spaces.
        """
        if v in self._vertices:
            vertex = self._vertices[v]
            for neighbour in vertex.edges:
                if not(neighbour is vertex):
                    neighbour.edges.remove(vertex)
            del self._vertices[v]

    def add_edge(self, v1, v2):
        """
        Add an edge between vertices v1 and v2.
        Adds the vertices v1 and v2 if they do not already exist.
        The vertex names v1 and v2 should be hashable and immutable,
        preferably integers or strings not containing spaces.
        """
        if v1 not in self._vertices:
            self.add_vertex(v1)
        if v2 not in self._vertices:
            self.add_vertex(v2)
        self._vertices[v1].edges.add(self._vertices[v2])

    def del_edge(self, v1, v2):
        """
        Delete the edge between vertices v1 and v2.
        The vertex names v1 and v2 should be hashable and immutable,
        preferably integers or strings not containing spaces.
        """
        if v1 not in self._vertices:
            return
        if v2 not in self._vertices:
            return
        self._vertices[v1].edges.discard(self._vertices[v2])

    def write_dot(self, file):
        """
        Write the graph into a file in the graphviz dot format.
        """
        file.write("graph g {\n")
        for v,vertex in iteritems(self._vertices):
            file.write("\""+str(v)+"\" [label="+str(vertex.color)+"];\n")
        for v,vertex in iteritems(self._vertices):
            for neighbour in vertex.edges:
                file.write("\""+str(v)+"\" -- \""+str(neighbour.name)+"\";\n")
        file.write("}\n")

    def _make_bliss_graph(self):
        g = intpybliss.create()
        bliss_map = {}
        bliss_map_inv = {}
        for v,vertex in iteritems(self._vertices):
            bliss_map[v] = intpybliss.add_vertex(g, vertex.color)
            bliss_map_inv[bliss_map[v]] = v
        for name,vertex in iteritems(self._vertices):
            for neighbour in vertex.edges:
                intpybliss.add_edge(g,
                                    bliss_map[name],
                                    bliss_map[neighbour.name])
        return (g,bliss_map,bliss_map_inv)

    def find_automorphisms(self,
                           reporter_function,
                           reporter_function_arg = None):
        """
        Find a generating set for the automorphism group of the graph.
        The generators are reported by calling the hook function
        'reporter_function' with two arguments: the generator and
        the user-defined argument 'reporter_function_arg'.
        """
        if not((reporter_function is None) or
               isinstance(reporter_function, types.FunctionType)):
            raise TypeError("the 'reporter_function' argument of "
                            "canonical_labeling() should be None or "
                            "a function")
        (g,bliss_map,bliss_map_inv) = self._make_bliss_graph()
        report_args = [reporter_function,reporter_function_arg, bliss_map_inv]
        intpybliss.find_automorphisms(g, _report, report_args)

    def canonical_labeling(self,
                           reporter_function = None,
                           reporter_function_arg = None):
        """
        Find a canonical labeling for the graph.
        The labeling is a mapping from the vertices to the set {0,...,N-1},
        where N is the number of vertices in the graph.
        The canonical form of the graph can be found by applying
        the 'relabel' method with the returned labeling.
        """
        if not((reporter_function is None) or
               isinstance(reporter_function, types.FunctionType)):
            raise TypeError("the 'reporter_function' argument of "
                            "canonical_labeling() should be None or "
                            "a function")
        (g,bliss_map,bliss_map_inv) = self._make_bliss_graph()
        report_args = [reporter_function,reporter_function_arg, bliss_map_inv]
        cl = intpybliss.canonical_form(g, _report, report_args)
        assert type(cl) is list
        canlab = {}
        for (index,image) in enumerate(cl):
            canlab[bliss_map_inv[index]] = image
        return canlab

    def relabel(self, lab):
        """
        Apply the argument labeling 'lab' to the graph, returning a copy of
        the graph where the vertices have been relabeled according to 'lab'.
        The labeling is given as a dictionary mapping each vertex
        in the graph into some distinct immutable name.
        """
        assert type(lab) is dict
        assert len(lab) == self.nof_vertices()
        g2 = Graph()
        for name,vertex in iteritems(self._vertices):
            inserted = g2.add_vertex(lab[name],vertex.color)
            if not inserted:
                raise RuntimeError("'lab' is not a bijection")
        for name,vertex in iteritems(self._vertices):
            for neighbour in vertex.edges:
                g2.add_edge(lab[name], lab[neighbour.name])
        return g2

    def copy(self):
        """
        Return a copy of this graph.
        """
        g2 = Graph()
        for name,vertex in iteritems(self._vertices):
            g2.add_vertex(name, vertex.color)
        for name,vertex in iteritems(self._vertices):
            for neighbour in vertex.edges:
                g2.add_edge(name, neighbour.name)
        return g2

    def is_equal(self, g):
        """
        Check whether the graph is equal to the graph 'g'.
        Returns True if this is the case and False otherwise.
        """
        if len(self._vertices) != len(g._vertices):
            return False
        for name,vertex in iteritems(self._vertices):
            if name not in g._vertices:
                return False
            if vertex.color != g._vertices[name].color:
                return False
            my_edges = set([v.name for v in vertex.edges])
            g_edges = set([v.name for v in g._vertices[name].edges])
            if my_edges != g_edges:
                return False
        return True

    def get_isomorphism(self, g):
        """
        Return an isomorphism from this graph to the graph 'g' or
        None is thegraphs are not isomorphic.
        Observe that for empty graphs an empty dictionary is returned and
        thus one should not check isomorphism by 'if g1.get_isomorphism(g2)'
        but by 'if g1.get_isomorphism(g2) != None'.
        The returned isomorphism 'lab' is a labeling such that
        self.relabel(lab).is_equal(g) == True.
        """
        if len(self._vertices) != len(g._vertices):
            return None
        # Simple color and degree check
        my_colors = {}
        g_colors = {}
        my_degrees = {}
        g_degrees = {}
        for name,vertex in iteritems(self._vertices):
            if vertex.color in my_colors: my_colors[vertex.color] += 1
            else: my_colors[vertex.color] = 1;
            degree = len(vertex.edges)
            if degree in my_degrees: my_degrees[degree] += 1
            else: my_degrees[degree] = 1;
        for name,vertex in iteritems(g._vertices):
            if vertex.color in g_colors: g_colors[vertex.color] += 1
            else: g_colors[vertex.color] = 1;
            degree = len(vertex.edges)
            if degree in g_degrees: g_degrees[degree] += 1
            else: g_degrees[degree] = 1;
        if my_colors != g_colors:
            return None
        if my_degrees != g_degrees:
            return None
        my_canlab = self.canonical_labeling()
        my_canform = self.relabel(my_canlab)
        g_canlab = g.canonical_labeling()
        g_canform = g.relabel(g_canlab)
        if not my_canform.is_equal(g_canform):
            return None
        g_canlab_inv = {}
        for name,image in iteritems(g_canlab):
            g_canlab_inv[image] = name
        isomorphism = {}
        for name,image in iteritems(my_canlab):
            isomorphism[name] = g_canlab_inv[image]
        assert self.relabel(isomorphism).is_equal(g) == True

        return isomorphism

    def isomorphic(self, g):
        return self.get_isomorphism(g) is not None

