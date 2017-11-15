'''
    graph.py

Module graph contains the definition of the Graph class
and utilities dealing with graph objects.
'''
from __future__ import absolute_import

__LICENSE__     = '''
Copyright (c) 2015 Peter Dobsan

This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 3 of the License, or (at your
option) any later version.  This program is distributed in the hope that
it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
'''

__all__ = [
    'Graph',
    'autgrp',
    'isomorphic',
    'certificate',
    'delete_random_edge',
]

from . import nautywrap
import random


class Graph(object):
    '''
    Graph instantiates an adjacency dictionary based graph object.
    It can represent vertex colored, directed or undirected graphs.
    '''

    def __init__(self, number_of_vertices, directed=False,
                 adjacency_dict={},
                 vertex_coloring=[]):
        '''
        *number_of_vertices*
            The number of vertices of the graph; the vertices are
            labeled from zero.  Mandatory argument.

        *directed*
            Indicate wether the grap is directed or not.  Optional,
            default is False.

        *adjacency_dict*
            key: a vertex, value: a list of vertices linked to the
            key vertex.  Optional, default is an empty dictionary.

        *vertex_coloring*
            A list of disjoint sets of vertices representing a
            partition of the vertex set; vertices not listed are
            placed into a single additional part.  Optional, default
            is no coloring.
        '''
        self.number_of_vertices = number_of_vertices
        self.directed = directed
        self.set_adjacency_dict(adjacency_dict)
        self.set_vertex_coloring(vertex_coloring)

    def _check_vertices(self, vs):
        for v in vs:
            if not (0 <= v and v < self.number_of_vertices):
                raise ValueError(
                'vertex %d conflicts with number_of_vertices=%d' %
                                 (v, self.number_of_vertices))

    def _get_adjacency_dict(self):
        return self._adjacency_dict

    adjacency_dict = property(_get_adjacency_dict)

    def set_adjacency_dict(self, adjacency_dict):
        '''
        Set the adjacency relations of the Graph.

        *adjacency_dict*
            key: a vertex, value: a list of vertices linked to the
            key vertex.
        '''
        for v, vs in adjacency_dict.items():
            self._check_vertices([v])
            self._check_vertices(vs)
        self._adjacency_dict = adjacency_dict.copy()

    def connect_vertex(self, v, neighbors):
        '''
        Connect a vertex to some other vertices.

        *v*
            A vertex of the Graph. The *tail* of the arcs if the Graph
            is directed.
        *neighbors*
            A vertex or a list of vertices to which *v* should be connected.
            The *heads* of the arcs if the Graph is directed.

        '''
        self._check_vertices([v])
        if isinstance(neighbors, list):
            self._check_vertices(neighbors)
            self._adjacency_dict[v] = neighbors
        else:
            self._check_vertices([neighbors])
            self._adjacency_dict.setdefault(v, [])
            self._adjacency_dict[v].append(neighbors)

    def _get_vertex_coloring(self):
        return self._vertex_coloring

    vertex_coloring = property(_get_vertex_coloring)

    def set_vertex_coloring(self, vertex_coloring):
        '''
        Define a vertex coloring of the Graph.

        *vertex_coloring*
            A list of disjoint sets of vertices representing a
            partition of the vertex set; vertices not listed are
            placed into a single additional part.
        '''
        self._vertex_coloring = []
        if vertex_coloring:
            vs = set(range(self.number_of_vertices))
            for p in vertex_coloring:
                if p <= vs:
                    self._vertex_coloring.append(p)
                    vs -= p
                else:
                    raise ValueError('Invalid partition: %s' % vertex_coloring)
            if vs:
                self._vertex_coloring.append(vs)
            if len(self._vertex_coloring) == 1:
                self._vertex_coloring = []

    def __repr__(self):
        s = ['Graph(number_of_vertices=%d, directed=%s,' %
             (self.number_of_vertices, self.directed)]
        s.append(' adjacency_dict = {')
        for k,v in self._adjacency_dict.items():
            v.sort()
            s.append('  %d: %s,' % (k,v))
        s.append(' },')
        s.append(' vertex_coloring = [')
        for x in self._vertex_coloring:
            s.append('  set(%s),' % list(x))
        s.append(' ],')
        s.append(')')
        return '\n'.join(s)


def autgrp(g):
    '''
    Compute the automorphism group of a graph.

    *g*
        A Graph object.

    return -> (generators, grpsize1, grpsize2, orbits, numorbits)
        For the detailed description of the returned components, see
        Nauty's documentation.
    '''
    if not isinstance(g, Graph):
        raise TypeError
    return nautywrap.graph_autgrp(g)


def certificate(g):
    '''
    Compute a certificate based on the canonical labeling of vertices.

    *g*
        A Graph object.

    return ->
        The certificate as a byte string.
    '''
    if not isinstance(g, Graph):
        raise TypeError
    return nautywrap.graph_cert(g)


def isomorphic(a, b):
    '''
    Determine if two graphs are isomorphic.

    *a,b*
        Two Graph objects.

    return ->
        True if *a* and *b* are isomorphic graphs, False otherwise,
    '''
    if a.number_of_vertices != b.number_of_vertices:
        return False
    elif list(map(len, a.vertex_coloring)) != list(map(len, b.vertex_coloring)):
        return False
    else:
        return certificate(a) == certificate(b)


def delete_random_edge(g):
    '''
    Delete a random edge from a graph.

    *g*
        A Graph object.

    return ->
        The deleted edge as a tuple or (None, None) if no edge is left.
    '''
    if g.adjacency_dict:
        # pick a random vertex 'x' which is connected
        x = random.sample(list(g.adjacency_dict),1)[0]
        # remove a random edge connected to 'x'
        xs = g.adjacency_dict[x]
        y = xs.pop(random.randrange(len(xs)))
        if not xs:
            g.adjacency_dict.pop(x)
        # if g is not directed make sure to remove edge completely
        if (not g.directed) and y in g.adjacency_dict:
            ys = g.adjacency_dict[y]
            if x in ys:
                ys.remove(x)
    else:
        # the graph has no edges
        x, y = None, None
    return (x, y)

