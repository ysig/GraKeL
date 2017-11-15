'''
pynauty  --  isomorphism testing and automorphism groups of graphs

A Python extension module to Brendan McKay's "nauty" C procedures
for determining the automorphism group of a vertex colored graph
and producing canonical labeling for isomorphism testing.

Classes:

    Graph   - An adjacency dictionary based graph object.
        Graph can represent vertex colored, directed or undirected graphs.

Functions:

    autgrp      - Compute the automorphism group of a graph.
    isomorphic  - Compare two graphs for isomorphism.
    certificate - Compute a "certificate" based on the canonical labeling
                  of the graph's vertices.
'''
from __future__ import absolute_import

__version__ = '0.6'

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
from .graph import *

del graph
del nautywrap

