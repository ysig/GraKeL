"""
=========================================================
Example of transforming NetworkX graphs to GraKeL graphs.
=========================================================
"""
from __future__ import print_function
print(__doc__)

import numpy as np
import networkx as nx

from grakel.utils import graph_from_networkx

# Creates a list of two simple graphs
G1 = nx.Graph()
G1.add_nodes_from([0,1,2])
G1.add_edges_from([(0,1), (1,2)])

G2 = nx.Graph()
G2.add_nodes_from([0,1,2])
G2.add_edges_from([(0,1), (0,2), (1,2)])

G_nx = [G1, G2]

# Transforms list of NetworkX graphs into a list of GraKeL graphs
G = graph_from_networkx(G_nx)
print("1 - Simple graphs transformed\n")


# Creates a list of two node-labeled graphs
G1 = nx.Graph()
G1.add_nodes_from([0,1,2])
G1.add_edges_from([(0,1), (1,2)])
nx.set_node_attributes(G1, {0:'a', 1:'b', 2:'a'}, 'label')

G2 = nx.Graph()
G2.add_nodes_from([0,1,2])
G2.add_edges_from([(0,1), (0,2), (1,2)])
nx.set_node_attributes(G2, {0:'a', 1:'b', 2:'c'}, 'label')

G_nx = [G1, G2]

# Transforms list of NetworkX graphs into a list of GraKeL graphs
G = graph_from_networkx(G_nx, node_labels_tag='label')
print("2 - Node-labeled graphs transformed\n")


# Creates a list of two node-attributed graphs
G1 = nx.Graph()
G1.add_nodes_from([0,1,2])
G1.add_edges_from([(0,1), (1,2)])
nx.set_node_attributes(G1, {0:np.array([1.1, 0.8]), 
	1:np.array([0.2, -0.3]), 2:np.array([0.9, 1.0])}, 'attributes')

G2 = nx.Graph()
G2.add_nodes_from([0,1,2])
G2.add_edges_from([(0,1), (0,2), (1,2)])
nx.set_node_attributes(G2, {0:np.array([1.8, 0.5]), 
	1:np.array([-0.1, 0.2]), 2:np.array([2.3, 1.2])}, 'attributes')

G_nx = [G1, G2]

# Transforms list of NetworkX graphs into a list of GraKeL graphs
G = graph_from_networkx(G_nx, node_labels_tag='attributes')
print("3 - Node-attributed graphs transformed")