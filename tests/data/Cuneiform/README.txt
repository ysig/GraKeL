README for dataset Cuneiform


=== Usage ===

This folder contains the following comma separated text files 
(replace DS by the name of the dataset):

n = total number of nodes
m = total number of edges
N = number of graphs

(1) 	DS_A.txt (m lines) 
	sparse (block diagonal) adjacency matrix for all graphs,
	each line corresponds to (row, col) resp. (node_id, node_id)

(2) 	DS_graph_indicator.txt (n lines)
	column vector of graph identifiers for all nodes of all graphs,
	the value in the i-th line is the graph_id of the node with node_id i

(3) 	DS_graph_labels.txt (N lines) 
	class labels for all graphs in the dataset,
	the value in the i-th line is the class label of the graph with graph_id i

(4) 	DS_node_labels.txt (n lines)
	column vector of node labels,
	the value in the i-th line corresponds to the node with node_id i

There are OPTIONAL files if the respective information is available:

(5) 	DS_edge_labels.txt (m lines; same size as DS_A_sparse.txt)
	labels for the edges in DS_A_sparse.txt 

(6) 	DS_edge_attributes.txt (m lines; same size as DS_A.txt)
	attributes for the edges in DS_A.txt 

(7) 	DS_node_attributes.txt (n lines) 
	matrix of node attributes,
	the comma seperated values in the i-th line is the attribute vector of the node with node_id i

(8) 	DS_graph_attributes.txt (N lines) 
	regression values for all graphs in the dataset,
	the value in the i-th line is the attribute of the graph with graph_id i


=== Description ===

The Cuneiform dataset contains graphs representing 29 different Hittite cuneiform signs.
The data was obtained from nine cuneiform tablets written by scholars of Hittitology in
the course of a study about individualistic characteristics of cuneiform hand writing.
After automated extraction of individual wedges, the affiliation of the wedges to the 
cuneiform signs were determined manually. The graph model is explained in detail in the
referenced publication.


=== References ===

Nils M. Kriege, Matthias Fey, Denis Fisseler, Petra Mutzel, Frank Weichert
Recognizing Cuneiform Signs Using Graph Based Methods. 2018. arXiv:1802.05908
https://arxiv.org/abs/1802.05908


=== Description of Labels ===

Node labels were converted to integer values using this map:

Component 0:
	0	depthPoint
	1	tailVertex
	2	leftVertex
	3	rightVertex

Component 1:
	0	vertical
	1	Winkelhaken
	2	horizontal



Edge labels were converted to integer values using this map:

Component 0:
	0	wedge
	1	arrangement



Class labels were converted to integer values using this map:

	0	tu
	1	ta
	2	ti
	3	nu
	4	na
	5	ni
	6	bu
	7	ba
	8	bi
	9	zu
	10	za
	11	zi
	12	su
	13	sa
	14	si
	15	hu
	16	ha
	17	hi
	18	du
	19	da
	20	di
	21	ru
	22	ra
	23	ri
	24	ku
	25	ka
	26	ki
	27	lu
	28	la
	29	li
