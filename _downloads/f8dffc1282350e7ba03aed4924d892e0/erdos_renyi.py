"""
===========================================================================
Graph classification on a randomly generated dataset of Erdos-Renyi graphs.
===========================================================================

Script makes use of :class:`grakel.Graph` and :class:`grakel.ShortestPath`
"""
from __future__ import print_function
print(__doc__)

import numpy as np

from random import random

from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score

from grakel import Graph
from grakel.kernels import ShortestPath

# Generates 3 sets of Erdos-Renyi graphs. Each edge is included in the graph with probability p
# independent from every other edge. The probability p is set equal to 0.25, 0.5 and 0.75 for 
# the graphs of the 1st, 2nd and 3rd set, respectivery
Gs = list()
y = list()
probs = [0.25, 0.5, 0.75]
for i in range(len(probs)):
	for j in range(5, 35):
		edges = list()
		for n1 in range(j):
			for n2 in range(n1+1, j):
				if random() <= probs[i]:
					edges.append((n1, n2))
					edges.append((n2, n1))

		Gs.append(Graph(edges))
		y.append(i)

# Splits the dataset into a training and a test set
G_train, G_test, y_train, y_test = train_test_split(Gs, y, test_size=0.1, random_state=42)

# Uses the shortest path kernel to generate the kernel matrices
gk = ShortestPath(normalize=True, with_labels=False)
K_train = gk.fit_transform(G_train)
K_test = gk.transform(G_test)

# Uses the SVM classifier to perform classification
clf = SVC(kernel="precomputed")
clf.fit(K_train, y_train)
y_pred = clf.predict(K_test)

# Computes and prints the classification accuracy
acc = accuracy_score(y_test, y_pred)
print("Accuracy:", str(round(acc*100, 2)) + "%")
