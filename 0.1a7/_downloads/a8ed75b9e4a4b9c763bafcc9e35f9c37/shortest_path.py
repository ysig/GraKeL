"""
=============================================================
Graph classification on MUTAG using the shortest path kernel.
=============================================================

Script makes use of :class:`grakel.ShortestPath`
"""
from __future__ import print_function
print(__doc__)

import numpy as np

from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score

from grakel.datasets import fetch_dataset
from grakel.kernels import ShortestPath

# Loads the MUTAG dataset
MUTAG = fetch_dataset("MUTAG", verbose=False)
G, y = MUTAG.data, MUTAG.target

# Splits the dataset into a training and a test set
G_train, G_test, y_train, y_test = train_test_split(G, y, test_size=0.1, random_state=42)

# Uses the shortest path kernel to generate the kernel matrices
gk = ShortestPath(normalize=True)
K_train = gk.fit_transform(G_train)
K_test = gk.transform(G_test)

# Uses the SVM classifier to perform classification
clf = SVC(kernel="precomputed")
clf.fit(K_train, y_train)
y_pred = clf.predict(K_test)

# Computes and prints the classification accuracy
acc = accuracy_score(y_test, y_pred)
print("Accuracy:", str(round(acc*100, 2)) + "%")
