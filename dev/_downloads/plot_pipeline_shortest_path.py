"""
==============================================================================
Cross-Validation scikit-learn Pipeline example on the MUTAG dataset using the Shortest Path Kernel.
==============================================================================

An example plot of :class:`grakel.GraphKernel`, :class:`grakel.ShortestPath`
"""
from __future__ import print_function
print(__doc__)

import numpy as np

from time import time

from sklearn import svm
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.metrics import accuracy_score

from grakel import datasets
from grakel import GraphKernel

# Loads the Mutag dataset from:
# https://ls11-www.cs.tu-dortmund.de/staff/morris/graphkerneldatasets
# the biggest collection of benchmark datasets for graph_kernels.
mutag = datasets.fetch_dataset("MUTAG", verbose=False)
G, y = mutag.data, mutag.target
C_grid = (10. ** np.arange(1,10,1) / len(G)).tolist()
n_folds = 10

estimator = make_pipeline(
    GraphKernel(kernel=dict(name="shortest_path"), normalize=True),
    GridSearchCV(svm.SVC(kernel='precomputed'), dict(C=C_grid),
                 scoring='accuracy'))

acc = accuracy_score(y, cross_val_predict(estimator, G, y, cv=n_folds))
print("Accuracy:", str(round(acc*100, 2)) + "%")
