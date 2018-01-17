"""
==============================================================================
Fit-Transform and classification on the MUTAG dataset using a WL-dirac kernel.
==============================================================================

An example plot of :class:`grakel.graph_kernels`
"""
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split
from sklearn import svm

import grakel.dataset.base as dataset
import grakel.graph_kernels as gkl

# Loads the Mutag dataset from:
# https://ls11-www.cs.tu-dortmund.de/staff/morris/graphkerneldatasets
# the biggest collection of benchmark datasets for graph_kernels.
G, C = dataset.load_dataset("MUTAG", verbose=True)

# Train-test split of graph data
GTr, GTe, ytr, yte = train_test_split(G, C, test_size=0.1)

# Initialise a weifeiler kernel, with a dirac base_kernel.
# gk = gkl.GraphKernel(kernel=[
#    {"name": "weisfeiler_lehman", "niter": 5},
#    {"name": "dirac"}], normalize=True)
gk = gkl.GraphKernel(kernel={"name": "shortest_path"}, normalize=True)

# Calculate the kernel matrix.
KTr = gk.fit_transform(GTr)
KTe = gk.transform(GTe)

# Initialise an SVM and fit.
clf = svm.SVC(kernel='precomputed', C=100)
clf.fit(KTr, ytr)

# Predict and test.
y_pred = clf.predict(KTe)

# Calculate accuracy of classification.
acc = accuracy_score(yte, y_pred)

print("Accuracy:", str(round(acc*100, 2)), "%")
