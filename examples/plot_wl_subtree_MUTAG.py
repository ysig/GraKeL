"""
=============================================================
Classification on the MUTAG dataset using a WL-subtree kernel
=============================================================

An example plot of :class:`grakel.GraphKernel`, :class:`grakel.weisfeiler_lehman`, :class:`grakel.subtree_wl`
"""
print(__doc__)

from time import time

from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split
from sklearn import svm

from grakel import datasets
from grakel import GraphKernel

# Loads the Mutag dataset from:
# https://ls11-www.cs.tu-dortmund.de/staff/morris/graphkerneldatasets
# the biggest collection of benchmark datasets for graph_kernels.
mutag = datasets.fetch_dataset("MUTAG", verbose=False)
G, y = mutag.data, mutag.target

# Train-test split of graph data
G_train, G_test, y_train, y_test = train_test_split(G, y, test_size=0.1, random_state=42)

start = time()
# Initialise a weifeiler kernel, with a dirac base_kernel.
gk = GraphKernel(kernel=[{"name": "weisfeiler_lehman", "niter": 5},
                         {"name": "subtree_wl"}], normalize=True)

# Calculate the kernel matrix.
K_train = gk.fit_transform(G_train)
K_test = gk.transform(G_test)
end = time()

# Initialise an SVM and fit.
clf = svm.SVC(kernel='precomputed', C=1)
clf.fit(K_train, y_train)

# Predict and test.
y_pred = clf.predict(K_test)

# Calculate accuracy of classification.
acc = accuracy_score(y_test, y_pred)

print("Accuracy:", str(round(acc*100, 2)), "% | Took:",
      str(round(end - start, 2)), "s")
