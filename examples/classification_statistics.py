"""
==============================================================================
Classification statistics on the MUTAG, ENZYMES datasets.
==============================================================================

An example plot of :class:`grakel.graph_kernels`
"""
import time

import numpy as np
import matplotlib.pyplot as plt

from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split
from sklearn import svm

import grakel.dataset.base as dataset
import grakel.graph_kernels as gkl


def sec_to_time(sec):
    """Print time in a correct format."""
    dt = list()
    days = sec // 86400
    if days > 0:
        sec -= 86400*days
        dt.append(str(days) + " d")

    hrs = sec // 3600
    if hrs > 0:
        sec -= 3600*hrs
        dt.append(str(hrs) + " h")

    mins = sec // 60
    if mins > 0:
        sec -= 60*mins
        dt.append(str(mins) + " m")

    if sec > 0:
        dt.append(str(sec) + " s")
    return " ".join(dt)

# Loads the MUTAG, ENZYMES dataset from:
# https://ls11-www.cs.tu-dortmund.de/staff/morris/graphkerneldatasets
# the biggest collection of benchmark datasets for graph_kernels.


datasets = ["MUTAG", "ENZYMES"]
kernels = {
    "Shortest Path": [{"name": "shortest_path"}],
    "Graphlet Sampling": [{"name": "graphlet_sampling",
                          "n_samples": 150}],
    "Random Walk": [{"name": "random_walk", "lambda": 10**(-3)}],
    "Weisfeiler-Lehman/Subtree": [{"name": "weisfeiler_lehman", "niter": 5},
                                  {"name": "subtree_wl"}],
    "Weisfeiler-Lehman/Shortest-Path": [{"name": "weisfeiler_lehman",
                                         "niter": 5},
                                        {"name": "shortest_path"}]

}

columns = datasets
rows = sorted(list(kernels.keys()))
data_acc = np.empty(shape=(len(rows), len(columns)))
data_time = np.empty(shape=(len(rows), len(columns)))
for (j, d) in enumerate(columns):
    print(d)
    G, C = dataset.load_dataset(d, verbose=False)

    # Train-test split of graph data
    GTr, GTe, ytr, yte = train_test_split(G, C, test_size=0.1)

    for (i, k) in enumerate(rows):
        print(k, end=" ")
        gk = gkl.GraphKernel(kernel=kernels[k], normalize=True, concurrency=-1)
        print("", end=".")
        # Calculate the kernel matrix.
        start = time.time()
        KTr = gk.fit_transform(GTr)
        KTe = gk.transform(GTe)
        end = time.time()
        print("", end=".")
        # Initialise an SVM and fit.
        clf = svm.SVC(kernel='precomputed')
        clf.fit(KTr, ytr)
        print("", end=". ")

        # Predict and test.
        y_pred = clf.predict(KTe)

        # Calculate accuracy of classification.
        data_acc[i, j] = round(accuracy_score(yte, y_pred)*100, 2)
        data_time[i, j] = round(end - start, 2)
        print(sec_to_time(data_time[i, j]))
    print("")

acc_table = plt.table(cellText=data_acc, rowLabels=rows, colLabels=columns,
                      loc='left')
time_table = plt.table(cellText=data_time, rowLabels=rows, colLabels=columns,
                       loc='right')
plt.show()
