"""
=========================================================
Classification statistics on the MUTAG, ENZYMES datasets.
=========================================================

An example plot of :class:`grakel.GraphKernel`
"""
print(__doc__)

import time
import matplotlib.pyplot as plt

from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split
from sklearn import svm

from grakel import datasets
from grakel import GraphKernel


def sec_to_time(sec):
    """Print time in a correct format."""
    dt = list()
    days = int(sec // 86400)
    if days > 0:
        sec -= 86400*days
        dt.append(str(days) + " d")

    hrs = int(sec // 3600)
    if hrs > 0:
        sec -= 3600*hrs
        dt.append(str(hrs) + " h")

    mins = int(sec // 60)
    if mins > 0:
        sec -= 60*mins
        dt.append(str(mins) + " m")

    if sec > 0:
        dt.append(str(round(sec, 2)) + " s")
    return " ".join(dt)


# Loads the MUTAG, ENZYMES dataset from:
# https://ls11-www.cs.tu-dortmund.de/staff/morris/graphkerneldatasets
# the biggest collection of benchmark datasets for graph_kernels.

datasets = ["MUTAG", "MSRC_21C"]

kernels = {
    "Shortest Path": [{"name": "shortest_path"}],
    "Graphlet Sampling": [{"name": "graphlet_sampling",
                           "n_samples": 150}],
    "Weisfeiler-Lehman/Subtree": [{"name": "weisfeiler_lehman", "niter": 5},
                                  {"name": "subtree_wl"}],
    "Weisfeiler-Lehman/Shortest-Path": [{"name": "weisfeiler_lehman",
                                         "niter": 5},
                                        {"name": "shortest_path"}]
}

columns = datasets
rows = sorted(list(kernels.keys()))
data_dataset = list()
for (j, d) in enumerate(columns):
    print(d)
    data_kernel = list()
    dataset_d = datasets.fetch_dataset(d, verbose=False)
    G, y = dataset_d.data, dataset_d.target

    # Train-test split of graph data
    G_train, G_test, y_train, y_test = train_test_split(G, y, test_size=0.1)

    for (i, k) in enumerate(rows):
        print(k, end=" ")
        gk = GraphKernel(kernel=kernels[k], normalize=True)
        print("", end=".")

        # Calculate the kernel matrix.
        start = time.time()
        K_train = gk.fit_transform(G_train)
        K_test = gk.transform(G_test)
        end = time.time()
        print("", end=".")

        # Initialise an SVM and fit.
        clf = svm.SVC(kernel='precomputed')
        clf.fit(K_train, y_train)
        print("", end=". ")

        # Predict and test.
        y_pred = clf.predict(K_test)

        # Calculate accuracy of classification.
        data_kernel.append(
            sec_to_time(round(end - start, 2)) +
            " ~ " + str(round(accuracy_score(y_test, y_pred)*100, 2)) + "%")
        print(data_kernel[-1])
    data_dataset.append(data_kernel)
    print("")


# Print results on a table using pyplot
bbox = [0.45, 0.25, 0.6, 0.6]
table = plt.table(cellText=[list(q) for q in zip(*data_dataset)],
                  rowLabels=rows, colLabels=columns, cellLoc = 'center',
                  rowLoc = 'center', loc='center', bbox=bbox)

_ = plt.axis('off')

plt.show()
