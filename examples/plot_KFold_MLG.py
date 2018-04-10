"""
=========================================================================================
KFold classification of a dataset [Cuneiform] using the fast Multiscale Laplacian kernel.
=========================================================================================

An example plot of :class:`grakel.GraphKernel`, :class:`grakel.multiscale_laplacian_fast`
"""
print(__doc__)

import argparse

from grakel import GraphKernel
from grakel import datasets

# Create an argument parser for the installer of pynauty
parser = argparse.ArgumentParser(
    description='Measuring classification accuracy '
                ' on multiscale_laplacian_fast')

parser.add_argument(
    '--dataset',
    help='chose the datset you want the tests to be executed',
    type=str,
    default="Cuneiform",
)

# Get the dataset name
dataset_name = parser.parse_args().dataset

# Check the dataset provided by the user
dinfo = datasets.get_dataset_info(dataset_name)
if dinfo is None:
    raise TypeError('Dataset not found!')
elif not dinfo["nl"]:
    raise TypeError('Dataset must have contain node attributes.')


# The baseline dataset for node/edge-attributes
dataset_attr = datasets.fetch_dataset(dataset_name,
                                      with_classes=True,
                                      prefer_attr_nodes=True,
                                      prefer_attr_edges=True,
                                      verbose=True)

import numpy as np

from tqdm import tqdm
from time import time

from sklearn.metrics import accuracy_score
from sklearn.model_selection import KFold
from sklearn import svm

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

# Loads the Mutag dataset from:
# https://ls11-www.cs.tu-dortmund.de/staff/morris/graphkerneldatasets
# the biggest collection of benchmark datasets for graph_kernels.
G, y = dataset_attr.data, dataset_attr.target
C_grid = (10. ** np.arange(4, 10, 1) / len(G)).tolist()

stats = {"acc": list(), "time": list()}
kf = KFold(n_splits=10, random_state=42, shuffle=True)
niter = kf.get_n_splits(y)

for train_index, test_index in tqdm(kf.split(G, y),
                                    total=niter):
    # Train-test split of graph data
    tri = train_index.tolist()
    tei = test_index.tolist()

    G_train, G_test = list(), list()
    y_train, y_test = list(), list()
    for (i, (g, t)) in enumerate(zip(G, y)):
        if len(tri) and i == tri[0]:
            G_train.append(g)
            y_train.append(t)
            tri.pop(0)
        elif len(tei) and i == tei[0]:
            G_test.append(g)
            y_test.append(t)
            tei.pop(0)

    start = time()
    gk = GraphKernel(kernel={"name": "multiscale_laplacian", "which": "fast"})

    # Calculate the kernel matrix.
    K_train = gk.fit_transform(G_train)
    K_test = gk.transform(G_test)
    end = time()

    # Cross validation on C, variable
    acc = 0
    for c in C_grid:
        # Initialise an SVM and fit.
        clf = svm.SVC(kernel='precomputed', C=c)

        # Fit on the train Kernel
        clf.fit(K_train, y_train)

        # Predict and test.
        y_pred = clf.predict(K_test)

        # Calculate accuracy of classification.
        acc = max(acc, accuracy_score(y_test, y_pred))

    stats["acc"].append(acc)
    stats["time"].append(end-start)

print("Mean values of ", niter, "folds:")
print("MLG [Fast] > Accuracy:",
      str(round(np.mean(stats["acc"])*100, 2)),
      "% | Took:", sec_to_time(np.mean(stats["time"])))
