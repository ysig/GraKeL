"""
=========================================================
Classification on MUTAG using a lovasz, svm-theta kernels
=========================================================

An example plot of :class:`grakel.GraphKernel`, :c;ass:`grakel.lovasz_theta`, :class:`grakel.svm_theta`
"""
print(__doc__)
import numpy as np

from time import time

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

# Loads the Mutag dataset from:
# https://ls11-www.cs.tu-dortmund.de/staff/morris/graphkerneldatasets
# the biggest collection of benchmark datasets for graph_kernels.
mutag = datasets.fetch_dataset("MUTAG", verbose=False)
G, y = mutag.data, mutag.target
C_grid = (10. ** np.arange(4,10,1) / len(G)).tolist()

niter = 10
kernel_names = ["lovasz_theta", "svm_theta"]
stats = {k: {"acc": list(), "time": list()} for k in kernel_names}

for i in range(niter):
    # Train-test split of graph data
    G_train, G_test, y_train, y_test = train_test_split(G, y, test_size=0.1)


    for kernel_name in kernel_names:
        start = time()
        # Initialise a weifeiler kernel, with a dirac base_kernel.
        gk = GraphKernel(kernel={"name": kernel_name}, normalize=True)

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

        stats[kernel_name]["acc"].append(acc)
        stats[kernel_name]["time"].append(end-start)

print("Mean values of", niter, "iterations:")
for k in kernel_names:
    print(k, "> Accuracy:", str(round(np.mean(stats[k]["acc"])*100, 2)), "% | Took:",
          sec_to_time(np.mean(stats[k]["time"])))
