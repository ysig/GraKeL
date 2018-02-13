"""The svm theta kernel as defined in :cite:`Johansson2015LearningWS`."""
import numpy as np

from sklearn.svm import OneClassSVM
from grakel.kernels import lovasz_theta


class svm_theta(lovasz_theta):
    """Calculate the SVM theta kernel.

    See :cite:`Johansson2015LearningWS`.

    Parameters
    ----------
    X,Y : *valid-graph-format*
        The pair of graphs on which the kernel is applied.

    n_samples : int, default=50
        Number of samples.

    subsets_size_range : tuple, len=2, default=(2,8)
        (min, max) size of the vertex set of sampled subgraphs.

    metric : function (number, number -> number), default=:math:`f(x,y)=x*y`
        The applied metric between the svm_theta numbers of the two graphs.

    Attributes
    ----------
    _internal_metric : function
        A metric function for adjacency matrices.
        Inside this class this metric is the svm_theta.

    """

    def _internal_metric(self, A):
        """Calculate the internal metric.

        Parameters
        ----------
        A : np.array, ndim=2
            The adjacency matrix.

        Returns
        -------
        metric: Number
            Returns the metric number.

        """
        return _calculate_svm_theta_(A)


def _calculate_svm_theta_(A):
    """Calculate the svm theta for the given graph.

    Parameters
    ----------
    A: np.array, ndim=2
        A square numpy array corresponding to the adjacency matrix.

    Returns
    -------
    svm_theta: float
        Returns the svm theta number.

    """
    K = A > 0
    np.fill_diagonal(K, False)
    K.astype(int)

    svm = OneClassSVM(kernel="precomputed")
    svm.fit(K)

    return np.sum(np.abs(svm.dual_coef_[0]), axis=0)
