"""The svm theta kernel as defined in :cite:`Johansson2015LearningWS`."""

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
    _metric_type : str, fixed="svm"
        The type of metric calculated from graphs.

    """

    _metric_type = "svm"
