from sklearn.utils.estimator_checks import check_estimator
from grakel import GraphKernel

def test_transformer():
    return check_estimator(GraphKernel)
