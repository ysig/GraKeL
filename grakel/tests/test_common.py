from sklearn.utils.estimator_checks import check_estimator
from grakel.graph_kernels import GraphKernel

def test_transformer():
    return check_estimator(GraphKernel)

print(test_transformer())
