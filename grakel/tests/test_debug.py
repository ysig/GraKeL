import pytest
from cvxopt.base import matrix, spmatrix
from cvxopt.solvers import sdp


@pytest.mark.parametrize(
    "nv, ne, e_list, x_list",
    [
        (
            6,
            6,
            [0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 6, 6, 6, 6],
            [1, 1, 2, 2, 8, 8, 10, 10, 15, 15, 16, 16, 0, 7, 14, 21, 28, 35],
        )
    ],
)
def test_windows_sdp(nv, ne, e_list, x_list) -> None:
    # initialise g sparse (to values -1, based on two list that
    # define index and one that defines shape
    g_sparse = spmatrix(-1, x_list, e_list, (nv * nv, ne + 1))

    # Initialise optimization parameters
    h = matrix(-1.0, (nv, nv))
    c = matrix([0.0] * ne + [1.0])

    # Solve the convex optimization problem
    # Should raise here on windows
    sol = sdp(c, Gs=[g_sparse], hs=[h])
    assert sol is not None
