"""Import datasets related with graph kernels, from a large collection."""
from grakel.datasets.base import fetch_dataset
from grakel.datasets.base import get_dataset_info
from grakel.datasets.testing import generate_dataset

__all__ = [
    "get_dataset_info",
    "fetch_dataset",
    "generate_dataset"
]
