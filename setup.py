from __future__ import print_function
import sys
from setuptools import setup, find_packages

with open('requirements.txt') as f:
    INSTALL_REQUIRES = [l.strip() for l in f.readlines() if l]

try:
    import numpy
except ImportError:
    print('numpy is required during installation')
    sys.exit(1)

try:
    import scipy
except ImportError:
    print('scipy is required during installation')
    sys.exit(1)

try:
    import cvxopt
except ImportError:
    print('cvxopt is required during installation')
    sys.exit(1)

try:
    import pynauty
except ImportError as ie:
    print('Import Error [pynauty]:',ie)
    sys.exit(1)

setup(name='grakel',
      version='0.0.1',
      description='A scikit-learn compatible library for graph kernels',
      author='Ioannis Siglidis [LiX / DaSciM]',
      packages=find_packages(),
      install_requires=INSTALL_REQUIRES,
      author_email='',
      )
