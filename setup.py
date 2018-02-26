"""The general `setup.py`  file."""
from __future__ import print_function
import sys
from setuptools import setup, find_packages, Extension
from Cython.Distutils import build_ext

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
    print('Import Error [pynauty]:', ie)
    sys.exit(1)

# Add the _c_functions extension on kernels
ext_address = "./grakel/kernels/_c_functions/"
ext = Extension(name="grakel.kernels._c_functions",
                sources=[ext_address + "functions.pyx",
                         ext_address + "src/ArashPartov.cpp",
                         ext_address + "src/sm_core.cpp"],
                include_dirs=[ext_address + "include", numpy.get_include()],
                depends=[ext_address + "include/functions.hpp"],
                language="c++",
                extra_compile_args=["-O3", "-std=c++11"])

setup(name='grakel',
      version='0.0.1',
      description='A scikit-learn compatible library for graph kernels',
      author='Ioannis Siglidis [LiX / DaSciM]',
      packages=find_packages(),
      install_requires=INSTALL_REQUIRES,
      author_email='ioannis.siglidis@inria.fr',
      ext_modules=[ext],
      cmdclass={'build_ext': build_ext}
      )
