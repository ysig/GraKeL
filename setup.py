"""The general `setup.py`  file."""
# Author: Ioannis Siglidis <y.siglidis@gmail.com>
# License: BSD 3 clause
# Author: Ioannis Siglidis <y.siglidis@gmail.com>
# License: BSD 3 clause
from __future__ import print_function
import sys
import importlib
import warnings
from platform import system
from setuptools import setup, find_packages, Extension


# Import/install setup dependencies
def install_and_import(package):
    try:
        importlib.import_module(package)
    except ImportError:
        from pip._internal import main as pip_main
        warnings.warn('package ' + package +
                      ' is required through installation: trying to install it with pip')
        try:
            pip_main(['install', package])
        except Exception:
            raise

    globals()[package] = importlib.import_module(package)


install_and_import('numpy')
from numpy import get_include
install_and_import('Cython')
from Cython.Build import build_ext

# Compile extensions

# Set optimization arguments for compilation
OS = system()
if OS == 'Windows':
    extra_compile_args = ["/O2"]
elif OS in ['Linux', 'Darwin']:
    extra_compile_args = ["-O3"]

# Add the _c_functions extension on kernels
ext_address = "./grakel/kernels/_c_functions/"
ext = Extension(name="grakel.kernels._c_functions",
                sources=[ext_address + "functions.pyx",
                         ext_address + "src/ArashPartov.cpp",
                         ext_address + "src/sm_core.cpp"],
                include_dirs=[ext_address + "include", get_include()],
                depends=[ext_address + "include/functions.hpp"],
                language="c++",
                extra_compile_args=extra_compile_args)

# Add the bliss library extension for calculating isomorphism
isodir = "./grakel/kernels/_isomorphism/"
blissdir = isodir + 'bliss-0.50/'

# The essential bliss source files
blisssrcs = ['graph.cc', 'heap.cc', 'orbit.cc', 'partition.cc', 'uintseqhash.cc']
blisssrcs = [blissdir + src for src in blisssrcs]
pn = str(sys.version_info[0])

# Compile intpybliss
intpybliss = Extension(name="grakel.kernels._isomorphism.intpybliss",
                       define_macros=[('MAJOR_VERSION', '0'),
                                      ('MINOR_VERSION', '50beta')],
                       include_dirs=[blissdir],
                       language="c++",
                       sources=[isodir + 'intpyblissmodule_' + pn + '.cc']+blisssrcs
                       )

# Make bliss extension
bliss = Extension(name="grakel.kernels._isomorphism.bliss",
                  include_dirs=[isodir],
                  language="c++",
                  sources=[isodir + 'bliss.pyx']
                  )

# Add readme pypi
with open("README.md", "r") as fh:
    long_description = fh.read()
    long_description = '\n'.join(s for i, s in enumerate(
            [s for s in long_description.split('\n')
             if not (len(s) >= 2 and s[:2] == "[!")]) if i != 2)

# Package requierements
with open('requirements.txt') as f:
    INSTALL_REQUIRES = [l.strip() for l in f.readlines() if l]


setup(name='grakel-dev',
      version='0.1a5',
      description='A scikit-learn compatible library for graph kernels',
      long_description=long_description,
      long_description_content_type='text/markdown',
      project_urls={
        'Documentation': 'https://ysig.github.io/GraKeL/dev/',
        'Send us Feedback!': 'http://www.lix.polytechnique.fr/dascim/contact/',
        'Source': 'https://github.com/ysig/GraKeL/tree/develop',
        'Tracker': 'https://github.com/ysig/GraKeL/issues',
        },
      author='Ioannis Siglidis [LiX / DaSciM]',
      author_email='y.siglidis@gmail.com',
      url='https://ysig.github.io/GraKeL/dev/',
      license="BSD",
      classifiers=['Intended Audience :: Science/Research',
                   'Intended Audience :: Developers',
                   'License :: OSI Approved',
                   'Programming Language :: C',
                   'Programming Language :: Python',
                   'Topic :: Software Development',
                   'Topic :: Scientific/Engineering',
                   'Operating System :: POSIX',
                   'Operating System :: Unix',
                   'Operating System :: MacOS',
                   'Programming Language :: Python :: 2.7',
                   'Programming Language :: Python :: 3.5',
                   'Programming Language :: Python :: 3.6',
                   'Programming Language :: Python :: 3.7',
                   ],
      python_requires='>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, !=3.4.*, <4',
      packages=find_packages(),
      package_data={'grakel.tests': ['data/Cuneiform/*.txt', 'data/MUTAG/*.txt']},
      install_requires=INSTALL_REQUIRES,
      extras_require={
        'lovasz': ["cvxopt>=1.2.0"]
      },
      ext_modules=[intpybliss, bliss, ext],
      cmdclass={'build_ext': build_ext},
      )
