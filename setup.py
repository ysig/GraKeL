"""The general `setup.py`  file."""
from __future__ import print_function
import sys
from warnings import warn
from platform import system
from setuptools import setup, find_packages, Extension

if 'setuptools' in sys.modules:
    try:
        from setuptools.command.build_ext import build_ext as _build_ext
    except ImportError:
        # We may be in the process of importing setuptools, which tries
        # to import this.
        from distutils.command.build_ext import build_ext as _build_ext
else:
    from distutils.command.build_ext import build_ext as _build_ext

with open('requirements.txt') as f:
    INSTALL_REQUIRES = [l.strip() for l in f.readlines() if l]

OS = system()

if OS == 'Windows':
    extra_compile_args = ["/O2"]
elif OS in ['Linux', 'Darwin']:
    extra_compile_args = ["-O3"]

try:
    import numpy
except ImportError:
    raise ImportError('numpy is required during installation')

try:
    from Cython.Build import build_ext
except ImportError:
    raise ImportError('build_ext from Cython.Build is required during installation')

# Add the _c_functions extension on kernels
ext_address = "./grakel/kernels/_c_functions/"
ext = Extension(name="grakel.kernels._c_functions",
                sources=[ext_address + "functions.pyx",
                         ext_address + "src/ArashPartov.cpp",
                         ext_address + "src/sm_core.cpp"],
                include_dirs=[ext_address + "include", numpy.get_include()],
                depends=[ext_address + "include/functions.hpp"],
                language="c++",
                extra_compile_args=extra_compile_args)

# Add the bliss library extension for calculating isomorphism
isodir = "./grakel/kernels/_isomorphism/"
blissdir = isodir + 'bliss-0.50/'
# The essential bliss source files
blisssrcs = ['graph.cc','heap.cc','orbit.cc','partition.cc','uintseqhash.cc']
blisssrcs = [blissdir + src for src in blisssrcs]
pn = str(sys.version_info[0])

intpybliss = Extension(name="grakel.kernels._isomorphism.intpybliss",
                  define_macros = [('MAJOR_VERSION', '0'),
                                   ('MINOR_VERSION', '50beta')],
                  include_dirs = [blissdir],
                  language="c++",
                  sources = [isodir + 'intpyblissmodule_' + pn + '.cc']+blisssrcs
                  )

bliss = Extension(name="grakel.kernels._isomorphism.bliss",
                  include_dirs = [isodir],
                  language="c++",
                  sources = [isodir + 'bliss.pyx']
                  )

setup(name='grakel-dev',
      version='0.1a3',
      description='A scikit-learn compatible library for graph kernels',
      long_description='A scikit-learn compatible library for graph kernels',
      project_urls={
        'Documentation': 'https://ysig.github.io/GraKeL/dev/',
        'Send us Feedback!': 'http://www.lix.polytechnique.fr/dascim/contact/',
        'Source': 'https://github.com/ysig/GraKeL/tree/develop',
        'Tracker': 'https://github.com/ysig/GraKeL/issues',
        },
      author='Ioannis Siglidis [LiX / DaSciM]',
      author_email='ioannis.siglidis@inria.fr',
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
