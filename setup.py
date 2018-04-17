"""The general `setup.py`  file."""
from __future__ import print_function
import sys
from warnings import warn
from platform import system
from setuptools import setup, find_packages, Extension

with open('requirements.txt') as f:
    INSTALL_REQUIRES = [l.strip() for l in f.readlines() if l]

OS = system()

if OS == 'Windows':
    warn('Installation in Windows is incomplete..'
         'To see why please check the documentation.')
    extra_compile_args = ["/O2"]
elif OS in ['Linux', 'Darwin']:
    extra_compile_args = ["-O3"]

try:
    from Cython.Distutils import build_ext
except ImportError:
    print('build_ext from Cython.Distutils is required during installation',
          file=sys.stderr)
    sys.exit(1)

"""
try:
    import pynauty
except ImportError as ie:
    if OS in ['Linux', 'Darwin']:
        from install_pynauty import main_unix
        from subprocess import CalledProcessError
        try:
            main_unix()
        except CalledProcessError:
            print('The automatic script-failed to install pynauty..\n'
                  'Try to install it on your own.', file=sys.stderr)
            sys.exit(1)
    else:
        print('problematic OS -- pynauty should be installed manually', file=sys.stderr)
        sys.exit(1)
"""

try:
    import numpy
except ImportError:
    print('numpy is required during installation', file=sys.stderr)
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
                extra_compile_args=extra_compile_args)

setup(name='grakel',
      version='0.1a2',
      description='A scikit-learn compatible library for graph kernels',
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
      python_requires='>=2.7,>=3.5',
      packages=find_packages(),
#      entry_points = {
#        'console_scripts': ['install-pynauty=install_pynauty:main']
#      },
      install_requires=INSTALL_REQUIRES,
      extras_require={
        'graphlet':  ["pynauty>=0.6.0"],
        'lovasz': ["cvxopt>=1.1.9"]
      },
      ext_modules=[ext],
      cmdclass={'build_ext': build_ext}
      )
