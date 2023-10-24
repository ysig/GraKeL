"""The general `setup.py` file."""
# Author: Ioannis Siglidis <y.siglidis@gmail.com>
# License: BSD 3 clause
import sys
from platform import system

from setuptools import Extension, find_packages, setup
from numpy import get_include
from Cython.Build import build_ext

# Compile extensions

# Set optimization arguments for compilation
OS = system()
if OS == "Windows":
    extra_compile_args = ["/O2", "/w"]
elif OS in ["Linux", "Darwin"]:
    extra_compile_args = ["-O3", "-w"]

# Remove the "-Wstrict-prototypes" compiler option, which isn't valid for C++.
import distutils.sysconfig

cfg_vars = distutils.sysconfig.get_config_vars()
for key, value in cfg_vars.items():
    if type(value) == str:
        cfg_vars[key] = value.replace("-Wstrict-prototypes", "")

# Add the _c_functions extension on kernels
ext_address = "./grakel/kernels/_c_functions/"
ext = Extension(
    name="grakel.kernels._c_functions",
    sources=[
        ext_address + "functions.pyx",
        ext_address + "src/ArashPartov.cpp",
        ext_address + "src/sm_core.cpp",
    ],
    include_dirs=[ext_address + "include", get_include()],
    depends=[ext_address + "include/functions.hpp"],
    language="c++",
    extra_compile_args=extra_compile_args,
)

# Add the bliss library extension for calculating isomorphism
isodir = "./grakel/kernels/_isomorphism/"
blissdir = isodir + "bliss-0.50/"

# The essential bliss source files
blisssrcs = ["graph.cc", "heap.cc", "orbit.cc", "partition.cc", "uintseqhash.cc"]
blisssrcs = [blissdir + src for src in blisssrcs]
pn = str(sys.version_info[0])

# Compile intpybliss
intpybliss = Extension(
    name="grakel.kernels._isomorphism.intpybliss",
    define_macros=[("MAJOR_VERSION", "0"), ("MINOR_VERSION", "50beta")],
    include_dirs=[blissdir],
    language="c++",
    extra_compile_args=extra_compile_args,
    sources=[isodir + "intpyblissmodule_" + pn + ".cc"] + blisssrcs,
)

# Make bliss extension
bliss = Extension(
    name="grakel.kernels._isomorphism.bliss",
    include_dirs=[isodir],
    language="c++",
    extra_compile_args=extra_compile_args,
    sources=[isodir + "bliss.pyx"],
)

setup(
    packages=find_packages(),
    package_data={"grakel.tests": ["data/Cuneiform/*.txt", "data/MUTAG/*.txt"]},
    ext_modules=[intpybliss, bliss, ext],
    cmdclass={"build_ext": build_ext},
)
