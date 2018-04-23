# GraKeL: A library for graph kernels

[![Travis Status](https://travis-ci.org/ysig/GraKeL.svg?branch=develop)](https://travis-ci.org/ysig/GraKeL)
[![Coverage Status](https://coveralls.io/repos/github/ysig/GraKeL/badge.svg?branch=develop)](https://coveralls.io/github/ysig/GraKeL?branch=develop)
[![CircleCI Status](https://circleci.com/gh/ysig/GraKeL/tree/develop.svg?style=shield)](https://circleci.com/gh/ysig/GraKeL/tree/develop)
[![Appveyor status](https://ci.appveyor.com/api/projects/status/sss4lpfxwgejn6de/branch/develop?svg=true)](https://ci.appveyor.com/project/ysig/grakel)

**grakel** is a library compatible with the project of [scikit-learn](http://scikit-learn.org/)

Installing grakel
=================

The grakel library requires:

* Python [>=2.7, >= 3.5]
* NumPy [>= 1.8.2]
* SciPy [>= 0.13.3]
* Cython [>= 0.27.3]
* cvxopt [>= 1.1.9] [optional: lovasz]

For installing NumPy, SciPy, cvxopt and Cython the procedure
is the well known:
```shell
(sudo) pip install extension>=extension_version
```
depending on if `pip` has superuser privilages.

**A Windows user concern about dependencies.**
If the windows user wants to use the *lovasz-theta* kernel, which is the only one that requires cvxopt,
she/he shoud install it manually. This is restrictive for Python3, because Cython developers have not
supported yet versions bigger than 3.4. As a result installation of this optional dependency should be
met either by installing prebuilt libraries of anaconda (or of Christofer Gholke).

To install the module execute:
```shell
$ (sudo) pip install grakel-dev
```

Contributing
============
For learning how to read to integrate your own kernel, please read section *Write your own kernel* inside
the package [documentation][doc]. 
For contributing to the GraKeL project, please read section *contributing* inside the package [documentation][doc].

Usage
=====
To learn how to use the GraKeL api **as a user**, please read the [documentation][doc] on sections *Introduction* and *A longer introduction*
(in case your are full of curiosity).

Testing
=======
To test the package, execute:
```shell
$ nosetests
```

for executing unit_tests or use a testing-interface for testing the `kernel` module:
```shell
$ python  grakel/tests/test_kernels.py --help
usage: test_kernels.py [-h] [--verbose] [--problematic] [--slow]
                       [--ignore_warnings] [--dataset DATASET] [--normalize]
                       [--develop | --all | --main]

A test file for all kernels

optional arguments:
  -h, --help         show this help message and exit
  --verbose          print kernels with their outputs on stdout
  --problematic      allow execution of problematic test cases in development
  --slow             allow execution of slow test cases in development
  --ignore_warnings  ignore warnings produced by kernel executions
  --dataset DATASET  chose the datset you want the tests to be executed
  --normalize        normalize the kernel output
  --develop          execute only tests connected with current development
  --all              execute all tests
  --main             execute the main tests [default]

```

for testing `graph_kernels`:
```shell
$ python grakel/tests/test_graph_kernel.py --help
usage: test_graph_kernels.py [-h] [--verbose] [--problematic] [--slow]
                             [--normalize] [--ignore_warnings]
                             [--dataset DATASET] [--develop | --all | --main]

A test file for all kernels

optional arguments:
  -h, --help         show this help message and exit
  --verbose          print kernels with their outputs on stdout
  --problematic      allow execution of problematic test cases in development
  --slow             allow execution of slow test cases in development
  --normalize        normalize the kernel output
  --ignore_warnings  ignore warnings produced by kernel executions
  --dataset DATASET  chose the datset you want the tests to be executed
  --develop          execute only tests connected with current development
  --all              execute all tests
  --main             execute the main tests [default]

```

and for testing the `Graph` class:
```shell
$ python grakel/tests/test_graph.py --help
usage: test_graph.py [-h] [--verbose] [--ignore_warnings]

A test file for all `Graph` type objects

optional arguments:
  -h, --help         show this help message and exit
  --verbose          verbose outputs on stdout
  --ignore_warnings  ignore warnings produced by kernel executions
```
You can also execute the kernel test locally through a *test-main-function* as
```shell
$ python -m grakel.tests
```
but in order for this to work you would need first to build the package cython extension
locally by executing:
```shell
$ python setup.py build_ext -i
```

[doc]: https://ysig.github.io/GraKeL/dev/
