# GraKeL: A library for graph kernels
**!!!Warning: Project still under construction**

[![Travis Status](https://travis-ci.org/ysig/GraKeL.svg?branch=develop)](https://travis-ci.org/ysig/GraKeL)
[![Coverage Status](https://coveralls.io/repos/github/ysig/GraKeL/badge.svg?branch=develop)](https://coveralls.io/github/ysig/GraKeL?branch=develop)
[![CircleCI Status](https://circleci.com/gh/ysig/GraKeL/tree/develop.svg?style=shield)](https://circleci.com/gh/ysig/GraKeL/tree/develop)
**project-template** A compatible library based on a template of the project
[scikit-learn](http://scikit-learn.org/) 

This library implements all the basic and advanced graph-kernels.

## Installation and Usage
To install the package one needs to have installed the packages
`numpy` and `scipy` and `pynauty`. For the last one a bash script
is provide inside the project folder, in order to make the installation
procedure easier.

To install the module execute:
```shell
$ python setup.py install
```

If the installation is successful, and `grakel` is correctly installed,
you should be able to execute the following in Python:
```python
>>> import grakel

>>> X = np.array([[0,1,2,1,0,0.5,0], [1,0,0,0,1,2,0.5], [2,0,0,3,0,0,2], [1,0,3,0,0,0,0], [0,1,0,0,0,3,1], [0.5,2,0,0,3,0,0], [0,0.5,2,0,1,0,0]])
>>> L = {0:'banana', 1:'cherry', 2:'banana', 3:'cherry', 4:'peach',5:'cherry',6:'lime'}
>>> result = grakel.kernels.dirac(X, X, L, L)
>>> print("dirac:",result)
```

