#GraKeL: A library for graph kernels
##!!!Warning: Project still under construction

[![Travis Status](https://travis-ci.org/scikit-learn-contrib/project-template.svg?branch=master)](https://travis-ci.org/scikit-learn-contrib/project-template)
[![Coveralls Status](https://coveralls.io/repos/scikit-learn-contrib/project-template/badge.svg?branch=master&service=github)](https://coveralls.io/r/scikit-learn-contrib/project-template)
[![CircleCI Status](https://circleci.com/gh/scikit-learn-contrib/project-template.svg?style=shield&circle-token=:circle-token)](https://circleci.com/gh/scikit-learn-contrib/project-template/tree/master)

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

