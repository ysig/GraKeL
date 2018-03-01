.. _installation:

=================
Installing grakel
=================
The grakel library requires:

* Python [>=2.7 or >= 3.5]
* NumPy [>= 1.8.2]
* SciPy [>= 0.13.3]
* Cython [>= 0.27.3]
* cvxopt [>= 1.1.9]
* pynauty [>= 0.6.0]

For installing NumPy, SciPy, cvxopt and Cython the procedure
is the well known: `pip install extension==extension_version`.

----------------------
Installing **pynauty**
----------------------
The package `pynauty`_ is a very interesting python module by P. Dobcsányi, that wraps
in the c++ code of the `nauty`_ library of B. D. McKay, used for solving the decisional
graph isomorphism problem in a very good complexity.

Unfortunately the `pynauty`_ package, currently does not support the *pip* installation
procedure and should be installed manually as follows the installation procedure on the
package's web-page. To make things easier for the user, we have created a python script
under the name `install_pynauty.py` on the `official git repo`_ that automates the 
installation procedure. Allthough this script allows user configuration and is used
for the package build in `travis`_, `circle_ci`_, `appveyor`_ platforms, used
for automatic integration of code and documentation, it may not execute correctly
in all users and all environment. Currently the script supports Linux and Windows
installation (although for the second you should check the Windows Compatibility 
Issues) section, in the current page.

.. _pynauty: https://web.cs.dal.ca/~peter/software/pynauty/html/
.. _nauty: http://users.cecs.anu.edu.au/~bdm/nauty/
.. _official git repo: https://github.com/ysig/GraKeL
.. _travis: https://travis-ci.com/
.. _circle_ci: https://circleci.com/
.. _appveyor: https://www.appveyor.com/

---------------------------
Windows Comatibility Issues
---------------------------
There are two major compatibility issues concerning windows. 
The first one concerning `cvxopt`_ has to do with the fact that `cvxopt developers`_
use for their `python build a package called Mingwpy`_ which constrains the valid
Python 3 version for Windows to 3.4. This huge drawback restricts our project
compatibility to Python-3.4 except if user installs the `Python prebuilt libraries of Christoph Gohlke`_.

The other compatibility issue concerns pynauty. As mentioned above `pynauty`_ wraps
the `nauty`_ package, which uses a totally *gcc-oriented* approach and development.
As a result, the pynauty package does not build on the visual studio - c++ native of 
Windows, while it only could be bypassed if user uses only the *make* and *gcc/g++*
from the *minimalist gnu-like* project for Windows, called `MinGW`_. Of course this
project is covered through a bigger set of tools for linux-windows native command
compatibility known as `msys`_. We should note here, that for the automatic installation
of the pynauty script to succeed on windows, through the use of the pynauty script 
the user should have `curl`_ and `tar`_ gnu-like executables discoverable inside his
Windows Path (which both are covered in `msys`_). In comparison with the similar and 
more general project known as `cygwin`_, for `MinGW`_ there seems to be a native
`Python compiler compatibility`_. Although this seems to be it, more problems arise
as we dive through, mainly an issue concerning an installation on Python-3.4 which
was firstly tested, because of the `cvxopt`_ compatibility on Windows. Although
in such case the project builds normally, using as compiler the mingw, upon import
of the pynauty package by the simple command: `import pynauty` on a python shell
or inside a script, the system instantly loads a bulk of 4GB ram, which is `a known issue`_
and continues to remain a 1/4 fraction after the whole import, possibly restricting
execution speed.
As for the installation using pre-built packages, we haven't yet been able to see if 
this error reproduces on Python-3.5. So as for now there is no recommended installation
procedure that will produce an end result, as in a unix like environment.

Contributions in overriding or overcoming compatibility issues, would be widely accepted
and for doing so see the section :ref:`contributing`.

.. _cvxopt: http://cvxopt.org/
.. _cvxopt developers: http://cvxopt.org/copyright.html
.. _cygwin: https://www.cygwin.com/
.. _python build a package called Mingwpy: http://cvxopt.org/install/index.html#windows
.. _Python prebuilt libraries of Christoph Gohlke: https://www.lfd.uci.edu/~gohlke/pythonlibs/#cvxopt
.. _MinGW: http://www.mingw.org/
.. _msys: http://www.msys2.org/
.. _curl: https://curl.haxx.se/
.. _tar: https://www.gnu.org/software/tar/
.. _Python compiler compatibility: https://wiki.python.org/moin/WindowsCompilers
.. _a known issue: https://github.com/ContinuumIO/anaconda-issues/issues/271

---------------------
Why so many packages?
---------------------
The field of computationally efficient `Graph Kernels`_ can be considered
more as a collection of techniques calculating `PSD`_ similarity matrices, between
graph objects, than a field (with the compact sense of the term). The big diversity
of methods and ideas implied to extract a similarity score between to graphs, leads
to the use of very precise and studied algorithms as graph isomorphism (for determining isomorphic graphlets on the graphlet-sampling kernel) or convex optimization (through
the use of semidefinite-programming for calculating lovasz-theta embeddings, on lovasz-theta kernel) which on the other hand will not appear in any other kernel at all.
As a result import of external libraries that have studied and optimized the solutions
of such problems in detail, gives a complexity - implementation "standard" to refer to, 
while follows the idea of standing on the shoulders of the scientific community, which
can be viable through the use of `free software`_.

.. _Graph Kernels: https://en.wikipedia.org/wiki/Graph_kernel
.. _PSD: https://en.wikipedia.org/wiki/Positive-definite_matrix
.. _free software: https://en.wikipedia.org/wiki/Free_software