<p align="center">
  <img width="50%" src="https://raw.githubusercontent.com/ysig/GraKeL/0.1a7/doc/_figures/logo.svg?sanitize=true" />
</p>

--------------------------------------------------------------------------------

[![Pypi Versions](https://img.shields.io/pypi/pyversions/grakel.svg)](https://pypi.org/pypi/grakel/)
[![Coverage Status](https://codecov.io/gh/ysig/GraKeL/branch/master/graph/badge.svg)](https://codecov.io/gh/ysig/GraKeL)
[![Travis Status](https://travis-ci.org/ysig/GraKeL.svg?branch=master)](https://travis-ci.org/ysig/GraKeL)
[![Appveyor status](https://ci.appveyor.com/api/projects/status/sss4lpfxwgejn6de?svg=true)](https://ci.appveyor.com/project/ysig/grakel)
[![CircleCI Status](https://circleci.com/gh/ysig/GraKeL.svg?style=svg)](https://circleci.com/gh/ysig/GraKeL)

**[Documentation](https://ysig.github.io/GraKeL/)** | **[Paper](http://jmlr.org/papers/volume21/18-370/18-370.pdf)**

*GraKeL* is a library that provides implementations of several well-established graph kernels. The library unifies these kernels into a common framework. Furthermore, it provides implementations of some frameworks that work on top of graph kernels. Specifically, GraKeL contains 15 kernels and 2 frameworks. The library is compatible with the [scikit-learn](http://scikit-learn.org/) pipeline allowing easy and fast integration inside machine learning algorithms.

--------------------------------------------------------------------------------

In detail, the following kernels and frameworks are currently implemented:

* **[Vertex histogram kernel](https://ysig.github.io/GraKeL/latest/generated/grakel.VertexHistogram.html)**
* **[Edge histogram kernel](https://ysig.github.io/GraKeL/latest/generated/grakel.EdgeHistogram.html)**
* **[Shortest path kernel](https://ysig.github.io/GraKeL/latest/generated/grakel.ShortestPath.html)** from Borgwardt and Kriegel: [Shortest-path kernels on graphs](https://www.dbs.ifi.lmu.de/~borgward/papers/BorKri05.pdf) (ICDM 2005)
* **[Graphlet kernel](https://ysig.github.io/GraKeL/latest/generated/grakel.GraphletSampling.html)** from Shervashidze *et al.*: [Efficient graphlet kernels for large graph comparison](http://proceedings.mlr.press/v5/shervashidze09a/shervashidze09a.pdf) (AISTATS 2009)
* **[Random walk kernel](https://ysig.github.io/GraKeL/latest/generated/grakel.RandomWalk.html)** from Vishwanathan *et al.*: [Graph Kernels](http://www.jmlr.org/papers/volume11/vishwanathan10a/vishwanathan10a.pdf) (JMLR 11(Apr))
* **[Neighborhood hash graph kernel](https://ysig.github.io/GraKeL/latest/generated/grakel.NeighborhoodHash.html)** from Hido and Kashima: [A Linear-time Graph Kernel](https://ieeexplore.ieee.org/abstract/document/5360243) (ICDM 2009)
* **[Weisfeiler-Lehman framework](https://ysig.github.io/GraKeL/latest/generated/grakel.WeisfeilerLehman.html)** from Shervashidze *et al.*: [Weisfeiler-Lehman Graph Kernels](http://www.jmlr.org/papers/volume12/shervashidze11a/shervashidze11a.pdf) (JMLR 12(Sep))
* **[Neighborhood subgraph pairwise distance kernel](https://ysig.github.io/GraKeL/latest/generated/grakel.NeighborhoodSubgraphPairwiseDistance.html)** from Costa and De Grave: [Fast Neighborhood Subgraph Pairwise Distance Kernel](https://pdfs.semanticscholar.org/7a10/f6a406b664d1159e7c4fefbdd6ac275aee53.pdf) (ICML 2010)
* **[Lovasz-theta kernel](https://ysig.github.io/GraKeL/latest/generated/grakel.LovaszTheta.html)** from Johansson *et al.*: [Global graph kernels using geometric embeddings](http://proceedings.mlr.press/v32/johansson14.pdf) (ICML 2014)
* **[SVM-theta kernel](https://ysig.github.io/GraKeL/latest/generated/grakel.SvmTheta.html)** from Johansson *et al.*: [Global graph kernels using geometric embeddings](http://proceedings.mlr.press/v32/johansson14.pdf) (ICML 2014)
* **[Ordered decompositional DAG kernel](https://ysig.github.io/GraKeL/latest/generated/grakel.OddSth.html)** from Da San Martino *et al.*: [A Tree-Based Kernel for Graphs](https://pdfs.semanticscholar.org/69ee/18dd7a214d4d656b5b95742212f050dabeac.pdf) (SDM 2012)
* **[GraphHopper kernel](https://ysig.github.io/GraKeL/latest/generated/grakel.GraphHopper.html)** from Feragen *et al.*: [Scalable kernels for graphs with continuous attributes](https://papers.nips.cc/paper/5155-scalable-kernels-for-graphs-with-continuous-attributes.pdf) (NIPS 2013)
* **[Propagation kernel](https://ysig.github.io/GraKeL/latest/generated/grakel.Propagation.html)** from Neumann *et al.*: [Propagation kernels: efficient graph kernels from propagated information](https://link.springer.com/content/pdf/10.1007/s10994-015-5517-9.pdf) (Machine Learning 102(2))
* **[Pyramid match kernel](https://ysig.github.io/GraKeL/latest/generated/grakel.PyramidMatch.html)** from Nikolentzos *et al.*: [Matching Node Embeddings for Graph Similarity](https://www.aaai.org/ocs/index.php/AAAI/AAAI17/paper/view/14494/14426) (AAAI 2017)
* **[Subgraph matching kernel](https://ysig.github.io/GraKeL/latest/generated/grakel.SubgraphMatching.html)** from Kriege and Mutzel: [Subgraph Matching Kernels for Attributed Graphs](https://arxiv.org/ftp/arxiv/papers/1206/1206.6483.pdf) (ICML 2012)
* **[Multiscale Laplacian kernel](https://ysig.github.io/GraKeL/latest/generated/grakel.MultiscaleLaplacian.html)** from Kondor and Pan: [The Multiscale Laplacian Graph Kernel](https://papers.nips.cc/paper/6135-the-multiscale-laplacian-graph-kernel.pdf) (NIPS 2016)
* **[Core framework](https://ysig.github.io/GraKeL/latest/generated/grakel.CoreFramework.html)** from Nikolentzos *et al.*: [A Degeneracy Framework for Graph Similarity](https://www.ijcai.org/proceedings/2018/0360.pdf) (IJCAI 2018)
* **[Weisfeiler-Lehman optimal assignment kernel](https://ysig.github.io/GraKeL/latest/generated/grakel.WeisfeilerLehmanOptimalAssignment.html)** from Kriege *et al.*: [On Valid Optimal Assignment Kernels and Applications to Graph Classification](http://papers.nips.cc/paper/6166-on-valid-optimal-assignment-kernels-and-applications-to-graph-classification.pdf) (NIPS 2016)

--------------------------------------------------------------------------------

To learn how to install and use GraKeL, and to find out more about the implemented kernels and frameworks, please read our [documentation](https://ysig.github.io/GraKeL/). To learn about the functionality of the library and about example applications, check out our [examples](https://github.com/ysig/GraKeL/tree/master/examples) in the `examples/` directory and our [tutorials](https://github.com/ysig/GraKeL/tree/master/tutorials) in the `tutorials/` directory.

In case you find a bug, please open an [issue](https://github.com/ysig/GraKeL/issues). To propose a new kernel, you can open a [feature request](https://github.com/ysig/GraKeL/issues).

## Installation

The GraKeL library requires the following packages to be installed:

* Python (>=2.7, >=3.5)
* NumPy (>=1.8.2)
* SciPy (>=0.13.3)
* Cython (>=0.27.3)
* cvxopt (>=1.2.0) [optional]
* future (>=0.16.0) (for python 2.7)

To install the package, run:

```sh
$ pip install grakel
```

## Running tests

To test the package, execute:
```sh
$ nosetests grakel
```

## Running examples

```
$ cd examples
$ python shortest_path.py
```

## Cite

If you use GraKeL in a scientific publication, please cite our paper (http://jmlr.org/papers/volume21/18-370/18-370.pdf):

```
@article{JMLR:v21:18-370,
  author  = {Giannis Siglidis and Giannis Nikolentzos and Stratis Limnios and Christos Giatsidis and Konstantinos Skianis and Michalis Vazirgiannis},
  title   = {GraKeL: A Graph Kernel Library in Python},
  journal = {Journal of Machine Learning Research},
  year    = {2020},
  volume  = {21},
  number  = {54},
  pages   = {1-5}
}
```

## License

GraKeL is distributed under the __BSD 3-clause__ license. The library makes use of the C++ source code of [BLISS](http://www.tcs.hut.fi/Software/bliss) (a tool for computing automorphism groups and canonical labelings of graphs) which is __LGPL__ licensed. Futhermore, the [cvxopt](https://cvxopt.org/) package (a software package for convex optimization) which is an optional dependency of GraKeL is __GPL__ licensed.
