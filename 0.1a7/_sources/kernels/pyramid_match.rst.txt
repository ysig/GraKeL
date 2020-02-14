.. _pyramid_match:

Pyramid Match Kernel
====================

The pyramid match kernel is a very popular algorithm in Computer Vision, and has proven useful for many applications including object recognition and image retrieval :cite:`grauman2007pyramid`, :cite:`lazebnik2006beyond`.
The pyramid match graph kernel extends its applicability to graph-structured data :cite:`nikolentzos2017matching`.
The kernel can handle both unlabeled graphs and graphs that contains discrete node labels.

The pyramid match graph kernel first embedds the vertices of each graph into a low-dimensional vector space using the eigenvectors of the :math:`d` largest in magnitude eigenvalues of the adjacency matrix of the graph.
Since the signs of these eigenvectors are arbitrary, it replaces all their components by their absolute values.
Each vertex is thus a point in the :math:`d`-dimensional unit hypercube.
To find an approximate correspondence between the sets of vertices of two graphs, the kernel maps these points to multi-resolution histograms, and compares the emerging histograms with a weighted histogram intersection function.

Initially, the kernel partitions the feature space into regions of increasingly larger size and takes a weighted sum of the matches that occur at each level.
Two points match with each other if they fall into the same region.
Matches made within larger regions are weighted less than those found in smaller regions.
The kernel repeatedly fits a grid with cells of increasing size to the :math:`d`-dimensional unit hypercube.
Each cell is related only to a specific dimension and its size along that dimension is doubled at each iteration, while its size along the other dimensions stays constant and equal to :math:`1`.
Given a sequence of levels from :math:`0` to :math:`L`, then at level :math:`l`, the :math:`d`-dimensional unit hypercube has :math:`2^l` cells along each dimension and :math:`D = 2^{l}d` cells in total.
Given a pair of graphs :math:`G,G'`, let :math:`H_G^l` and :math:`H_{G'}^l` denote the histograms of :math:`G` and :math:`G'` at level :math:`l` and :math:`H_G^l(i)`, :math:`H_{G'}^l(i)`, the number of vertices of :math:`G`, :math:`G'` that lie in the :math:`i^{th}` cell.
The number of points in two sets which match at level $l$ is then computed using the histogram intersection function

.. math::

  I(H_G^l,H_{G'}^l) = \sum_{i=1}^D \min\big(H_G^l(i),H_{G'}^l(i)\big)

The matches that occur at level :math:`l` also occur at levels :math:`0, \ldots, l-1`.
We are interested in the number of new matches found at each level which is given by :math:`I(H_{G_1}^l,H_{G_2}^l) - I(H_{G_1}^{l+1},H_{G_2}^{l+1})` for :math:`l=0,\ldots,L-1`.
The number of new matches found at each level in the pyramid is weighted according to the size of that level's cells.
Matches found within smaller cells are weighted more than those made in larger cells.
Specifically, the weight for level :math:`l` is set equal to :math:`\frac{1}{2^{L-l}}`.
Hence, the weights are inversely proportional to the length of the side of the cells that varies in size as the levels increase.
The pyramid match kernel is then defined as follows

.. math::

  k(G,G') = I(H_G^L,H_{G'}^L) + \sum_{l=0}^{L-1} \frac{1}{2^{L-l}}\big(I(H_G^l,H_{G'}^l) - I(H_G^{l+1},H_{G'}^{l+1})\big)

The complexity of the pyramid match kernel is :math:`\mathcal{O}(dnL)` where $n$ is the number of nodes of the graphs under comparison.

In the case of labeled graphs, the kernel restricts matchings to occur only between vertices that share same labels.
It represents each graph as a set of sets of vectors, and matches pairs of sets of two graphs corresponding to the same label using the pyramid match kernel.
The emerging kernel for labeled graphs corresponds to the sum of the separate kernels

.. math::

    k(G, G') = \sum_{i=1}^c k^i(G,G')

where :math:`c` is the number of distinct labels and :math:`k^i(G_1,G_2)` is the pyramid match kernel between the sets of vertices of the two graphs which are assigned the label :math:`i`.

The above kernel is implemented below

.. currentmodule:: grakel

.. autosummary::

   PyramidMatch

Bibliography
------------
.. bibliography:: graph_kernels.bib
   :filter: docname in docnames
