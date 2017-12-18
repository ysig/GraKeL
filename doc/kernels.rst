Kernels (between graphs)
========================

The documentation of the `kernels` submodule.

.. automodule:: grakel.kernels
   :members:

A list of default functions
---------------------------

:math:`l_{v}^{default}(x) = \begin{cases} 1 & x\neq \O \\ 0 & \text{else} \\ \end{cases}`
   
:math:`k_{v}^{default}(x,y,L_{x},L_{y}) = \begin{cases} 1 & L_{x}(x)= L_{y}(y) \\ 0 & \text{else} \\ \end{cases}`

:math:`k_{e}^{default}(x, y, E_{x}, E_{y}, L^{E}_{x}, L^{E}_{y})= \begin{cases} 1 & (x\in E_{x} \wedge y\in E_{y} \wedge L^{E}_{x}(x)= L^{E}_{y}) \vee (x\not\in E_{x} \wedge y\not\in E_{y}) \\ 0 & \text{else} \\ \end{cases}`

Bibliography
------------    
.. bibliography:: graph_kernels.bib
