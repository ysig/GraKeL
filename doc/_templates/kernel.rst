:mod:`{{module}}`.{{objname}}
{{ underline }}==============

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}

   {% block methods %}
   .. automethod:: __init__
   {% endblock %}


Bibliography
------------
.. bibliography:: ../kernels/graph_kernels.bib
   :filter: docname in docnames

.. include:: {{module}}.{{objname}}.examples

.. raw:: html

    <div class="clearer"></div>
