Tutorials
=========

.. toctree::
	:hidden:
	
	 Home <self>

Detailed examples
-----------------

A collection of `Jupyter notebooks <https://jupyter.org/>`_, with examples and detailed instructions on how to run the code:

1. :doc:`Manipulating a trajectory <tutorials/notebooks/1_trajectory>` (:download:`download <tutorials/notebooks/1_trajectory.ipynb>`)
2. :doc:`Setting and running a workflow <tutorials/notebooks/2_workflow>` (:download:`download <tutorials/notebooks/2_workflow.ipynb>`)
3. :doc:`Setting and computing descriptors <tutorials/notebooks/3_descriptors>` (:download:`download <tutorials/notebooks/3_descriptors.ipynb>`)
4. :doc:`More advanced applications <tutorials/notebooks/4_going_further>` (:download:`download <tutorials/notebooks/4_going_further.ipynb>`)

.. important::
	If you download the notebooks, be sure to download the content of the `data <https://github.com/jorisparet/partycls/tree/master/data>`_ folder from GitHub as well before executing them. The data folder should be placed one level above the folder containing the notebooks (*i.e.* :file:`../data/`) for the paths to remain correct.

.. toctree::
	:glob:
	:hidden:
	:maxdepth: 2
	:numbered:

	tutorials/notebooks/*

Structural descriptors
----------------------

An overview and a detailed presentation of the built-in structural descriptors:	

- :doc:`tutorials/descriptors/overview`
- :doc:`tutorials/descriptors/gr`
- :doc:`tutorials/descriptors/tetra`
- :doc:`tutorials/descriptors/ba`
- :doc:`tutorials/descriptors/sba`
- :doc:`tutorials/descriptors/bo`
- :doc:`tutorials/descriptors/sbo`
- :doc:`tutorials/descriptors/ld`
- :doc:`tutorials/descriptors/rbo`
- :doc:`tutorials/descriptors/compact`

You can also check out the API of the :py:mod:`partycls.descriptor` subpackage.

.. toctree::
	:glob:
	:hidden:
	:maxdepth: 1

	tutorials/descriptors/overview
	tutorials/descriptors/gr
	tutorials/descriptors/tetra
	tutorials/descriptors/ba
	tutorials/descriptors/sba
	tutorials/descriptors/bo
	tutorials/descriptors/sbo
	tutorials/descriptors/ld
	tutorials/descriptors/rbo
	tutorials/descriptors/compact
