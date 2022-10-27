partycls
========

.. image:: https://img.shields.io/pypi/l/partycls.svg
	:target: https://en.wikipedia.org/wiki/GNU_General_Public_License
	:alt: license badge
.. image:: https://img.shields.io/badge/python-3.6+-blue.svg
	#:target:
	:alt:  python version badge
.. image:: https://img.shields.io/pypi/v/partycls.svg
	:target: https://pypi.python.org/pypi/partycls/
	:alt: pyPI version badge
.. image:: https://img.shields.io/badge/coverage-89%25-green
	#:target:
	:alt: coverage badge
.. image:: https://joss.theoj.org/papers/10.21105/joss.03723/status.svg
	:target: https://doi.org/10.21105/joss.03723
	:alt: joss badge
.. image:: https://img.shields.io/badge/GitHub-100?logo=github&logoColor=white
	:target: https://github.com/jorisparet/partycls/
	:alt: github badge

**partycls** is a Python package for cluster analysis of systems of interacting particles. By grouping particles that share similar structural or dynamical properties, partycls enables rapid and unsupervised exploration of the system's relevant features. It provides descriptors suitable for applications in condensed matter physics, such as structural analysis of disordered or partially ordered materials, and integrates the necessary tools of unsupervised learning into a streamlined workflow.

You can find a more detailed presentation in our open-access article: `“partycls: A Python package for structural clustering” <https://joss.theoj.org/papers/10.21105/joss.03723>`_ .

Main features
-------------

Trajectory formats
~~~~~~~~~~~~~~~~~~

partycls accepts several trajectory formats (including custom ones) either through its built-in trajectory reader or via third-party packages, such as `MDTraj <www.mdtraj.org>`_ and `atooms <https://framagit.org/atooms/atooms>`_. The code is currently optimized for small and medium system sizes (of order :math:`10^4` particles). Multiple trajectory frames can be analyzed to extend the structural dataset.

Structural descriptors
~~~~~~~~~~~~~~~~~~~~~~

partycls implements various structural descriptors:

- :doc:`tutorials/descriptors/gr`
- :doc:`tutorials/descriptors/tetra`
- :doc:`tutorials/descriptors/ba`
- :doc:`tutorials/descriptors/sba`
- :doc:`tutorials/descriptors/bo`
- :doc:`tutorials/descriptors/sbo`
- :doc:`tutorials/descriptors/labo`
- :doc:`tutorials/descriptors/rbo`
- :doc:`tutorials/descriptors/compact`
- :doc:`tutorials/descriptors/coord`

On top of these native descriptors, partycls supports additional structural descriptors via `DScribe <https://singroup.github.io/dscribe>`_.

Unsupervised machine-learning
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

partycls performs feature scaling, dimensionality reduction and cluster analysis using the `scikit-learn <https://scikit-learn.org>`_ package and additional built-in algorithms.

Documentation
-------------

Check the :doc:`tutorials <tutorials>` to see various examples and detailed instructions on how to run the code, as well as an in-depth presentation of the built-in structural descriptors.

For a more detailed documentation, you can check the :doc:`API <api>`.

.. toctree::
	:hidden:

	Home <self>
	install
	tutorials
	api
	support
	citing
	changelog
	about
