Overview
========

Each structural descriptors inherits from the base class :py:class:`StructuralDescriptor <partycls.descriptors.descriptor.StructuralDescriptor>` and is computed for a given :py:class:`Trajectory <partycls.trajectory.Trajectory>`.

Most descriptors require information about the neighbors of each of the considered particles. Nearest neighbors, *i.e.* particles inside of the first coordination shell, are computed at the level of the trajectory itself or directly read from the input trajectory file. This information is then propagated to the descriptors that rely on it to perform the computations.

There are several method to identify the neighbors of the central particles, and there are different ways for this information to be processed by the descriptors.

Nearest neighbors
-----------------

:py:class:`Trajectory <partycls.trajectory.Trajectory>` instances have several attributes and methods to identify nearest neighbors and determine the cutoffs delimiting the first coordination shell.

Nearest neighbors can be computed using the following method:

.. automethod:: partycls.trajectory.Trajectory.compute_nearest_neighbors

Nearest neighbors cutoffs can be automatically computed for each pair of species in the trajectory on the basis of the first minimum of the partial radial distribution functions :math:`g_{\alpha\beta}(r)`, where :math:`\alpha` and :math:`\beta` are species indices. This is done using the following method:

.. automethod:: partycls.trajectory.Trajectory.compute_nearest_neighbors_cutoffs

Alternatively, both the nearest neighbors method and cutoffs can be set as instance attributes of :py:class:`Trajectory <partycls.trajectory.Trajectory>`:

.. code-block:: python

	traj = Trajectory('trajectory.xyz')
	traj.nearest_neighbors_method = 'fixed'
	traj.nearest_neighbors_cutoffs = [1.5, 1.4, 1.4, 1.3]
	traj.compute_nearest_neighbors()

Once computed, the lists of nearest neighbors are stored as :py:class:`Particle <partycls.particle.Particle>` attributes. These can be accessed through the :py:meth:`Trajectory.get_property() <partycls.trajectory.Trajectory.get_property>` method, or by iterating over the particles:

.. code-block:: python

	neighbors = traj.get_property('neighbors')
	print(neighbors)

.. code-block:: python

	for system in traj:
		for particle in system.particle:
			print(particle.nearest_neighbors)

We now present in more details the various available methods to identify the nearest neighbors.

1. From the trajectory
~~~~~~~~~~~~~~~~~~~~~~

Nearest neighbors can be directly from the input trajectory file if it contains an additional column such as ``neighbors`` or ``nearest_neighbors``. Below is an example of such a file in XYZ format: 

.. code-block::

	100
	columns:id,pos,neighbors cell:5.000,5.000,5.000
	A -1.100 -2.166 -0.629 9,12,54,74
	A -1.754  0.583  1.231 2,27,63
	A -0.338  1.957 -1.365 4,45,56,78,81
	B  1.030 -0.220 -1.256 14,31,35
	B  1.322 -1.556  2.134 41,63,70,92
	...

This, however, must be specified when reading the input file through the ``additional_fields`` parameter.

.. warning::
	Currently, this only works for trajectory files in XYZ format.

Example:

.. code-block:: python

	traj = Trajectory('trajectory.xyz', additional_fields=['neighbors'])

2. Fixed-cutoffs
~~~~~~~~~~~~~~~~

Set using one of:

- ``Trajectory.nearest_neighbors_method = 'fixed'``
- ``Trajectory.compute_nearest_neighbors(method='fixed')``

Nearest neighbors are defined on the basis of a fixed cutoff distance :math:`r_{\alpha\beta}^c`, where :math:`\alpha` and :math:`\beta` are species indices. The cutoff distance is equal to the first minimum of the corresponding partial radial distribution function,  :math:`g_{\alpha\beta}(r)`.

Example:

.. code-block:: python

	traj = Trajectory('trajectory.xyz')
	traj.compute_nearest_neighbors(method='fixed', cutoffs=[1.5, 1.4, 1.4, 1.3])

3. Solid-angle based nearest neighbors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Set using one of:

- ``Trajectory.nearest_neighbors_method = 'sann'``
- ``Trajectory.compute_nearest_neighbors(method='sann')``

The *solid-angle based nearest neighbors* algorithm (SANN) :cite:`van_meel_2012` is a parameter-free algorithm for the identification of nearest neighbors. It attributes to each possible neighbor of a particle a solid angle and determines the cutoff radius by the requirement that the sum of the solid angles is :math:`4 \pi`.

.. important ::
	This method requires cutoffs (or computes them automatically if not provided) to use as a first guess to identify the possible nearest neighbors. However, cutoffs do not play a role in the algorithm itself. A good choice for these cutoffs is the first minima of the partial radial distribution functions :math:`g_{\alpha\beta}(r)`.

Example:

.. code-block:: python

	traj = Trajectory('trajectory.xyz')
	traj.compute_nearest_neighbors(method='sann', cutoffs=[1.5, 1.4, 1.4, 1.3])

4. Radical Voronoi neighbors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Set using one of:

- ``Trajectory.nearest_neighbors_method = 'voronoi'``
- ``Trajectory.compute_nearest_neighbors(method='voronoi')``

Voronoi tessellation can be used in molecular simulations to identify nearest neighbors by construction of Voronoi polyhedra :cite:`bernal_1959`, which consists in drawing orthogonal planes at the mid-points between the central particle and each of its neighbors (*i.e.* its `Wigner-Seitz cell <https://en.wikipedia.org/wiki/Wigner%E2%80%93Seitz_cell>`_). In particular, the *radical* variant of Voronoi tessellation :cite:`gellatly_1982`, which accounts for the relative sizes of the particles to determine of the positions of the intersecting planes, provides better results for multi-components systems.

.. warning::
	This method uses the effective particles' radii. This information must thus be provided either from an additional field in the input trajectory file, or directly at the level of the :py:class:`Trajectory <partycls.trajectory.Trajectory>` instance using the :py:meth:`Trajectory.set_property() <partycls.trajectory.Trajectory.set_property>` method. Otherwise, default values will be used.

Examples:

1. Particles' radii are read from the input trajectory file:

.. code-block:: python

	traj = Trajectory('trajectory.xyz', additional_fields=['radius'])
	traj.compute_nearest_neighbors(method='voronoi')

2. Particles' radii are set in the :py:class:`Trajectory <partycls.trajectory.Trajectory>`:

.. code-block:: python

	traj = Trajectory('trajectory.xyz')
	traj.set_property('radius', 0.5, subset="species == 'A'")
	traj.set_property('radius', 0.4, subset="species == 'B'")
	traj.compute_nearest_neighbors(method='voronoi')

Neighbors & structural descriptors
----------------------------------

The table below shows the requirements of each structural descriptors in terms of neighbors and cutoffs. Note that these requirements are automatically satisfied when computing the descriptors if not explicitly set by the user.

.. list-table:: Required types of neighbors & cutoffs
	:widths: 20 20 20 20
	:header-rows: 1

	* - 
	  - Nearest [1]_
	  - Extended [2]_
	  - Cutoffs [3]_
	* - :doc:`gr`
	  - ✕
	  - ✕
	  - ✕
	* - :doc:`tetra`
	  - √
	  - ✕
	  - ✕
	* - :doc:`ba`
	  - √
	  - ✕
	  - ✕
	* - :doc:`bo`
	  - √
	  - ✕
	  - ✕
	* - :doc:`labo`
	  - √
	  - ✕
	  - ✕
	* - :doc:`compact`
	  - √
	  - ✕
	  - ✕
	* - :doc:`coord`
	  - √
	  - ✕
	  - ✕
	* - :doc:`rbo`
	  - ✕
	  - √
	  - ✕
	* - :doc:`sba`
	  - ✕
	  - √
	  - √
	* - :doc:`sbo`
	  - ✕
	  - √
	  - √

.. [1] Nearest neighbors, computed in the :py:class:`Trajectory <partycls.trajectory.Trajectory>` on the basis of one the nearest neighbors methods presented above.
.. [2] Extended neighbors at **fixed distances**, *i.e.* beyond the first coordination shell. These are computed directly in the descriptor.
.. [3] Nearest neighbors cutoffs :math:`\{ r_{\alpha\beta}^c \}` for normalization. These are computed automatically if not provided in the :py:class:`Trajectory <partycls.trajectory.Trajectory>`.

References
----------

.. bibliography:: ../../references.bib
	:style: unsrt
	:filter: docname in docnames
