Radial descriptor
=================

Definition
----------

Let :math:`\mathbf{r}_i` be the position of particle :math:`i` and define :math:`r_{ij} = |\mathbf{r}_j - \mathbf{r}_i|` as the distance between particle :math:`i` and its neighbors :math:`j`.

We define :math:`n_i(r_m)` as the number of neighbors :math:`j` of particle :math:`i` for which :math:`r_{ij}` is between :math:`r_m = r_\mathrm{min} + m \times \Delta r` and :math:`r_{m+1} = r_\mathrm{min} + (m+1) \times \Delta r` :cite:`paret_2020`. Here, :math:`\Delta r` has the interpration of a bin width in a histogram.

We then consider :math:`n_i(r_m)` for a set of distances :math:`\{ d_n \}` separated by :math:`\Delta r`, :math:`\{ d_n \} = \{ r_0, r_1, \dots, r_{n_\mathrm{max}} \}`. The resulting feature vector for particle :math:`i` is given by

.. math::
	X^\mathrm{R}(i) = (\: n_i(r_0) \;\; n_i(r_1) \;\; \dots \;\; n_i(r_{n_\mathrm{max}}) \:) .

Setup
-----

Instantiating this descriptor on a ``Trajectory`` can be done as follows:

.. code-block:: python

	from partycls import Trajectory
	from partycls.descriptor import RadialDescriptor

	traj = Trajectory("trajectory.xyz")
	D = RadialDescriptor(traj)

The constructor takes the following parameters:

.. automethod:: partycls.descriptor.gr.RadialDescriptor.__init__

Demonstration
-------------

We consider an input trajectory file :file:`trajectory.xyz` in XYZ format that contains two particle types ``"A"`` and ``"B"``:

.. code-block:: python

	from partycls import Trajectory

	# open the trajectory
	traj = Trajectory("trajectory.xyz")

We now instantiate a ``RadialDescriptor`` on this trajectory and restrict the analysis to type-B particles only. We set :math:`\Delta r = 0.1` and :math:`(r_\mathrm{min},r_\mathrm{max}) = (1, 1.5)` to set the grid of distances :math:`\{d_n\}`:

.. code-block:: python

	from partycls.descriptor import RadialDescriptor

	# instantiation
	D = RadialDescriptor(traj, dr=0.1, bounds=(1.0, 1.5))

	# print the grid of angles
	print("grid:\n", D.grid)

	# restrict the analysis to type-B particles
	D.add_filter("species == 'B'", group=0)

	# compute the descriptor's data matrix
	X = D.compute()

	# print the first three feature vectors
	print("feature vectors:\n", X[0:3])
	
Output:

.. code-block:: litteral

	grid:
	 [1.05 1.15 1.25 1.35 1.45]
	feature vectors:
	 [[2 2 1 1 3]
	  [5 1 0 1 2]
	  [4 2 1 0 1]]

- ``grid`` shows the grid of distances :math:`\{ d_n \}`, where :math:`\Delta r = 0.1`.
- ``feature vectors`` shows the first three feature vectors :math:`X^\mathrm{R}(1)`, :math:`X^\mathrm{R}(2)` and :math:`X^\mathrm{R}(3)` corresponding to the grid.

References
----------

.. bibliography:: ../../references.bib
	:style: unsrt
	:filter: docname in docnames
