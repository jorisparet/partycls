Bond-angle descriptor
=====================

Definition
----------

The *empirical* distribution of bond angles :math:`N_i(\theta_n)` around a central particle :math:`i` is obtained by counting, for all the possible pairs of it nearest neighbors :math:`(j,k)`, the number of bond angles :math:`\theta_{jik}` between :math:`\theta_n = n \times \Delta \theta` and :math:`\theta_{n+1} = (n+1) \times \Delta \theta`, where :math:`\Delta \theta` has the interpration of a bin width in a histogram :cite:`paret_2020`.

We then consider :math:`N_i(\theta_n)` for a set of angles :math:`\{ \theta_n \}` that go from :math:`\theta_0 = 0^\circ` to :math:`\theta_{n_\mathrm{max}}=180^\circ` by steps of :math:`\Delta \theta`. The resulting feature vector for particle :math:`i` is given by

.. math::
	X^\mathrm{BA}(i) = (\: N_i(\theta_0) \;\; N_i(\theta_1) \;\; \dots \;\; N_i(\theta_{n_\mathrm{max}}) \:) .

Setup
-----

Instantiating this descriptor on a :py:class:`Trajectory <partycls.trajectory.Trajectory>` can be done as follows:

.. code-block:: python

	from partycls import Trajectory
	from partycls.descriptor import BondAngleDescriptor

	traj = Trajectory("trajectory.xyz")
	D = BondAngleDescriptor(traj)

The constructor takes the following parameters:

.. automethod:: partycls.descriptor.ba.BondAngleDescriptor.__init__

Requirements
------------

The computation of this descriptor relies on:

- **Lists of nearest neighbors**. These can either be read from the input trajectory file, computed in the :py:class:`Trajectory <partycls.trajectory.Trajectory>`, or computed from inside the descriptor using a default method.

Demonstration
-------------

We consider an input trajectory file :file:`trajectory.xyz` in XYZ format that contains two particle types ``"A"`` and ``"B"``. We compute the lists of nearest neighbors using the fixed-cutoffs method:

.. code-block:: python

	from partycls import Trajectory

	# open the trajectory
	traj = Trajectory("trajectory.xyz")

	# compute the neighbors using pre-computed cuttofs
	traj.nearest_neighbors_cuttofs = [1.45, 1.35, 1.35, 1.25]
	traj.compute_nearest_neighbors(method='fixed')
	nearest_neighbors = traj.get_property("nearest_neighbors")
	
	# print the first three neighbors lists for the first trajectory frame
	print("neighbors:\n",nearest_neighbors[0][0:3])

.. code-block:: none
	:caption: **Output:**

	neighbors:
	 [list([16, 113, 171, 241, 258, 276, 322, 323, 332, 425, 767, 801, 901, 980])
	  list([14, 241, 337, 447, 448, 481, 496, 502, 536, 574, 706, 860, 951])
	  list([123, 230, 270, 354, 500, 578, 608, 636, 639, 640, 796, 799, 810, 826, 874, 913])]

We now instantiate a :py:class:`BondAngleDescriptor <partycls.descriptor.ba.BondAngleDescriptor>` on this trajectory and restrict the analysis to type-B particles only. We set :math:`\Delta \theta = 18^\circ`:

.. code-block:: python

	from partycls.descriptor import BondAngleDescriptor

	# instantiation
	D = BondAngleDescriptor(traj, dtheta=18.0)

	# print the grid of angles (in degrees)
	print("grid:\n", D.grid)

	# restrict the analysis to type-B particles
	D.add_filter("species == 'B'", group=0)

	# compute the descriptor's data matrix
	X = D.compute()

	# print the first three feature vectors
	print("feature vectors:\n", X[0:3])

.. code-block:: none
	:caption: **Output:**

	grid:
	 [  9.  27.  45.  63.  81.  99. 117. 135. 153. 171.]
	feature vectors:
	 [[ 0  0  4 44 12 18 28 14  6  6]
	  [ 0  0  6 44 12 16 26 16  2 10]
	  [ 0  0 16 42  6 34 26 10 18  4]]

- ``grid`` shows the grid of angles :math:`\{ \theta_n \}` in degrees, where :math:`\Delta \theta = 18^\circ`.
- ``feature vectors`` shows the first three feature vectors :math:`X^\mathrm{BA}(1)`, :math:`X^\mathrm{BA}(2)` and :math:`X^\mathrm{BA}(3)` corresponding to the grid.

References
----------

.. bibliography:: ../../references.bib
	:style: unsrt
	:filter: docname in docnames
