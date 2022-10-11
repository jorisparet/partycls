Tetrahedral descriptor
======================

Definition
----------

The degree of tetrahedrality of a particle :math:`i` is the average deviation of the bond angles :math:`\{ \theta_{jik} \}` between :math:`i` and all the possible pairs of its nearest neighbors :math:`(j,k)` from the ideal angle in a tetrahedron, :math:`\theta_\mathrm{tetra} = 109.5^\circ`:

.. math::
	T(i) = \frac{1}{N_\mathrm{ba}(i)} \sum_{j=1}^{N_b(i)} \sum_{\substack{k=1 \\ k \neq j}}^{N_b(i)} | \cos(\theta_{jik}) - \cos(\theta_\mathrm{tetra}) | ,

where :math:`N_\mathrm{ba}(i)` is the total number of bond angles (*i.e.* the number of pairs) around particle :math:`i` and :math:`N_b(i)` is the number of its nearest neighbors. The resulting feature vector for particle :math:`i` is given by

.. math::
	X^\mathrm{T}(i) = (\: T(i) \:) .

.. note::
	Unlike most descriptors, this descriptor is scalar. Its feature vector :math:`X^\mathrm{T}(i)` is thus composed of a single feature.

Setup
-----

Instantiating this descriptor on a :py:class:`Trajectory <partycls.trajectory.Trajectory>` can be done as follows:

.. code-block:: python

	from partycls import Trajectory
	from partycls.descriptor import TetrahedralDescriptor

	traj = Trajectory("trajectory.xyz")
	D = TetrahedralDescriptor(traj)

The constructor takes the following parameters:

.. automethod:: partycls.descriptor.tetrahedrality.TetrahedralDescriptor.__init__

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

We now instantiate a :py:class:`TetrahedralDescriptor <partycls.descriptor.tetrahedrality.TetrahedralDescriptor>` on this trajectory and restrict the analysis to type-B particles only:

.. code-block:: python

	from partycls.descriptor import TetrahedralDescriptor

	# instantiation
	D = TetrahedralDescriptor(traj)

	# restrict the analysis to type-B particles
	D.add_filter("species == 'B'", group=0)

	# compute the descriptor's data matrix
	X = D.compute()

	# print the first three feature vectors
	print("feature vectors:\n", X[0:3])

.. code-block:: none
	:caption: **Output:**

	feature vectors:
	 [[0.48286880]
	  [0.48912898]
	  [0.47882811]]

-  ``feature vectors`` shows the first three feature vectors :math:`X^\mathrm{T}(1)`, :math:`X^\mathrm{R}(2)` and :math:`X^\mathrm{R}(3)`.
