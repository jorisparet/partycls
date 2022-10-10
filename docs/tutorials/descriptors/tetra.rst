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

Instantiating this descriptor on a ``Trajectory`` can be done as follows:

.. code-block:: python

	from partycls import Trajectory
	from partycls.descriptor import TetrahedralDescriptor

	traj = Trajectory("trajectory.xyz")
	D = TetrahedralDescriptor(traj)

The constructor takes the following parameters:

.. automethod:: partycls.descriptor.tetrahedrality.TetrahedralDescriptor.__init__

Demonstration
-------------

We consider an input trajectory file ``"trajectory.xyz"`` in XYZ format that contains two particle types ``"A"`` and ``"B"``:

.. code-block:: python

	from partycls import Trajectory

	# open the trajectory
	traj = Trajectory("trajectory.xyz")

	# compute the neighbors using Voronoi tessellation
	traj.compute_nearest_neighbors(method='voronoi')
	nearest_neighbors = traj.get_property("nearest_neighbors")
	
	# print the first three neighbors lists for the first trajectory frame
	print("neighbors:\n",nearest_neighbors[0][0:3])

Output:

.. code-block:: litteral

	neighbors:
	 [list([767, 113, 323, 276, 322, 332, 980, 425, 801, 16, 171, 258, 901, 241])
	  list([448, 481, 496, 574, 706, 536, 337, 241, 951, 16, 14, 258, 502, 447, 860])
	  list([799, 123, 608, 913, 500, 826, 230, 636, 796, 772, 810, 639, 270, 578, 874, 397, 354])]

We now instantiate and compute a ``TetrahedralDescriptor`` on this trajectory:

.. code-block:: python

	from partycls.descriptor import TetrahedralDescriptor

	# instantiation
	D = TetrahedralDescriptor(traj, dtheta=3.0)

	# restrict the analysis to type-B particles
	D.add_filter("species == 'B'", group=0)

	# compute the descriptor's data matrix
	X = D.compute()

	# print the first three feature vectors
	print("feature vectors:\n", X[0:3])
	
Output:

.. code-block:: litteral

	feature vectors:
	 [[0.48286880]
	  [0.48912898]
	  [0.47882811]]
