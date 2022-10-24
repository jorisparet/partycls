Coordination descriptor
=======================

Definition
----------

The coordination number :math:`n_\alpha(i)` of a particle :math:`i` is given by the number of its nearest neighbors whose chemical species is :math:`\alpha`, where  :math:`\alpha` is either one of the :math:`n` chemical species :math:`\{ \alpha_i \}_{i=1 \dots n}` in the trajectory (*i.e.* partial coordination number, :math:`\alpha = \alpha_i`) or all species at once (*i.e.* total coordination number, :math:`\alpha = \mathrm{all}`).

The resulting **full** feature vector for particle :math:`i` is given by

.. math::
	X^\mathrm{N}(i) = (\: n_\mathrm{all}(i) \;\;  n_{\alpha_1}(i) \;\; \dots \;\; n_{\alpha_n}(i) \:) ,

but its size depends on whether the user requests the total coordination number, the partial ones, or both.

.. note::

	By applying a filter on ``group=1``, the total/partial coordination number(s) can also be computed by considering and **arbitrary** subset of particles (*e.g.* particles whose radius is smaller than a certain value).

Setup
-----

Instantiating this descriptor on a :py:class:`Trajectory <partycls.trajectory.Trajectory>` can be done as follows:

.. code-block:: python

	from partycls import Trajectory
	from partycls.descriptors import CoordinationDescriptor

	traj = Trajectory("trajectory.xyz")
	D = CoordinationDescriptor(traj)

The constructor takes the following parameters:

.. automethod:: partycls.descriptors.coordination.CoordinationDescriptor.__init__

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

We now instantiate a :py:class:`CoordinationDescriptor <partycls.descriptors.coordination.CoordinationDescriptor>` on this trajectory, restrict the analysis to type-B particles only, and compute both the total and partial coordination numbers:

.. code-block:: python

	from partycls.descriptors import CoordinationDescriptor

	# instantiation
	D = CoordinationDescriptor(traj, total=True, partial=True)

	# print the grid of species for coordination numbers
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
	 ['all' 'A' 'B']
	feature vectors:
	 [[12  5  7]
	  [12  5  7]
	  [13  7  6]]

- ``grid`` shows the grid of considered species :math:`\alpha`. Here, we chose total and partial coordination numbers.
-  ``feature vectors`` shows the first three feature vectors :math:`X^\mathrm{N}(1)`, :math:`X^\mathrm{N}(2)` and :math:`X^\mathrm{N}(3)`. The sum of columns 2 and 3 must equal column 1.
