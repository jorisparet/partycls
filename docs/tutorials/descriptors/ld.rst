Lechner-Dellago descriptor
==========================

.. Important::
	See :doc:`bo` first.

Definition
----------

The complex coefficients :math:`q_{lm}(i)` of particle :math:`i` can be averaged over its :math:`N_b(i)` nearest neighbors, as suggested by Lechner and Dellago :cite:`lechner_2008`,

.. math::
	\bar{q}_{lm}(i) = \frac{1}{N_b(i)+1} \left[ q_{l m}(i) + \sum_{j=1}^{N_b(i)} q_{l m}(j) \right],

and then made invariant,

.. math::
	\bar{Q}_{l}(i) = \sqrt{ \frac{4\pi}{2l + 1}\sum_{m=-l}^l |\bar{q}_{lm(i)}|^2 } ,

to provide an improved descriptor for crystal structure detection.

We then consider :math:`\bar{Q}_l(i)` for a sequence of orders :math:`\{ l_n \} = \{ l_\mathrm{min}, \dots, l_\mathrm{max} \}`. The resulting feature vector for particle :math:`i` is given by

.. math::
	X^\mathrm{LD}(i) = (\: \bar{Q}_{l_\mathrm{min}}(i) \;\; \dots \;\; \bar{Q}_{l_\mathrm{max}}(i) \:) .

Setup
-----

Instantiating this descriptor on a :py:class:`Trajectory <partycls.trajectory.Trajectory>` can be done as follows:

.. code-block:: python

	from partycls import Trajectory
	from partycls.descriptor import LechnerDellagoDescriptor

	traj = Trajectory("trajectory.xyz")
	D = LechnerDellagoDescriptor(traj)

The constructor takes the following parameters:

.. automethod:: partycls.descriptor.bo.LechnerDellagoDescriptor.__init__

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

We now instantiate a :py:class:`LechnerDellagoDescriptor <partycls.descriptor.bo.LechnerDellagoDescriptor>` on this trajectory and restrict the analysis to type-B particles only. We set set the grid of orders :math:`\{l_n\} = \{2,4,6,8\}`:

.. code-block:: python

	from partycls.descriptor import LechnerDellagoDescriptor

	# instantiation
	D = LechnerDellagoDescriptor(traj, orders=[2,4,6,8])

	# print the grid of orders
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
	 [2 4 6 8]
	feature vectors:
	 [[0.03366521 0.04034078 0.08648078 0.1120834 ]
	  [0.01483751 0.03889963 0.16849717 0.11150705]
	  [0.02312734 0.02640117 0.11722934 0.11053876]]

- ``grid`` shows the grid of orders :math:`\{ l_n \}`.
- ``feature vectors`` shows the first three feature vectors :math:`X^\mathrm{LD}(1)`, :math:`X^\mathrm{LD}(2)` and :math:`X^\mathrm{LD}(3)` corresponding to the grid.

References
----------

.. bibliography:: ../../references.bib
	:style: unsrt
	:filter: docname in docnames
