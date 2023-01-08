Smoothed bond-orientational descriptor
======================================

.. Important::
	See :doc:`bo` first.

Definition
----------

This is a smooth version of the :doc:`bo`, in which the coefficients :math:`q_{lm}(i)` are multiplied by a weighting function :math:`f(r)` that depends on the radial distance :math:`r` between the central particle :math:`i` and other surrounding particles :math:`j`, where :math:`j` can be any particle in the system (*i.e.* not necessarily a nearest neighbors of :math:`i`) :cite:`coslovich_2022`.

The smoothed complex coefficients are given by

.. math::
	q_{lm}^{S}(i) = \frac{1}{Z(i)} \sum_{j=1}^{N} f({r}_{ij}) Y_{lm}(\hat{\mathbf{r}}_{ij}) ,


where :math:`Z(i)=\sum_{j=1}^{N} f({r}_{ij})` is a normalization constant and the superscript :math:`S` indicates the smooth nature of the descriptor. We use

.. math::
	f(r_{ij}) = \exp \left[- (r_{ij} / r_{\alpha\beta}^c)^\gamma \right] H(R_{\alpha\beta}^c - r_{ij}) ,


where :math:`r_{\alpha\beta}^c` is the first minimum of the corresponding partial radial distribution function for the pair :math:`(i,j)` and :math:`\gamma` is an integer.
Also, :math:`H` is the `Heaviside step function <https://en.wikipedia.org/wiki/Heaviside_step_function>`_, which ensures, for efficiency reasons, that the descriptor only has contributions from particles within a distance :math:`R_{\alpha\beta}^c = \xi \times r_{\alpha\beta}^c` from the central one, where :math:`\xi > 1` is a scaling factor.

The rotational invariants are defined similarly to the :doc:`bo`.

We then consider :math:`Q_l^S(i)` for a sequence of orders :math:`\{ l_n \} = \{ l_\mathrm{min}, \dots, l_\mathrm{max} \}`. The resulting feature vector for particle :math:`i` is given by

.. math::
	X^\mathrm{SBO}(i) = (\: Q_{l_\mathrm{min}}^S(i) \;\; \dots \;\; Q_{l_\mathrm{max}}^S(i) \:) .

Setup
-----

Instantiating this descriptor on a :py:class:`Trajectory <partycls.trajectory.Trajectory>` can be done as follows:

.. code-block:: python

	from partycls import Trajectory
	from partycls.descriptors import SmoothedBondOrientationalDescriptor

	traj = Trajectory("trajectory.xyz")
	D = SmoothedBondOrientationalDescriptor(traj)

The constructor takes the following parameters:

.. automethod:: partycls.descriptors.smoothed_bo.SmoothedBondOrientationalDescriptor.__init__

Requirements
------------

The computation of this descriptor relies on:

- **Nearest neighbors cutoffs**. These can either be set in the :py:class:`Trajectory <partycls.trajectory.Trajectory>` prior to the computation of the descriptor, or computed from inside the descriptor using a default method.

Demonstration
-------------

We consider an input trajectory file :file:`trajectory.xyz` in XYZ format that contains two particle types ``"A"`` and ``"B"``. We can either compute or set the nearest neighbors cutoffs :math:`\{ r_{\alpha\beta}^c \}` for the smoothing directly in :py:class:`Trajectory <partycls.trajectory.Trajectory>`:

.. code-block:: python

	from partycls import Trajectory

	# open the trajectory
	traj = Trajectory("trajectory.xyz")

	# compute the nearest neighbors cutoffs
	traj.compute_nearest_neighbors_cutoffs(dr=0.1)
	print("computed cuttofs\n", traj.nearest_neighbors_cutoffs)

	# set the nearest neighbors cuttofs
	traj.nearest_neighbors_cutoffs = [1.45, 1.35, 1.35, 1.25]
	print("manually set cuttofs\n", traj.nearest_neighbors_cutoffs)

.. code-block:: none
	:caption: **Output:**

	computed cutoffs:
	 [1.4500000000000004, 1.3500000000000003, 1.3500000000000003, 1.2500000000000002]
	manually set cutoffs:
	 [1.45, 1.35, 1.35, 1.25]

.. note::
	If not computed in :py:class:`Trajectory <partycls.trajectory.Trajectory>` or manually set, the cutoffs :math:`\{ r_{\alpha\beta}^c \}` will be computed from inside the descriptor.

We now instantiate a :py:class:`SmoothedBondOrientationalDescriptor <partycls.descriptors.smoothed_bo.SmoothedBondOrientationalDescriptor>` on this trajectory and restrict the analysis to type-B particles only. We set the grid of orders :math:`\{l_n\} = \{2,4,6,8\}`, :math:`\xi=1.3` and :math:`\gamma=8`:

.. code-block:: python

	from partycls.descriptors import SmoothedBondOrientationalDescriptor

	# instantiation
	D = SmoothedBondOrientationalDescriptor(traj,
						orders=[2,4,6,8],
						cutoff_enlargement=1.3,
						exponent=8)

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
	 [[0.03156284 0.11095233 0.4112718  0.23829355]
	  [0.0698711  0.08107918 0.47678647 0.26671868]
	  [0.06221017 0.09806095 0.39152213 0.16630718]]


- ``grid`` shows the grid of orders :math:`\{ l_n \}`.
- ``feature vectors`` shows the first three feature vectors :math:`X^\mathrm{SBO}(1)`, :math:`X^\mathrm{SBO}(2)` and :math:`X^\mathrm{SBO}(3)` corresponding to the grid.

References
----------

.. bibliography:: ../../references.bib
	:style: unsrt
	:filter: docname in docnames
