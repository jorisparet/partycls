Smoothed bond-angle descriptor
==============================

.. Important::
	See :doc:`ba` first.

Definition
----------

This is a smooth version of the :doc:`ba` in which the bond angles :math:`\theta_{jik}` are multiplied by a weighting function :math:`f(r_{ij}, r_{ik})` that depends on the radial distances :math:`r_{ij}` and :math:`r_{ik}` between :math:`(i,j)` and :math:`(i,k)` respectively, where :math:`j` and :math:`k` can be any particle in the system (*i.e.* not necessarily a nearest neighbors of :math:`i`).

Essentially, we define the (non-integer) number of bond angles around particle :math:`i` as

.. math::
	N_i^S(\theta_n) = \sum_{j=1}^N \sum_{\substack{k=1 \\ k \neq j}}^N f(r_{ij}, r_{ik}) \delta(\theta_n - \theta_{jik}) ,

where the superscript :math:`S` indicates the smooth nature of the descriptor, and :math:`f(r_{ij}, r_{ik})` is the following smoothing function:

.. math::
	f(r_{ij}, r_{ik}) = \exp \left[ - \left[ ( r_{ij} / r_{\alpha\beta}^c )^\gamma + ( r_{ik} / r_{\alpha\beta'}^c )^\gamma \right] \right] H( R_\mathrm{max}^c - r_{ij} ) H( R_\mathrm{max}^c - r_{ik}) ,

where:

- :math:`r_{\alpha\beta}^c` and :math:`r_{\alpha\beta'}^c` are the first minima of the corresponding partial radial distribution functions for the pairs :math:`(i,j)` and :math:`(i,k)`.
- :math:`\gamma` is an integer.
-  :math:`R_\mathrm{max}^c = \xi \times \max(\{ r_{\alpha\beta}^c \})` is the largest nearest neighbor cutoff rescaled by :math:`\xi > 1`.
- :math:`H` is the `Heavide step function <https://en.wikipedia.org/wiki/Heaviside_step_function>`_, which ensures for efficiency reasons, that the descriptor only has contributions from particles within a distance :math:`R_\mathrm{max}^c`.

We then consider :math:`N_i^S(\theta_n)` for a set of angles :math:`\{ \theta_n \}` that go from :math:`\theta_0 = 0^\circ` to :math:`\theta_{n_\mathrm{max}}=180^\circ` by steps of :math:`\Delta \theta`. The resulting feature vector for particle :math:`i` is given by

.. math::
	X^\mathrm{SBA}(i) = (\: N_i^S(\theta_0) \;\; N_i^S(\theta_1) \;\; \dots \;\; N_i^S(\theta_{n_\mathrm{max}}) \:) .

Setup
-----

Instantiating this descriptor on a ``Trajectory`` can be done as follows:

.. code-block:: python

	from partycls import Trajectory
	from partycls.descriptor import SmoothedBondAngleDescriptor

	traj = Trajectory("trajectory.xyz")
	D = SmoothedBondAngleDescriptor(traj)

The constructor takes the following parameters:

.. automethod:: partycls.descriptor.smoothed_ba.SmoothedBondAngleDescriptor.__init__

Demonstration
-------------

We consider an input trajectory file ``"trajectory.xyz"`` in XYZ format that contains two particle types ``"A"`` and ``"B"``:

.. code-block:: python

	from partycls import Trajectory

	traj = Trajectory("trajectory.xyz")

We compute the nearest neighbors cutoffs :math:`\{ r_{\alpha\beta}^c \}` for the smoothing:

.. code-block:: python

	traj.compute_nearest_neighbors_cutoffs(dr=0.1)
	print("computed cutoffs:\n", traj.nearest_neighbors_cutoffs)

Output:

.. code-block:: litteral

	computed cutoffs:
	 [1.45, 1.35, 1.35, 1.25]
	manually set cutoffs:
	 [1.45, 1.35, 1.35, 1.25]


Alternatively, we can set the cutoffs manually if they are known:

.. code-block:: python

	traj.nearest_neighbors_cutoffs = [1.45, 1.35, 1.35, 1.25]
	print("manually set cutoffs:\n", traj.nearest_neighbors_cutoffs)

Output:

.. code-block:: litteral

	manually set cutoffs:
	 [1.45, 1.35, 1.35, 1.25]

.. note::
	If not computed in the ``Trajectory`` or manually set, the cutoffs will be computed from inside the descriptor.

We now instantiate and compute a ``SmoothedBondAngleDescriptor`` on this trajectory:

.. code-block:: python

	from partycls.descriptor import SmoothedBondAngleDescriptor

	# instantiation
	D = SmoothedBondAngleDescriptor(traj,
					dtheta=18.0,
					cutoff_enlargement=1.3,
					exponent=8)

	# print the grid of angles (in degrees)
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
	 [  9.  27.  45.  63.  81.  99. 117. 135. 153. 171.]
	feature vectors:
	 [[ 0.          0.24735055  5.1652519  29.43498845 10.5325834  14.99235213
	  19.81940987 10.74915154  5.74995792  3.83545611]
	  [ 0.          0.16020613  4.79852719 35.17585892  9.27868908 14.30365693
	  20.88630866 12.92153832  2.269351    7.38748952]
	  [ 0.          0.08214317 11.23967682 32.2093987   4.0642088  24.10157113
	  19.94955473  7.72183504 12.2267004   3.29940419]]
