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

Instantiating this descriptor on a ``Trajectory`` can be done as follows:

.. code-block:: python

	from partycls import Trajectory
	from partycls.descriptor import BondAngleDescriptor

	traj = Trajectory("trajectory.xyz")
	D = BondAngleDescriptor(traj)

The constructor takes the following parameters:

.. automethod:: partycls.descriptor.ba.BondAngleDescriptor.__init__

Demonstration
-------------

We consider an input trajectory file ``"trajectory.xyz"`` in XYZ format that contains two particle types ``"A"`` and ``"B"``:

.. code-block:: python

	from partycls import Trajectory

	traj = Trajectory("trajectory.xyz")

We now instantiate and compute a ``BondAngleDescriptor`` on this trajectory:

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
	
Output:

.. code-block:: litteral

	grid:
	 [  9.  27.  45.  63.  81.  99. 117. 135. 153. 171.]
	feature vectors:
	 [[ 0  0  4 44 12 18 28 14  6  6]
	  [ 0  0  6 44 12 16 26 16  2 10]
	  [ 0  0 16 42  6 34 26 10 18  4]]

References
----------

.. bibliography:: ../../references.bib
	:style: unsrt
	:filter: docname in docnames
