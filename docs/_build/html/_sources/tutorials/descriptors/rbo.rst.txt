Radial bond-orientational descriptor
====================================

.. Important::
	See :doc:`bo` first.

Definition
----------

The radial bond-order descriptor captures the radial dependence of bond order in the most straightforward way :cite:`boattini_2021`. It uses a radial basis of Gaussian functions of width :math:`\delta` centered on a grid of points :math:`\{d_n\}_{n=1 \dots n_\mathrm{max}}`,

.. math::
	G_n(r) = \exp{\left(-\frac{(d_n - r)^2}{2\delta^2}\right)} .

The complex radial bond-order coefficients are defined as

.. math::
  q_{l m n}^{R}(i) = \frac{1}{Z(i)} \sum_{j=1}^{N}
  G_n(r_{ij}) Y_{l m}(\hat{\mathbf{r}}_{ij}) ,

where :math:`Z(i) = \sum_{j=1}^N G_n(r_{ij})` is a normalization constant and the superscript :math:`R` indicates the radial dependence of the descriptor.
In the following, we actually use

.. math::
	G_n(r_{ij}) = \exp{\left(-\frac{(d_n - r_{ij})^\gamma}{2\delta^2}\right)} H(R_\mathrm{max} - r_{ij}) ,

where :math:`\gamma` is an integer, :math:`H` is the `Heaviside step function <https://en.wikipedia.org/wiki/Heaviside_step_function>`_, which allows to neglect the contributions of particles further than a distance :math:`R_\mathrm{max} = d_{n_\mathrm{max}} + \sigma \times \delta` from the central particle, where :math:`d_{n_\mathrm{max}}` is the largest distance in the grid of points :math:`\{ d_n \}` and :math:`\sigma` is a skin distance.
Then, only the diagonal coefficients of the power spectrum, namely 

.. math::
	Q_{l,n}^R(i) = \sqrt{ \frac{4\pi}{2l + 1} \sum_{m=-l}^l |q_{l m n}^R(i)|^2 } ,

are retained to form the descriptor of particle :math:`i`.

We then consider :math:`Q^R_{l,n}(i)` for a sequence of orders :math:`\{ l_m \} = \{ l_\mathrm{min}, \dots, l_\mathrm{max} \}` and for a grid of distances :math:`\{ d_n \}`. The resulting feature vector for particle :math:`i` is given by

.. math::
	X^\mathrm{RBO}(i) = (\: \dots \;\; Q^R_{l,n}(i) \;\; Q^R_{l, n+1}(i) \;\; \dots \;\; Q^R_{l_\mathrm{max}, n_\mathrm{max}}(i) \:) .

Setup
-----

Instantiating this descriptor on a :py:class:`Trajectory <partycls.trajectory.Trajectory>` can be done as follows:

.. code-block:: python

	from partycls import Trajectory
	from partycls.descriptors import RadialBondOrientationalDescriptor

	traj = Trajectory("trajectory.xyz")
	D = RadialBondOrientationalDescriptor(traj)

The constructor takes the following parameters:

.. automethod:: partycls.descriptors.radial_bo.RadialBondOrientationalDescriptor.__init__

.. hint::

	The alias :py:class:`BoattiniDescriptor <partycls.descriptors.radial_bo.BoattiniDescriptor>` can be used in place of :py:class:`RadialBondOrientationalDescriptor <partycls.descriptors.radial_bo.RadialBondOrientationalDescriptor>`.

Demonstration
-------------

We consider an input trajectory file :file:`trajectory.xyz` in XYZ format that contains two particle types ``"A"`` and ``"B"``:

.. code-block:: python

	from partycls import Trajectory

	# open the trajectory
	traj = Trajectory("trajectory.xyz")

We now instantiate a :py:class:`RadialBondOrientationalDescriptor <partycls.descriptors.radial_bo.RadialBondOrientationalDescriptor>` on this trajectory and restrict the analysis to type-B particles only. We set set the grid of orders :math:`\{l_m\} = \{4,6\}`, the grid of distances :math:`\{d_n\} = \{1.0, 1.25, 1.50\}`, the Gaussian shell width :math:`\delta=0.2` and the skin width :math:`\sigma = 2.5`:

.. code-block:: python

	from partycls.descriptors import RadialBondOrientationalDescriptor

	# instantiation
	D = RadialBondOrientationalDescriptor(traj,
					      orders=[4,6],
					      distance_grid=[1.0, 1.25, 1.50],
					      delta=0.2,
					      skin=2.5)

	# print the grid of orders
	print("orders:\n", D.orders)
	# print the grid of distances
	print("distances:\n", D.distance_grid)
	# print the mixed grid of orders and distances
	print("mixed grid:\n", D.grid)

	# restrict the analysis to type-B particles
	D.add_filter("species == 'B'", group=0)

	# compute the descriptor's data matrix
	X = D.compute()

	# print the first three feature vectors
	print("feature vectors:\n", X[0:3])

.. code-block:: none
	:caption: **Output:**

	orders:
	 [4 6]
	distances:
	 [1.   1.25 1.5 ]
	mixed grid:
	 [(4, 1.0), (4, 1.25), (4, 1.5), (6, 1.0), (6, 1.25), (6, 1.5)]
	feature vectors:
	 [[0.10876497 0.08802174 0.08958731 0.41008552 0.18769869 0.17186267]
	  [0.09036718 0.05933577 0.07556535 0.49122192 0.2396785  0.17942327]
	  [0.09049389 0.11888597 0.06950139 0.3981864  0.18440161 0.18568249]]

- ``orders`` shows the grid of orders :math:`\{ l_m \}`.
- ``distances`` shows the grid of distances :math:`\{ d_n \}`.
- ``mixed grid`` shows the mixed grid of orders and distances.
- ``feature vectors`` shows the first three feature vectors :math:`X^\mathrm{RBO}(1)`, :math:`X^\mathrm{RBO}(2)` and :math:`X^\mathrm{RBO}(3)` corresponding to the (mixed) grid.

References
----------

.. bibliography:: ../../references.bib
	:style: unsrt
	:filter: docname in docnames
