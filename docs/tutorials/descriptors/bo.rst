Bond-orientational descriptor
=============================

Definition
----------

Bond-order parameters :cite:`steinhardt_1983` are standard measures of structure in the first coordination shell. Let :math:`\mathbf{r}_i` be the position of particle :math:`i` and define :math:`\mathbf{r}_{ij} = \mathbf{r}_j - \mathbf{r}_i` and :math:`r_{ij} = |\mathbf{r}_{ij}|`. Then consider the weighted microscopic density around particle :math:`i`:

.. math::
	\rho(\mathbf{r}; i) = \sum_{j=1}^{N_b(i)} w_j \delta(\mathbf{r} - \mathbf{r}_{ij})

where :math:`w_j` is a particle-dependent weight and the sum involves a set of :math:`N_b(i)` particles, which defines the coordination shell of interest for particle :math:`i`.

We project the microscopic density on a unit-radius sphere, that is, :math:`\hat{\rho}(\hat{\mathbf{r}}; i) = \sum_{j=1}^{N_b(i)} w_j \delta(\mathbf{r} - \hat{\mathbf{r}}_{ij})`,
where :math:`\hat{\mathbf{r}} = \mathbf{r} / |\mathbf{r}|` and similarly :math:`\hat{\mathbf{r}}_{ij} = \mathbf{r}_{ij}/|\mathbf{r}_{ij}|`. Expanding in spherical harmonics yields

.. math::
	\hat{\rho}(\hat{\mathbf{r}}; i) = \sum_{l=0}^\infty \sum_{m=-l}^l c_{l m}(i) Y_{l m}(\hat{\mathbf{r}}) ,

with coefficients

.. math::
	c_{l m}(i) =  \int d\mathbf{r} \rho(\mathbf{r}; i) Y_{l m}(\hat{\mathbf{r}}) .

In the conventional bond-order analysis, one sets the weights :math:`w_j` to unity and considers the normalized complex coefficients,

.. math::
	\begin{align}
	q_{lm}(i) & = \frac{1}{N_b(i)} \int d\mathbf{r} \rho(\mathbf{r}; i) Y_{l m}(\hat{\mathbf{r}}) 
	\nonumber \\ & = \frac{1}{N_b(i)} \sum_{j=1}^{N_b(i)} Y_{l m}(\hat{\mathbf{r}}_{ij}) .
	\end{align}

The rotational invariants,

.. math::
	Q_{l}(i) = \sqrt{ \frac{4\pi}{2l + 1}\sum_{m=-l}^l |q_{lm}(i)|^2 },

provide a detailed structural description of the local environment around particle :math:`i`.


We then consider :math:`Q_l(i)` for a sequence of orders :math:`\{ l_n \} = \{ l_\mathrm{min}, \dots, l_\mathrm{max} \}`. The resulting feature vector for particle :math:`i` is given by

.. math::
	X^\mathrm{BO}(i) = (\: Q_{l_\mathrm{min}}(i) \;\; \dots \;\; Q_{l_\mathrm{max}}(i) \:) .

Setup
-----

Instantiating this descriptor on a :py:class:`Trajectory <partycls.trajectory.Trajectory>` can be done as follows:

.. code-block:: python

	from partycls import Trajectory
	from partycls.descriptors import BondOrientationalDescriptor

	traj = Trajectory("trajectory.xyz")
	D = BondOrientationalDescriptor(traj)

The constructor takes the following parameters:

.. automethod:: partycls.descriptors.bo.BondOrientationalDescriptor.__init__

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

We now instantiate a :py:class:`BondOrientationalDescriptor <partycls.descriptors.bo.BondOrientationalDescriptor>` on this trajectory and restrict the analysis to type-B particles only. We set set the grid of orders :math:`\{l_n\} = \{2,4,6,8\}`:

.. code-block:: python

	from partycls.descriptors import BondOrientationalDescriptor

	# instantiation
	D = BondOrientationalDescriptor(traj, orders=[2,4,6,8])

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
	 [[0.06498973 0.10586717 0.46374576 0.22207796]
	  [0.12762569 0.09640384 0.49318559 0.29457554]
	  [0.08327171 0.11151433 0.37917788 0.17902556]]

- ``grid`` shows the grid of orders :math:`\{ l_n \}`.
- ``feature vectors`` shows the first three feature vectors :math:`X^\mathrm{BO}(1)`, :math:`X^\mathrm{BO}(2)` and :math:`X^\mathrm{BO}(3)` corresponding to the grid.

References
----------

.. bibliography:: ../../references.bib
	:style: unsrt
	:filter: docname in docnames
