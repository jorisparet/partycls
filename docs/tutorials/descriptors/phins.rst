Physics-inspired descriptor
===========================

Definition
----------

The physics-inspired descriptor is a collection of established structural order parameters (observables), first introduced by *Jung et al.* :cite:`jung_2023` to train *GlassMLP*, a deep neural network used to predict the long-time dynamics in deeply supercooled liquids.

For each particle :math:`i`, four distinct observables are computed and coarse-grained over a set of lengthscales :math:`\{L_n\}`:

1. Local density: 

.. math::
	\bar{\rho}_i(L_n;\beta_m) = \sum_{j \in N_{\beta_m}^i} \exp(-\frac{r_{ij}}{L_n}) ,

where the sum runs over all :math:`N_{\beta_m}^i` particles of type :math:`\beta_m` within distance :math:`r_{ij} < R_c` from particle :math:`i` (including :math:`i` itself if its type is :math:`\beta_m`). Particles :math:`j` beyond the cutoff distance :math:`R_c` will be ignored for efficiency reasons.

2. Potential energy:

.. math::
	\bar{E}_i(L_n;\beta_m) = \frac{1}{\bar{\rho}_i(L_n;\beta_m)} \sum_{j \in N_{\beta_m}^i} E_j \exp(-\frac{r_{ij}}{L_n}) ,
	
where :math:`E_j` is extracted from the pair potential :math:`V(r)`, such that :math:`E_j = \sum_{j \neq k} V(r_{jk})`.

3. Perimeter or surface of the Voronoi cell centered on particle :math:`i` (2- and 3-dimensional cases, respectively):

.. math::
	\bar{\mathcal{C}}_i(L_n;\beta_m) = \frac{1}{\bar{\rho}_i(L_n;\beta_m)} \sum_{j \in N_{\beta_m}^i} \mathcal{C}_j \exp(- \frac{r_{ij}}{L_n}) .

4. Local variance of the potential energy:

.. math::
	\mathrm{Var}( \bar{E}_i(L_n;\beta_m) ) = \frac{1}{\bar{\rho}_i(L_n;\beta_m)} \sum_{j \in N_{\beta_m}^i} \left( E_j - \bar{E}_i(L_n;\beta_m) \right)^2 \exp(-\frac{r_{ij}}{L_n}) .

The prefactor :math:`\frac{1}{\bar{\rho}_i(L_n;\beta_m)}` is used for normalization purposes. These four observables are then arranged in

.. math::
	\bar{S}_i(L_n;\beta_m) = \{ \bar{\rho}_i(L_n;\beta_m), \bar{E}_i(L_n;\beta_m), \bar{\mathcal{C}}_i(L_n;\beta_m),  \mathrm{Var}( \bar{E}_i(L_n;\beta_m) )\}
	
and computed for each distance in :math:`\{ L_n \}` and for each particle type in :math:`\{ \beta_m \}`, including one additional computation where **all** particles are considered at once regardless of their type. The total number :math:`M` of features produced is given by :math:`M = 4 \times N_L \times (N_\beta + 1)`, where :math:`N_L` is the total number of coarse-graining lengths and :math:`N_\beta` is the total number of distinct particle types in the system.

The resulting feature vector for particle :math:`i` is given by

.. math::
	X^\mathrm{P}(i) = (\: \bar{S}_i(L_1;\beta_1) \;\; \dots \;\; \bar{S}_i(L_n;\beta_m) \;\; \dots \;\;  \bar{S}_i(L_{N_L};\beta_{N_{\beta}+1} ) \:) .

Setup
-----

.. warning::
	This descriptor uses particles' properties for the potential energy :math:`E_i` and the Voronoi cell metric :math:`\mathcal{C}_i`. They must either be read from the input trajectory file using the ``additional_fields`` parameter or set using the method :py:class:`Trajectory.set_property <partycls.trajectory.Trajectory.set_property>`. They are later identified in the descriptor by the ``energy_field`` and ``cell_field`` parameters respectively.

Instantiating this descriptor on a :py:class:`Trajectory <partycls.trajectory.Trajectory>` can be done as follows:

.. code-block:: python

	from partycls import Trajectory
	from partycls.descriptors import PhysicsInspiredDescriptor

	traj = Trajectory("trajectory.xyz")
	D = PhysicsInspiredDescriptor(traj)

The constructor takes the following parameters:

.. automethod:: partycls.descriptors.physics_inspired.PhysicsInspiredDescriptor.__init__

.. hint::

	The alias :py:class:`JungDescriptor <partycls.descriptors.physics_inspired.JungDescriptor>` can be used in place of :py:class:`PhysicsInspiredDescriptor <partycls.descriptors.physics_inspired.PhysicsInspiredDescriptor>`.

Demonstration
-------------

We consider an input trajectory file :file:`trajectory.xyz` in XYZ format of a **2D system** that contains two particle types ``"A"`` and ``"B"``. The file contains two additional columns ``energy`` and ``perimeter`` that are used for the potential energy :math:`E_i` and the Voronoi cell perimeter :math:`C_i` of each particle, respectively.

.. code-block:: python

	from partycls import Trajectory

	# open the trajectory
	traj = Trajectory("trajectory.xyz", additional_fields=['energy', 'perimeter'])


We now instantiate a :py:class:`PhysicsInspiredDescriptor <partycls.descriptors.physics_inspired.PhysicsInspiredDescriptor>` on this trajectory and manually set the coarse-graining lengths :math:`\{ L_n \}` with the ``distance_grid`` parameter, the cutoff :math:`R_c` with the ``cutoff`` parameter, and set the fields used to store :math:`E_i` and :math:`C_i` with the parameters ``energy_field`` and ``cell_field`` respectively. We restrict the analysis to type-B particles only:

.. code-block:: python

	from partycls.descriptors import PhysicsInspiredDescriptor

	# instantiation
	D = PhysicsInspiredDescriptor(traj,
				       distance_grid=[0.0, 1.0, 2.0],
				       cutoff=10.0,
				       energy_field="energy",
				       cell_metric="perimeter")

	# print the grid of coarse-graining distances
	print("distances:\n", D.distance_grid)
	# print the grid of distances, particle types and observables
	print("grid of features:\n", D.grid)
	# print the total number of features (cross-checked with the expression for M)
	print("number of features:\n",
	      D.n_features,
	      '\n',
	      4 * len(D.distance_grid) * (len(traj[0].distinct_species) + 1))

	# restrict the analysis to type-B particles
	D.add_filter("species == 'B'", group=0)

	# compute the descriptor's data matrix
	X = D.compute()

	# print the first three feature vectors
	print("feature vectors:\n", X[0:3])

.. code-block:: none
	:caption: **Output:**
	
	distances:
	 [0. 1. 2.]
	grid of features:
	 [('length=0.0', 'type=A', 'obs=DEN'),
	  ('length=0.0', 'type=A', 'obs=EPOT'),
	  ('length=0.0', 'type=A', 'obs=CELL'),
	  ('length=0.0', 'type=A', 'obs=VAR_EPOT'),
	  ...
	 ]
	number of features:
	 48
	 48
	feature vectors:
	 ...
	 
- ``distances`` shows the coarse-graining lengths :math:`\{ L_n \}`.
- ``grid of features`` shows all the combinations of coarse-graining lengths, particle types to compute each of the four coarse-grained observables. ``DEN`` refers to local density, ``EPOT`` to local potential energy, ``CELL`` to the Voronoi cell metric, and ``VAR_EPOT`` to the variance of the potential energy.
- ``number of features`` shows the total number of features :math:`M`.
-  ``feature vectors`` shows the first three feature vectors :math:`X^\mathrm{P}(1)`, :math:`X^\mathrm{P}(2)` and :math:`X^\mathrm{P}(3)`.

References
----------

.. bibliography:: ../../references.bib
	:style: unsrt
	:filter: docname in docnames
