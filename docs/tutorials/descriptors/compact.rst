Compactness descriptor
======================

Definition
----------

The compactness descriptor quantifies the local packing efficiency of a particle :math:`i` by comparing it to a reference ideal packing configuration of its nearest neighbors :cite:`tong_2018`.

We consider a central particle :math:`i` surrounded by :math:`N_b(i)` nearest neighbors, usually identified by means of a radical Voronoi tessellation :cite:`gellatly_1982`. We then consider triplets :math:`(j,k,m)` of neighboring particles, for which all particles are simultaneously nearest neighbors of each other and of particle :math:`i`. Such a triplet is identified with the central particle :math:`i` to be a tetrahedron.

Particles in the tetrahedron :math:`\langle ijkm \rangle` have radii :math:`(r_i, r_j, r_k, r_m)` respectively, and the triplet :math:`(j,k,m)` are at distances :math:`(r_{ij},r_{ik},r_{ik})` from the central particle :math:`i`, which are the lengths of each edge of the tetrahedron :math:`\langle ijkm \rangle`.

The *reference* tetrahedron for these four particles is the configuration in which they are all perfectly in touch, *i.e.* the edge lengths :math:`(\sigma_{ij},\sigma_{ik},\sigma_{im})` of this reference tetrahedron are the sums of the corresponding particle radii, :math:`\sigma_{ij} = r_i + r_j`, etc.

The irregularity of the tetrahedron :math:`\langle ijkm \rangle` in the original configuration is measured as

.. math::
	\omega_{\langle ijkm \rangle} = \frac{ \sum_{\langle ab \rangle} | r_{ab} - \sigma_{ab} |}{\sum_{\langle ab \rangle} \sigma_{ab}} ,

where :math:`\langle a b \rangle` runs over the six edges of the tetrahedron :math:`\langle ijkm \rangle`. Finally, the compactness of particle :math:`i` is given by

.. math::
	\Omega(i) = \frac{1}{N_\mathrm{tetra}(i)} \sum_{\langle ijkm \rangle} \omega_{\langle ijkm \rangle} ,

where :math:`N_\mathrm{tetra}(i)` is the total number of tetrahedra surrounding particle :math:`i` and the summation is performed over all these tetrahedra. The resulting feature vector for particle :math:`i` is given by

.. math::
	X^\mathrm{C}(i) = (\: \Omega(i) \:) .

.. note::
	Unlike most descriptors, this descriptor is **scalar**. Its feature vector :math:`X^\mathrm{C}(i)` is thus composed of a single feature, and the inherited ``grid`` attribute is therefore not relevant.

Setup
-----

.. warning::
	This descriptor uses particles' radii. They must either be read from the input trajectory file using the ``additional_fields`` parameter or set using the method :py:class:`Trajectory.set_property <partycls.trajectory.Trajectory.set_property>`. Otherwise, default values will be used.

Instantiating this descriptor on a :py:class:`Trajectory <partycls.trajectory.Trajectory>` can be done as follows:

.. code-block:: python

	from partycls import Trajectory
	from partycls.descriptors import CompactnessDescriptor

	traj = Trajectory("trajectory.xyz")
	D = CompactnessDescriptor(traj)

The constructor takes the following parameters:

.. automethod:: partycls.descriptors.compactness.CompactnessDescriptor.__init__

.. hint::

	The alias :py:class:`TongTanakaDescriptor <partycls.descriptors.compactness.TongTanakaDescriptor>` can be used in place of :py:class:`CompactnessDescriptor <partycls.descriptors.compactness.CompactnessDescriptor>`.

Requirements
------------

The computation of this descriptor relies on:

- **Lists of nearest neighbors**. These can either be read from the input trajectory file, computed in the :py:class:`Trajectory <partycls.trajectory.Trajectory>`, or computed from inside the descriptor using a default method.

Demonstration
-------------

We consider an input trajectory file :file:`trajectory.xyz` in XYZ format that contains two particle types ``"A"`` and ``"B"``. We set the particle radii manually (0.5 and 0.4 for type-A and type-B particles, respectively) and compute the lists of nearest neighbors using the radical Voronoi tessellation method:

.. code-block:: python

	from partycls import Trajectory

	# open the trajectory
	traj = Trajectory("trajectory.xyz")

	# set the particles radii
	traj.set_property("radius", 0.5, subset="species == 'A'")
	traj.set_property("radius", 0.4, subset="species == 'B'")

	# compute the neighbors using Voronoi tessellation
	traj.compute_nearest_neighbors(method='voronoi')
	nearest_neighbors = traj.get_property("nearest_neighbors")
	
	# print the first three neighbors lists for the first trajectory frame
	print("neighbors:\n",nearest_neighbors[0][0:3])

.. code-block:: none
	:caption: **Output:**

	neighbors:
	 [list([323, 113, 322, 276, 767, 332, 980, 425, 16, 171, 258, 801, 901, 436, 241])
	  list([448, 951, 706, 337, 481, 536, 14, 16, 258, 241, 496, 574, 502, 447, 860])
	  list([639, 397, 799, 578, 230, 913, 636, 796, 640, 772, 500, 270, 354, 123, 874, 608, 826, 810])]

We now instantiate a :py:class:`CompactnessDescriptor <partycls.descriptors.compactness.CompactnessDescriptor>` on this trajectory and restrict the analysis to type-B particles only:

.. code-block:: python

	from partycls.descriptors import CompactnessDescriptor

	# instantiation
	D = CompactnessDescriptor(traj)

	# restrict the analysis to type-B particles
	D.add_filter("species == 'B'", group=0)

	# compute the descriptor's data matrix
	X = D.compute()

	# print the first three feature vectors
	print("feature vectors:\n", X[0:3])

.. code-block:: none
	:caption: **Output:**

	feature vectors:
	 [[1.38606558]
	  [2.62615932]
	  [1.71241472]]

-  ``feature vectors`` shows the first three feature vectors :math:`X^\mathrm{C}(1)`, :math:`X^\mathrm{C}(2)` and :math:`X^\mathrm{C}(3)`.

References
----------

.. bibliography:: ../../references.bib
	:style: unsrt
	:filter: docname in docnames
