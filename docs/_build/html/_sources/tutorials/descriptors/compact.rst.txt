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
	Unlike most descriptors, this descriptor is scalar. Its feature vector :math:`X^\mathrm{C}(i)` is thus composed of a single feature.

Setup
-----

.. warning::
	This descriptor uses particles' radii. They must either be read from the input trajectory file using the ``additional_fields`` parameter or set using the method ``trajectory.set_property``. Otherwise, default values will be used.

Instantiating this descriptor on a ``Trajectory`` can be done as follows:

.. code-block:: python

	from partycls import Trajectory
	from partycls.descriptor import CompactnessDescriptor

	traj = Trajectory("trajectory.xyz")
	D = CompactnessDescriptor(traj)

The constructor takes the following parameters:

.. automethod:: partycls.descriptor.compactness.CompactnessDescriptor.__init__

Demonstration
-------------

References
----------

.. bibliography:: ../../references.bib
	:style: unsrt
	:filter: docname in docnames
