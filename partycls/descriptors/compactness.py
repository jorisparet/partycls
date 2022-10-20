import numpy
from .descriptor import StructuralDescriptor
from .realspace_wrap import compute

class CompactnessDescriptor(StructuralDescriptor):
    """
    Compactness descriptor.
    
    The compactness descriptor quantifies the local packing efficiency of a 
    particle :math:`i` by comparing it to a reference ideal packing configuration 
    of its nearest neighbors :cite:`tong_2018`.

    We consider a central particle :math:`i` surrounded by :math:`N_b(i)` nearest 
    neighbors, usually identified by means of a radical Voronoi tessellation 
    :cite:`gellatly_1982`. We then consider triplets :math:`(j,k,m)` of 
    neighboring particles, for which all particles are simultaneously nearest 
    neighbors of each other and of particle :math:`i`. Such a triplet is 
    identified with the central particle :math:`i` to be a tetrahedron.

    Particles in the tetrahedron :math:`\langle ijkm \\rangle` have radii 
    :math:`(r_i, r_j, r_k, r_m)` respectively, and the triplet :math:`(j,k,m)` 
    are at distances :math:`(r_{ij},r_{ik},r_{ik})` from the central particle 
    :math:`i`, which are the lengths of each edge of the tetrahedron 
    :math:`\langle ijkm \\rangle`.

    The *reference* tetrahedron for these four particles is the configuration in 
    which they are all perfectly in touch, *i.e.* the edge lengths 
    :math:`(\\sigma_{ij},\\sigma_{ik},\\sigma_{im})` of this reference tetrahedron 
    are the sums of the corresponding particle radii, 
    :math:`\\sigma_{ij} = r_i + r_j`, etc.

    The irregularity of the tetrahedron :math:`\langle ijkm \\rangle` in the 
    original configuration is measured as

    .. math::
        \omega_{\langle ijkm \\rangle} = \\frac{ \\sum_{\langle ab \\rangle} | r_{ab} - \\sigma_{ab} |}{\\sum_{\langle ab \\rangle} \\sigma_{ab}} ,

    where :math:`\langle a b \\rangle` runs over the six edges of the tetrahedron 
    :math:`\langle ijkm \\rangle`. Finally, the compactness of particle :math:`i` 
    is given by

    .. math::
        \Omega(i) = \\frac{1}{N_\mathrm{tetra}(i)} \\sum_{\langle ijkm \\rangle} \omega_{\langle ijkm \\rangle} ,

    where :math:`N_\mathrm{tetra}(i)` is the total number of tetrahedra 
    surrounding particle :math:`i` and the summation is performed over all these 
    tetrahedra. The resulting feature vector for particle :math:`i` is given by

    .. math::
        X^\mathrm{C}(i) = (\: \Omega(i) \:) .

    .. note::
        Unlike most descriptors, this descriptor is **scalar**. Its feature vector
        :math:`X^\mathrm{C}(i)` is thus composed of a single feature, and the inherited
        ``grid`` attribute is therefore not relevant.
    
    See the tutorials for more details.

    Attributes
    ----------
    trajectory : Trajectory
        Trajectory on which the structural descriptor will be computed.
        
    active_filters : list
        All the active filters on both groups prior to the computation of the
        descriptor.
        
    dimension : int
        Spatial dimension of the descriptor (2 or 3).
        
    features : numpy.ndarray
        Array of all the structural features for the particles in group=0 in
        accordance with the defined filters (if any). This attribute is 
        initialized when the method ``compute`` is called (default value is ``None``).

    groups : tuple
        Composition of the groups: ``groups[0]`` and ``groups[1]`` contain lists of all
        the ``Particle`` instances in groups 0 and 1 respectively. Each element of 
        the tuple is a list of ``Particle`` in ``trajectory``, *e.g.* ``groups[0][0]``
        is the list of all the particles in the first frame of ``trajectory`` that 
        belong to group=0.

    verbose : bool
        Show progress information and warnings about the computation of the 
        descriptor when verbose is ``True``, and remain silent when verbose is 
        ``False``.

    neighbors_boost : float, default: 1.5
        Scaling factor to estimate the number of neighbors relative to a
        an ideal gas with the same density. This is used internally to set
        the dimensions of lists of neighbors. A too small number creates a
        risk of overfilling the lists of neighbors, and a too large number
        increases memory usage. This only works if the associated ``Trajectory``
        has valid cutoffs in the ``Trajectory.nearest_neighbors_cutoffs`` list
        attribute. This sets the value of the ``max_num_neighbors`` attribute
        during the computation of the descriptor.

    max_num_neighbors : int, default: 100
        Maximum number of neighbors. This is used internally to set the dimensions
        of lists of neighbors. This number is automatically adjusted to limit
        memory usage if the associated ``Trajectory`` has valid cutoffs in the 
        ``Trajectory.nearest_neighbors_cutoffs`` list attribute. The
        default value ``100`` is used if no cutoffs can be used to estimate a
        better value. The default value is sufficient in most cases, otherwise 
        this number can manually be increased **before** computing the descriptor.
    """

    name = 'compactness'
    symbol = 'compact'
    
    def __init__(self, trajectory, accept_nans=True, verbose=False):
        """
        Parameters
        ----------
        trajectory : Trajectory
            Trajectory on which the structural descriptor will be computed.

        accept_nans: bool, default: True
            If ``False``, discard any row from the array of features that contains a 
            `NaN` element. If ``True``, keep `NaN` elements in the array of features.

        verbose : bool, default: False
            Show progress information and warnings about the computation of the 
            descriptor when verbose is ``True``, and remain silent when verbose 
            is ``False``.
        """
        StructuralDescriptor.__init__(self,
                                      trajectory,
                                      accept_nans=accept_nans,
                                      verbose=verbose)
        self._dimension_check(dimension=3)
        self.grid = numpy.zeros(1, dtype=numpy.float64)
        
    def compute(self):
        """
        Compute the compactness for the particles in group=0.
        Returns the data matrix and also overwrites the ``features`` attribute.

        Returns
        -------
        features : numpy.ndarray
            Data matrix with compactness.
        """
        # set up
        self._set_up(dtype=numpy.float64)
        self._manage_nearest_neighbors()
        self._filter_neighbors()
        self._filter_subsidiary_neighbors()
        n_frames = len(self.trajectory)
        row = 0
        # all relevant arrays
        pos_all = self.trajectory.dump('position')
        radii = self.trajectory.dump('radius')
        box = self.trajectory.dump('cell.side')
        # computation
        for n in self._trange(n_frames):
            pos_all_n = pos_all[n].T
            for i in range(len(self.groups[0][n])):
                idx_i = self.groups[0][n][i]._index
                nn_i = self._neighbors_number[n][i]
                tetra_i = self.tetrahedra(idx_i, 
                                          list(self._neighbors[n][i][0:nn_i]),
                                          self._subsidiary_neighbors[n][i])
                omega_i = compute.compactness(pos_all_n, tetra_i.T, radii[n], box[n])
                self.features[row] = omega_i
                row += 1
        self._handle_nans()
        return self.features
    
    def tetrahedra(self, i, neigh_i, neigh_neigh_i):
        """
        Return a 2D numpy array that contains all the tetrahedra around
        particle ``i``. Each row is a 4-array with the indices of the
        particles forming the tetrahedron. The first index is ``i``
        by construction.

        Parameters
        ----------
        i : int
            Index of the central particle.

        neigh_i : list
            List of neighbors of the central particles.

        neigh_neigh_i: list
            Neighbors of the neighbors of the central particle.
            Each element is a list.

        Returns
        -------
        tetra_i : numpy.ndarray
            The array that contains the indices of all the tetrahedra
            around particle ``i``.
        """
        tetras = []
        for j_idx, j in enumerate(neigh_i):
            neigh_j = neigh_neigh_i[j_idx]
            common = set(neigh_i) & set(neigh_j)
            for k in common:
                for m in common:
                    if m <= k:
                        continue
                    k_idx = neigh_i.index(k)
                    if m in neigh_neigh_i[k_idx]:
                        tetras.append([i,j,k,m])
        return numpy.array(tetras, dtype=numpy.int64)

class TongTanakaDescriptor(CompactnessDescriptor):
    """
    Alias for the class ``CompactnessDescriptor``.
    """
    pass