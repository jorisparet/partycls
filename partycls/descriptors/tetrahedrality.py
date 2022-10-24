import numpy
from .descriptor import StructuralDescriptor
from .realspace_wrap import compute

class TetrahedralDescriptor(StructuralDescriptor):
    """
    Tetrahedral descriptor.
    
    The degree of tetrahedrality of a particle :math:`i` is the average deviation 
    of the bond angles :math:`\{ \\theta_{jik} \}` between :math:`i` and all the 
    possible pairs of its nearest neighbors :math:`(j,k)` from the ideal angle in
    a tetrahedron, :math:`\\theta_\mathrm{tetra} = 109.5^\circ`:

    .. math::
        T(i) = \\frac{1}{N_\mathrm{ba}(i)} \\sum_{j=1}^{N_b(i)} \\sum_{\\substack{k=1 \\ k \\neq j}}^{N_b(i)} | \cos(\\theta_{jik}) - \cos(\\theta_\mathrm{tetra}) | ,

    where :math:`N_\mathrm{ba}(i)` is the total number of bond angles (*i.e.* the
    number of pairs) around particle :math:`i` and :math:`N_b(i)` is the number 
    of its nearest neighbors. The resulting feature vector for particle 
    :math:`i` is given by

    .. math::
        X^\mathrm{T}(i) = (\: T(i) \:) .

    .. note::
        Unlike most descriptors, this descriptor is **scalar**. Its feature vector
        :math:`X^\mathrm{T}(i)` is thus composed of a single feature, and the 
        inherited ``grid`` attribute is therefore not relevant.
    
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
    
    name = 'tetrahedral'
    symbol = 'tetra'
    
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
        StructuralDescriptor.__init__(self, trajectory,
                                      accept_nans=accept_nans,
                                      verbose=verbose)
        self.grid = numpy.zeros(1, dtype=numpy.float64)
        
    def compute(self):
        """
        Compute the tetrahedrality for the particles in group=0.
        Returns the data matrix and also overwrites the ``features`` attribute.

        Returns
        -------
        features : numpy.ndarray
            Data matrix with the degree of tetrahedrality.
        """
        self._set_up(dtype=numpy.float64)
        self._manage_nearest_neighbors()
        self._filter_neighbors()
        n_frames = len(self.trajectory)
        # all relevant arrays
        pos_0 = self.dump('position', group=0)
        pos_all = self.trajectory.dump('position')
        idx_0 = self.dump('_index', group=0)
        box = self.trajectory.dump('cell.side')
        # computation
        start = 0
        for n in self._trange(n_frames):
            pos_0_n = pos_0[n].T
            pos_all_n = pos_all[n].T
            npart = len(self.groups[0][n])
            feat_n = compute.tetrahedrality_all(idx_0[n],
                                                pos_0_n, pos_all_n,
                                                self._neighbors[n],
                                                self._neighbors_number[n],
                                                box[n])
            self.features[start: start+npart, 0] = feat_n
            start += npart
        self._handle_nans()
        return self.features