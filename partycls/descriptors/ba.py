import numpy
from .descriptor import StructuralDescriptor
from .realspace_wrap import compute


class BondAngleDescriptor(StructuralDescriptor):
    """
    Bond-angle descriptor.

    The *empirical* distribution of bond angles :math:`N_i(\\theta_n)` around a 
    central particle :math:`i` is obtained by counting, for all the possible pairs
    of it nearest neighbors :math:`(j,k)`, the number of bond angles 
    :math:`\\theta_{jik}` between :math:`\\theta_n = n \\times \Delta \\theta` and 
    :math:`\\theta_{n+1} = (n+1) \\times \Delta \\theta`, where :math:`\Delta \\theta` 
    has the interpration of a bin width in a histogram (see Ref. 
    :cite:`paret_2020`).

    We then consider :math:`N_i(\\theta_n)` for a set of angles 
    :math:`\{ \\theta_n \}` that go from :math:`\\theta_0 = 0^\circ` to 
    :math:`\\theta_{n_\mathrm{max}}=180^\circ` by steps of :math:`\Delta \\theta`. 
    The resulting feature vector for particle :math:`i` is given by

    .. math::
        X^\mathrm{BA}(i) = (\: N_i(\\theta_0) \;\; N_i(\\theta_1) \;\; \dots \;\; N_i(\\theta_{n_\mathrm{max}}) \:) .
    
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
        
    grid : numpy.ndarray
        Grid of angles :math:`\{ \\theta_n \}`.
        
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

    name = 'bond-angle'
    symbol = 'ba'

    def __init__(self, trajectory, dtheta=3.0, accept_nans=True, verbose=False):
        """
        Parameters
        ----------
        trajectory : Trajectory
            Trajectory on which the structural descriptor will be computed.
            
        dtheta : float
            Bin width :math:`\Delta \\theta` in degrees.
        
        accept_nans: bool, default: True
            If ``False``, discard any row from the array of features that contains a 
            `NaN` element. If ``True``, keep `NaN` elements in the array of features.

        verbose : bool, default: False
            Show progress information and warnings about the computation of the 
            descriptor when verbose is ``True``, and remain silent when verbose 
            is ``False``.
        """
        StructuralDescriptor.__init__(self, trajectory,
                                      verbose=verbose,
                                      accept_nans=accept_nans)
        self._dtheta = dtheta
        self.grid = numpy.arange(dtheta / 2.0, 180.0, dtheta, dtype=numpy.float64)

    @property
    def dtheta(self):
        """
        Bin width :math:`\Delta \\theta` for the grid of angles 
        :math:`\{ \\theta_n \}`.
        """
        return self._dtheta

    @dtheta.setter
    def dtheta(self, value):
        self._dtheta = value
        self.grid = numpy.arange(value / 2.0, 180.0, value, dtype=numpy.float64)

    def compute(self):
        """
        Compute the bond-angle correlations for the particles in group=0
        for the grid of angles :math:`\{ \\theta_n \}`. Returns the data matrix and 
        also overwrites the ``features`` attribute.

        Returns
        -------
        features : numpy.ndarray
            Data matrix with bond-angle correlations.
        """
        # set up
        self._set_up(dtype=numpy.int64)
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
            feat_n = compute.angular_histogram_all(idx_0[n],
                                                   pos_0_n, pos_all_n,
                                                   self._neighbors[n],
                                                   self._neighbors_number[n],
                                                   box[n],
                                                   self.n_features,
                                                   self.dtheta)
            self.features[start: start+npart, :] = feat_n
            start += npart
        self._handle_nans()
        return self.features

    def normalize(self, distribution, method="sin"):
        """
        Normalize a bond-angle distribution.

        Parameters
        ----------
        distribution : numpy.ndarray
            Distribution to normalize.
            
        method : str, default: "sin"
            Normalization method:

            - ``method='sin'``: by construction, the probability density of has a sinusoidal enveloppe in **3D** for uniformly distributed points on a sphere (default)
            - ``method='pdf'`` : gives a flat probability density for uniformly distributed points on a sphere ;

        Raises
        ------
        ValueError
            If ``method`` is invalid.

        Returns
        -------
        numpy.ndarray
            Normalized distribution.
        """
        if method == "sin":
            return distribution / numpy.sum(distribution) / self.dtheta
        elif method == "pdf":
            norm_sin = numpy.sin(self.grid / 360.0 * 2 * numpy.pi)
            distribution = distribution / norm_sin
            return distribution / numpy.sum(distribution) / self.dtheta
        else:
            raise ValueError("unknown value {}".format(method))
