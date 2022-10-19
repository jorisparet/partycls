import numpy
from .ba import BondAngleDescriptor
from .realspace_wrap import compute


class SmoothedBondAngleDescriptor(BondAngleDescriptor):
    """
    Smoothed bond-angle descriptor.

    This is a smooth version of the bond-angle descriptor in which the bond angles 
    :math:`\\theta_{jik}` are multiplied by a weighting function 
    :math:`f(r_{ij}, r_{ik})` that depends on the radial distances :math:`r_{ij}` 
    and :math:`r_{ik}` between :math:`(i,j)` and :math:`(i,k)` respectively, where 
    :math:`j` and :math:`k` can be any particle in the system (*i.e.* not 
    necessarily a nearest neighbors of :math:`i`).

    Essentially, we define the *smoothed* number of bond angles around particle 
    :math:`i` as

    .. math::
        N_i^S(\\theta_n) = \\sum_{j=1}^N \\sum_{\\substack{k=1 \\ k \\neq j}}^N f(r_{ij}, r_{ik}) \delta(\\theta_n - \\theta_{jik}) ,

    where the superscript :math:`S` indicates the smooth nature of the descriptor,
    and :math:`f(r_{ij}, r_{ik})` is the following smoothing function:

    .. math::
        f(r_{ij}, r_{ik}) = \exp \left[ - \left[ ( r_{ij} / r_{\\alpha\\beta}^c )^\gamma + ( r_{ik} / r_{\\alpha\\beta'}^c )^\gamma \\right] \\right] H( R_\mathrm{max}^c - r_{ij} ) H( R_\mathrm{max}^c - r_{ik}) ,

    where:

    - :math:`r_{\\alpha\\beta}^c` and :math:`r_{\\alpha\\beta'}^c` are the first minima of the corresponding partial radial distribution functions for the pairs :math:`(i,j)` and :math:`(i,k)`.
    - :math:`\gamma` is an integer.
    - :math:`R_\mathrm{max}^c = \\xi \\times \max(\{ r_{\\alpha\\beta}^c \})` is the largest nearest neighbor cutoff rescaled by :math:`\\xi > 1`.
    - :math:`H` is the `Heavide step function <https://en.wikipedia.org/wiki/Heaviside_step_function>`_, which ensures for efficiency reasons, that the descriptor only has contributions from particles within a distance :math:`R_\mathrm{max}^c`.

    We then consider :math:`N_i^S(\\theta_n)` for a set of angles 
    :math:`\{ \\theta_n \}` that go from :math:`\\theta_0 = 0^\circ` to 
    :math:`\\theta_{n_\mathrm{max}}=180^\circ` by steps of :math:`\Delta \\theta`. 
    The resulting feature vector for particle :math:`i` is given by

    .. math::
        X^\mathrm{SBA}(i) = (\: N_i^S(\\theta_0) \;\; N_i^S(\\theta_1) \;\; \dots \;\; N_i^S(\\theta_{n_\mathrm{max}}) \:) .

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

    name = 'smoothed-bond-angle'
    symbol = 'sba'

    def __init__(self, trajectory, dtheta=3.0, cutoff_enlargement=1.3, exponent=8,
                 accept_nans=True, verbose=False):
        """
        Parameters
        ----------
        trajectory : Trajectory
            Trajectory on which the structural descriptor will be computed.
            
        dtheta : float
            Bin width :math:`\Delta \\theta` in degrees.
            
        cutoff_enlargement : float, default: 1.3
            Scaling factor :math:`\\xi` for the largest nearest neighbors cutoff 
            :math:`\max(\{ r_\mathrm{\\alpha\\beta}^c \})` to consider neighbors :math:`j`
            a distance :math:`R_\mathrm{max}^c = \\xi \\times \max(\{r_{\\alpha\\beta}^c\})`
            away from the central particle :math:`i`.
            
        exponent : int, default: 8
            Exponent :math:`\gamma` in the smoothing function
            :math:`f(r_{ij},r_{ik})`.
            
        accept_nans: bool, default: True
            If ``False``, discard any row from the array of features that contains a 
            `NaN` element. If ``True``, keep `NaN` elements in the array of features.

        verbose : bool, default: False
            Show progress information and warnings about the computation of the 
            descriptor when verbose is ``True``, and remain silent when verbose 
            is ``False``.
        """
        BondAngleDescriptor.__init__(self, trajectory,
                                     dtheta=dtheta,
                                     accept_nans=accept_nans,
                                     verbose=verbose)
        self.cutoff_enlargement = cutoff_enlargement
        self.exponent = exponent        

    def compute(self):
        """
        Compute the smoothed bond-angle correlations for the particles in group=0
        for the grid of angles :math:`\{ \\theta_n \}`. Returns the data matrix and 
        also overwrites the ``features`` attribute.

        Returns
        -------
        features : numpy.ndarray
            Data matrix with smoothed bond-angle correlations.
        """
        # set up
        self._set_up(dtype=numpy.float64)
        self._manage_nearest_neighbors_cutoffs()
        n_frames = len(self.trajectory)
        row = 0
        # all relevant arrays
        pos_0 = self.dump('position', group=0)
        pos_all = self.trajectory.dump('position')
        idx_0 = self.dump('_index', group=0)
        spe_0_id = self.dump('species_id', group=0)
        spe_all_id = self.trajectory.dump('species_id')
        box = self.trajectory.dump('cell.side')
        n_species = len(self.trajectory[0].distinct_species)
        # compute extended neighbors with extended cutoffs
        standard_cutoffs = numpy.asarray(self.trajectory.nearest_neighbors_cutoffs)
        extended_cutoffs = self.cutoff_enlargement * standard_cutoffs
        self._compute_extended_neighbors(extended_cutoffs)
        # computation
        standard_cutoffs = standard_cutoffs.reshape(n_species, n_species).T
        start = 0
        for n in self._trange(n_frames):
            pos_0_n = pos_0[n].T
            pos_all_n = pos_all[n].T
            npart = len(self.groups[0][n])
            feat_n = compute.smoothed_angular_histogram_all(idx_0[n],
                                                            pos_0_n, pos_all_n,
                                                            spe_0_id[n], spe_all_id[n],
                                                            self._extended_neighbors[n],
                                                            self._extended_neighbors_number[n],
                                                            standard_cutoffs,
                                                            self.exponent, box[n],
                                                            self.n_features,
                                                            self.dtheta)
            self.features[start: start+npart, :] = feat_n
            start += npart
        self._handle_nans()
        return self.features
