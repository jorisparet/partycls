import numpy
from .bo import BondOrientationalDescriptor
from .realspace_wrap import compute

class SmoothedBondOrientationalDescriptor(BondOrientationalDescriptor):
    """
    Smoothed bond-orientational descriptor.
    
    This is a smooth version of the bond-orientational descriptor, in which the 
    coefficients :math:`q_{lm}(i)` are multiplied by a weighting function 
    :math:`f(r)` that depends on the radial distance :math:`r` between the central 
    particle :math:`i` and other surrounding particles :math:`j`, where :math:`j` 
    can be any particle in the system (*i.e.* not necessarily a nearest neighbors 
    of :math:`i`).

    The smoothed complex coefficients are given by

    .. math::
        q_{lm}^{S}(i) = \\frac{1}{Z(i)} \\sum_{j=1}^{N} f({r}_{ij}) Y_{lm}(\hat{\mathbf{r}}_{ij}) ,

    where :math:`Z(i)=\\sum_{j=1}^{N} f({r}_{ij})` is a normalization constant and 
    the superscript :math:`S` indicates the smooth nature of the descriptor. We 
    use

    .. math::
        f(r_{ij}) = \exp \left[- (r_{ij} / r_{\\alpha\\beta}^c)^\gamma \\right] H(R_{\\alpha\\beta}^c - r_{ij}) ,

    where :math:`r_{\\alpha\\beta}^c` is the first minimum of the corresponding 
    partial radial distribution function for the pair :math:`(i,j)` and 
    :math:`\gamma` is an integer. Also, :math:`H` is the 
    `Heaviside step function <https://en.wikipedia.org/wiki/Heaviside_step_function>`_, 
    which ensures, for efficiency reasons, that the descriptor only has 
    contributions from particles within a distance 
    :math:`R_{\\alpha\\beta}^c = \\xi \\times r_{\\alpha\\beta}^c` from the central
    one, where :math:`\\xi > 1` is a scaling factor.

    The rotational invariants are defined similarly to the bond-orientational 
    descriptor.

    We then consider :math:`Q_l^S(i)` for a sequence of orders 
    :math:`\{ l_n \} = \{ l_\mathrm{min}, \dots, l_\mathrm{max} \}`. The resulting 
    feature vector for particle :math:`i` is given by

    .. math::
        X^\mathrm{SBO}(i) = (\: Q_{l_\mathrm{min}}^S(i) \;\; \dots \;\; Q_{l_\mathrm{max}}^S(i) \:) .
    
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
        Grid of orders :math:`\{ l_n \}`.
        
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

    name = 'smoothed bond-orientational'
    symbol = 'sbo'
    
    def __init__(self, trajectory, lmin=1, lmax=8, orders=None,
                cutoff_enlargement=1.3, exponent=8,
                accept_nans=True, verbose=False):
        """
        Parameters
        ----------
        trajectory : Trajectory
            Trajectory on which the structural descriptor will be computed.
            
        lmin : int, default: 1
            Minimum order :math:`l_\mathrm{min}`. This sets the lower bound of 
            the grid :math:`\{ l_n \}`.
            
        lmax : int, default: 8
            Maximum order :math:`l_\mathrm{max}`. This sets the upper bound of 
            the grid :math:`\{ l_n \}`. For numerical reasons, 
            :math:`l_\mathrm{max}` cannot be larger than 16.
            
        orders: list, default: None
            Sequence :math:`\{l_n\}` of specific orders to compute, *e.g.* 
            ``orders=[4,6]``. This has the priority over ``lmin`` and ``lmax``.
            
        cutoff_enlargement : float, default: 1.3
            Scaling factor :math:`\\xi` for the nearest neighbors cutoffs 
            :math:`r_{\\alpha\\beta}^c` to consider neighbors :math:`j` a distance
            :math:`R_{\\alpha\\beta}^c = \\xi \\times r_{\\alpha\\beta}^c` away from the
            central particle :math:`i`.
            
        exponent : int, default: 8
            Exponent :math:`\gamma` in the smoothing function
            :math:`f(r_{ij})`.

        accept_nans: bool, default: True
            If ``False``, discard any row from the array of features that contains a 
            `NaN` element. If ``True``, keep `NaN` elements in the array of features.

        verbose : bool, default: False
            Show progress information and warnings about the computation of the 
            descriptor when verbose is ``True``, and remain silent when verbose 
            is ``False``.
        """
        BondOrientationalDescriptor.__init__(self, trajectory,
                                             lmin=lmin, lmax=lmax,
                                             orders=orders,
                                             accept_nans=accept_nans,
                                             verbose=verbose)
        self.cutoff_enlargement = cutoff_enlargement
        self.exponent = exponent       

    def compute(self):
        """
        Compute the smoothed bond-orientational correlations for the particles in 
        group=0 for the grid of orders :math:`\{ l_n \}`. Returns the data matrix 
        and also overwrites the ``features`` attribute.

        Returns
        -------
        features : numpy.ndarray
            Data matrix with bond-orientational correlations.
        """
        # set up
        self._set_up(dtype=numpy.float64)
        self._manage_nearest_neighbors_cutoffs()
        n_frames = len(self.trajectory)
        # all relevant arrays
        pos_0 = self.dump('position', group=0)
        pos_all = self.trajectory.dump('position')
        spe_0_id = self.dump('species_id', group=0)
        spe_all_id = self.trajectory.dump('species_id')
        box = self.trajectory.dump('cell.side')
        n_species = len(self.trajectory[0].distinct_species)
        # compute extended neighbors with extended cutoffs
        standard_cutoffs = numpy.array(self.trajectory.nearest_neighbors_cutoffs)
        extended_cutoffs = self.cutoff_enlargement * standard_cutoffs
        self._compute_extended_neighbors(extended_cutoffs)
        # computation
        # TODO: it should not be necessary to transpose this
        standard_cutoffs = standard_cutoffs.reshape(n_species, n_species).T
        start = 0
        for n in self._trange(n_frames):
            # TODO: perhaps no transpose is better
            pos_0_n = pos_0[n].T
            pos_all_n = pos_all[n].T
            npart = len(self.groups[0][n])
            for ln, l in enumerate(self.grid):
                feat_n = compute.smoothed_ql_all(l,
                                                 self._extended_neighbors[n],
                                                 self._extended_neighbors_number[n],
                                                 pos_0_n, pos_all_n,
                                                 spe_0_id[n], spe_all_id[n],
                                                 box[n], standard_cutoffs,
                                                 self.exponent)
                self.features[start: start+npart, ln] = feat_n
            self._handle_nans()
            start += npart
        return self.features
