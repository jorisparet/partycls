import numpy
from .ba import BondAngleDescriptor
from .realspace_wrap import compute


class SmoothedBondAngleDescriptor(BondAngleDescriptor):
    """
    Smoothed bond-angle descriptor.
    
    Cutoffs `rc_ij` for the pair (i,j) of nearest neighbors are computed using 
    the corresponding partial RDF, g_ij(r), but more neighbors are considered 
    by looking further away from the central particle, using a 
    `cutoff_enlargement` parameter. The binning of the angle between the 
    central particle i and its neighbors (j,k) is then weighted by an 
    exponential decay w(r_ij, r_ik) that depends on the distances from i, 
    such that:
        
    w(r_ij, r_ik) = exp[ -( (r_ij/rc_ij)^n + (r_ik/rc_ik)^n ) ]
    
    See the parent class for more details.
    
    Parameters
    ----------
    
    trajectory : str or an instance of `Trajectory`.
        Trajectory on which the structural descriptor will be computed.
        
    dtheta : float
        Bin width in degrees.
        
    cutoff_enlargement : float, default: 1.3
        Consider neighbors j `cutoff_enlargement * rc_ij` away from the central
        particle i.
        
    exponent : int, default: 8
        Exponent `n` in the power law for the exponential decay in w(r).
        
    accept_nans: bool, default: True
        If False, discard any row from the array of features that contains a Nan
        element. If True, keep NaN elements in the array of features.

    verbose : bool, default: False
        Show progress information and warnings about the computation of the 
        descriptor when verbose is True, and remain silent when verbose is False.

    Attributes
    ----------
    
    trajectory : Trajectory
        Trajectory on which the structural descriptor will be computed.
        
    active_filters : list of str
        All the active filters on both groups prior to the computation of the
        descriptor.
        
    dimension : int
        Spatial dimension of the descriptor (2 or 3).
        
    grid : array
        Grid over which the structural features will be computed.
        
    features : ndarray
        Array of all the structural features for the particles in group=0 in
        accordance with the defined filters (if any). This attribute is 
        initialized when the method `compute` is called (default value is None).
        
    Examples:
    ---------
    
    >>> D = SmoothedBondAngleDescriptor('trajectory.xyz', cutoff_enlargement=1.3, exponent=8)
    >>> D.add_filter("species == 'A'", group=0)
    >>> D.compute()
    """

    name = 'smoothed-bond-angle'
    symbol = 'sba'

    def __init__(self, trajectory, dtheta=3.0, cutoff_enlargement=1.3, exponent=8,
                 accept_nans=True, verbose=False):
        BondAngleDescriptor.__init__(self, trajectory,
                                     dtheta=dtheta,
                                     accept_nans=accept_nans,
                                     verbose=verbose)
        self.cutoff_enlargement = cutoff_enlargement
        self.exponent = exponent        

    def compute(self):
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
