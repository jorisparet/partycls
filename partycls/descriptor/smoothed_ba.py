import numpy
from .descriptor import StructuralDescriptor, AngularStructuralDescriptor
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
        
    cutoff_enlargement : float
        Consider neighbors j `cutoff_enlargement * rc_ij` away from the central
        particle i.
        
    exponent : int
        Exponent `n` in the power law for the exponential decay in w(r).
        
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

    def __init__(self, trajectory, dtheta=3.0, cutoff_enlargement=1.3, exponent=8):
        BondAngleDescriptor.__init__(self, trajectory, dtheta=3.0)
        self.cutoff_enlargement = cutoff_enlargement
        self.exponent = exponent        

    def compute(self):
        # set up
        StructuralDescriptor._set_up(self, dtype=numpy.float64)
        AngularStructuralDescriptor._manage_nearest_neighbors(self)
        n_frames = len(self.trajectory)
        row = 0
        # all relevant arrays
        pos_0 = self.dump('position', group=0)
        pos_1 = self.dump('position', group=1)
        idx_0 = self.dump('index', group=0)
        spe_0_id = self.dump('species_id', group=0)
        spe_1_id = self.dump('species_id', group=1)
        box = self.trajectory.dump('cell.side')
        pairs = numpy.asarray(self.trajectory[0].pairs_of_species_id)
        # compute extended neighbors with extended cutoffs
        standard_cutoffs = numpy.asarray(self.trajectory.nearest_neighbors_cutoffs)
        extended_cutoffs = self.cutoff_enlargement * standard_cutoffs
        AngularStructuralDescriptor._compute_extended_neighbors(self, extended_cutoffs)
        # computation
        for n in range(n_frames):
            for i in range(len(idx_0[n])):
                hist_n_i = compute.smoothed_angular_histogram(idx_0[n][i],
                                                              pos_0[n][i], pos_1[n].T,
                                                              spe_0_id[n][i], spe_1_id[n],
                                                              self._extended_neighbors[n][i],
                                                              pairs, standard_cutoffs,
                                                              self.exponent, box[n],
                                                              self.n_features,
                                                              self.dtheta)
                self.features[row] = hist_n_i
                row += 1
        return self.features
