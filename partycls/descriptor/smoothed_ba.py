import numpy
from .descriptor import StructuralDescriptor
from .ba import BondAngleDescriptor
from .realspace_wrap import compute


class SmoothedBondAngleDescriptor(BondAngleDescriptor):
    """
    Smoothed bond angle descriptor.
    
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
        
    cutoffs : list of float
        List of enlarged cutoff distances to identify the nearest neighbors 
        using the fixed-cutoff ('FC') method.
        
    standard_cutoffs_FC : list of float
        List of standard cutoffs (i.e. not enlarged) with the fixed-cutoff 
        ('FC') method.
        
    nearest_neighbors_method : str, default: 'FC'
        Nearest neighbor method, 'FC' or 'SANN'.
        
    Examples:
    ---------
    
    >>> D = SmoothedBondAngleDescriptor('trajectory.xyz', cutoff_enlargement=1.3, exponent=8)
    >>> D.nearest_neighbors_method = 'SANN'
    >>> D.add_filter("species == 'A'")
    >>> D.compute()   
    """

    name = 'smoothed-bond-angle'
    symbol = 'sba'

    def __init__(self, trajectory, dtheta=3.0, cutoff_enlargement=1.3, exponent=8):
        BondAngleDescriptor.__init__(self, trajectory, dtheta=3.0)
        self.cutoff_enlargement = cutoff_enlargement
        self.exponent = exponent        

    def compute(self):
        StructuralDescriptor._sanity_checks(self)
        # all relevant arrays
        n_frames = len(self.groups[0])
        pos_0 = self.dump('position', 0)
        pos_1 = self.dump('position', 1)
        idx_0 = self.dump('index', 0)
        spe_0 = self.dump('species_id', 0)
        spe_1 = self.dump('species_id', 1)
        pairs = numpy.asarray(self.trajectory[0].pairs_of_species_id)
        features = numpy.empty((self.size, self.n_features), dtype=numpy.float64)
        row = 0
        
        # compute cutoffs for the descriptor if not provided
        if None in self.cutoffs:
            if self.nearest_neighbors_method != 'FC':
                print("Warning: using fixed cutoffs for the computation of \
                      the descriptor. Neighbors are determined using the \
                     `{}` method.".format(self.nearest_neighbors_method))
            self._compute_cutoffs()
        # store the standard FC cutoffs
        self.standard_cutoffs_FC = numpy.array(self.cutoffs)
        
        # compute nearest neighbors with enlarged cutoffs
        self.cutoffs = list(self.cutoff_enlargement * numpy.array(self.cutoffs))
        self.nearest_neighbors(method=self.nearest_neighbors_method)
       
        for n in range(n_frames):
            box = self.trajectory[n].cell.side
            for i in range(len(idx_0[n])):
                # individual bond-angle distribution using nearest-neighbors
                neigh_i = self.neighbors[n][i]
                hist_n_i = compute.smoothed_angular_histogram(idx_0[n][i],
                                                              pos_0[n][i], pos_1[n].T,
                                                              spe_0[n][i], spe_1[n],
                                                              neigh_i,
                                                              pairs, self.standard_cutoffs_FC,
                                                              self.exponent, box,
                                                              self.n_features,
                                                              self.dtheta)
                features[row] = hist_n_i
                row += 1
        self.features = features
        return features
