import numpy
from .descriptor import StructuralDescriptor
from .bo import BondOrientationalDescriptor
from .realspace_wrap import compute

class SmoothedBondOrientationalDescriptor(BondOrientationalDescriptor):
    """
    Smoothed bond orientational descriptor.
    
    Cutoffs `rc_ij` for the pair (i,j) of nearest neighbors are computed using 
    the corresponding partial RDF (or provided cutoffs), but more 
    neighbors are considered by looking further away from the central particle,
    using a `cutoff_enlargement` parameter. The value of q_lm between the 
    central particle i and one of its neighbors j is then weighted by an 
    exponential decay w(r_ij) that depends on the distance from i, such that:
        
    w(r_ij) = exp[ -(r_ij/rc_ij)^n) ]
    
    See the parent class for more details.
    
    Parameters
    ----------
    
    trajectory : str or an instance of `Trajectory`.
        Trajectory on which the structural descriptor will be computed.
        
    lmin : int, default: 0
        Minimum degree. This set the lower bound of the grid.
        
    lmax : int, default: 8
        Minimum degree. This set the upper bound of the grid.
        
    orders : list, default: None
        Specific values of orders to compute, e.g. orders=[4,6]. This has
        the priority over `lmin` and `lmax`.
        
    cutoff_enlargement : float, default: 1.3
        Consider neighbors j `cutoff_enlargement * rc_ij` away from the central
        particle i.
        
    exponent: int, default : 8
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
        
    nearest_neighbors_method : str, default: 'auto'
        Nearest neighbor method, 'FC' or 'SANN'. If method is 'auto', neighbors
        are read directly from the trajectory (if provided). If no neighbors 
        are found, it uses method='FC' instead.
    """

    name = 'smoothed bond-orientational'
    symbol = 'sbo'
    
    def __init__(self, trajectory, lmin=1, lmax=8, orders=None, cutoff_enlargement=1.3, exponent=8):
        BondOrientationalDescriptor.__init__(self, trajectory, lmin=lmin, lmax=lmax, orders=orders)
        self.cutoff_enlargement = cutoff_enlargement
        self.exponent = exponent       

    def compute(self):
        # set up
        StructuralDescriptor._set_up(self, dtype=numpy.float64)
        n_frames = len(self.trajectory)
        row = 0
        # all relevant arrays
        pos_0 = self.dump('position', group=0)
        pos_1 = self.dump('position', group=1)
        idx_0 = self.dump('index', group=0)
        spe_0 = self.dump('species_id', group=0)
        spe_1 = self.dump('species_id', group=1)
        pairs = numpy.asarray(self.trajectory[0].pairs_of_species_id)
        # Force 'FC' to compute neighbors
        if self.nearest_neighbors_method != 'FC':
            print("Warning: overriding the computation of nearest neighbors by using extended cutoffs.")
            self.nearest_neighbors_method = 'FC'
        
        # Recompute cutoffs if not provided
        if None in self.cutoffs:
            self._compute_cutoffs()
        # store the standard FC cutoffs
        self.standard_cutoffs_FC = numpy.array(self.cutoffs)
        
        # compute nearest neighbors with enlarged cutoffs
        self.cutoffs = list(self.cutoff_enlargement * numpy.array(self.cutoffs))
        self.nearest_neighbors(method=self.nearest_neighbors_method)
        
        for n in range(n_frames):
            box = self.trajectory[n].cell.side
            for i in range(len(idx_0[n])):
                # compute BO parameters for particle `i`
                neigh_i = self.neighbors[n][i]
                hist_n_i = numpy.empty_like(self.grid, dtype=numpy.float64)
                for ln, l in enumerate(self.grid):
                    hist_n_i[ln] = compute.smoothed_ql(l, neigh_i, pos_0[n][i], pos_1[n].T,
                                                       spe_0[n][i], spe_1[n], pairs,
                                                       box, self.standard_cutoffs_FC,
                                                       self.exponent)
                self.features[row] = hist_n_i
                row += 1
        return self.features