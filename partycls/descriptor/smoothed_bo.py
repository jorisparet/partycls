import numpy
from .descriptor import StructuralDescriptor, AngularStructuralDescriptor
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
        
    lmin : int, default: 1
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

    Examples:
    ---------
    
    >>> D = SmoothedBondOrientationalDescriptor('trajectory.xyz', cutoff_enlargement=1.3, exponent=8)
    >>> D.add_filter("species == 'A'", group=0)
    >>> D.compute()   
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
        AngularStructuralDescriptor._manage_nearest_neighbors(self)
        n_frames = len(self.trajectory)
        row = 0
        # all relevant arrays
        pos_0 = self.dump('position', group=0)
        pos_all = self.trajectory.dump('position')
        spe_0_id = self.dump('species_id', group=0)
        spe_all_id = self.trajectory.dump('species_id')
        box = self.trajectory.dump('cell.side')
        pairs = numpy.array(self.trajectory[0].pairs_of_species_id)
        # compute extended neighbors with extended cutoffs
        standard_cutoffs = numpy.array(self.trajectory.nearest_neighbors_cutoffs)
        extended_cutoffs = self.cutoff_enlargement * standard_cutoffs
        AngularStructuralDescriptor._compute_extended_neighbors(self, extended_cutoffs)
        # computation
        for n in range(n_frames):
            for i in range(len(self.groups[0][n])):
                hist_n_i = numpy.empty_like(self.grid, dtype=numpy.float64)
                for ln, l in enumerate(self.grid):
                    hist_n_i[ln] = compute.smoothed_ql(l, self._extended_neighbors[n][i], 
                                                       pos_0[n][i], pos_all[n].T,
                                                       spe_0_id[n][i], spe_all_id[n], pairs,
                                                       box[n], standard_cutoffs,
                                                       self.exponent)
                self.features[row] = hist_n_i
                row += 1
        return self.features