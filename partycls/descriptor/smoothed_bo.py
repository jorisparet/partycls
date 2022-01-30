import numpy
from .descriptor import StructuralDescriptor
from .bo import BondOrientationalDescriptor
from .realspace_wrap import compute

class SmoothedBondOrientationalDescriptor(BondOrientationalDescriptor):
    """
    Smoothed bond orientational descriptor.
    
    Cutoffs `rc_ij` for the pair (i,j) of nearest neighbors are computed using 
    the corresponding partial RDF, g_ij(r), but more neighbors are considered 
    by looking further away from the central particle, using a 
    `cutoff_enlargement` parameter. The value of q_lm between the central 
    particle i and one of its neighbors j is then weighted by an exponential
    decay w(r_ij) that depends on the distance from i, such that:
        
    w(r_ij, r_ik) = exp[ -(r_ij/rc_ij)^n) ]
    
    See the parent class for more details.
    
    Parameters
    ----------
    
    trajectory : str or an instance of `Trajectory`.
        Trajectory on which the structural descriptor will be computed.
        
    lmin : int, default: 0
        Minimum degree. This set the lower bound of the grid.
        
    lmax : int, default: 8
        Minimum degree. This set the upper bound of the grid.
        
    orders: list, default: None
        Specific values of orders to compute, e.g. orders=[4,6]. This has
        the priority over `lmin` and `lmax`.
        
    cutoff_enlargement : float
        Consider neighbors j `cutoff_enlargement * rc_ij` away from the central
        particle i.
        
    power_law : int
        Value `n` for the power law decay w(r).
    """

    name = 'smoothed bond-orientational'
    symbol = 'sbo'
    
    def __init__(self, trajectory, lmin=1, lmax=8, orders=None, cutoff_enlargement=1.3, power_law=8):
        BondOrientationalDescriptor.__init__(self, trajectory, lmin=lmin, lmax=lmax)
        self.cutoff_enlargement = cutoff_enlargement
        self.power_law = power_law       

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
        # override the computation of cutoffs from nearest_neighbors()
        self._compute_cutoffs()
        cutoffs = numpy.array(self.cutoffs)
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
                                                       box, cutoffs,
                                                       self.power_law)
                features[row] = hist_n_i
                row += 1
        self.features = features
        return features                