import numpy
from .descriptor import StructuralDescriptor
from .bo import BondOrientationalDescriptor
from .realspace_wrap import compute

class RadialBondOrientationalDescriptor(BondOrientationalDescriptor):
    """
    Distance-dependent bond orientational descriptor.
    
    Parameters
    ----------
    
    ...
    
    """
    
    name = 'radial bond-orientational'
    symbol = 'rbo'
    
    
    def __init__(self, trajectory, lmin=1, lmax=8, orders=None, bounds=(1,2.5), dr=0.1, delta=0.1):
        BondOrientationalDescriptor.__init__(self, trajectory, lmin=lmin, lmax=lmax)
        # Set the grid of distances
        self._set_bounds(dr, bounds)
        self.delta = delta
        
    @property
    def bounds(self):
        """
        Lower and upper bounds of the distance grid.
        """
        return self._bounds

    @bounds.setter
    def bounds(self, value):
        self._set_bounds(self._dr, value)

    @property
    def dr(self):
        """
        Grid spacing.
        """
        return self._dr

    @dr.setter
    def dr(self, value):
        self._set_bounds(value, self._bounds)
    
    @property
    def n_features(self):
        """
        Number of features of the descriptor.
        """
        return len(self.grid) * len(self.distance_grid)
    
    def compute(self):
        StructuralDescriptor._sanity_checks(self)
        # all relevant arrays
        n_frames = len(self.groups[0])
        pos_0 = self.dump('position', 0)
        pos_1 = self.dump('position', 1)
        idx_0 = self.dump('index', 0)
        features = numpy.empty((self.size, self.n_features), dtype=numpy.float64)
        row = 0
        # compute the neighbors up to the largest distance in the grid
        R_cut = self.distance_grid[-1]
        print("Warning: overriding the computation of nearest neighbors by using fixed cutoffs at a distance R={rcut:.6f}".format(rcut=R_cut))
        n_pairs = len(self.trajectory[0].pairs_of_species)
        self.cutoffs = numpy.array([R_cut for i in range(n_pairs)])
        self.nearest_neighbors(method='FC')
        for n in range(n_frames):
            box = self.trajectory[n].cell.side
            for i in range(len(idx_0[n])):
                # compute BO parameters for particle `i`
                neigh_i = self.neighbors[n][i]
                hist_n_i = numpy.empty_like(features[0], dtype=numpy.float64)
                for ln, l in enumerate(self.grid):
                    for rn, r in enumerate(self.distance_grid):
                        hist_n_i[ln+rn] = compute.radial_ql(l, r, self.delta, 
                                                         neigh_i, 
                                                         pos_0[n][i], 
                                                         pos_1[n].T,
                                                         box)
                features[row] = hist_n_i
                row += 1
        self.features = features
        return features
        
    def _set_bounds(self, dr, bounds):
        # take the smallest side as maximal upper bound for the distance grid
        sides = numpy.array(self.trajectory.get_property('cell.side'))
        L = numpy.min(sides)
        
        self._dr = dr
        if len(bounds) == 2 and bounds[0] < bounds[1] and bounds[1] <= L / 2:
                rmin, rmax = bounds
                r = numpy.arange(rmin + (self._dr / 2), rmax, self._dr, dtype=numpy.float64)
                # set grid and bounds
                self.distance_grid = r
                self._bounds = (r[0], r[-1])
        else:
            raise ValueError('`bounds` is not correctly defined.')
