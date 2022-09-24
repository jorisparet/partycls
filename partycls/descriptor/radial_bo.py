import numpy
from .descriptor import StructuralDescriptor
from .bo import BondOrientationalDescriptor
from .realspace_wrap import compute

class RadialBondOrientationalDescriptor(BondOrientationalDescriptor):
    """
    Distance-dependent bond orientational descriptor.
    
    The bond orientational parameters q_l are computed over a grid of distances
    [r_min, ..., r_max] instead of considering all the neighbors of the central
    particle at once. The value of q_lm between the central particle i and one 
    of its neighbors j is thus weighted by
    
    w(r_ij, r, delta) = exp[-0.5 * ((r_ij-r)^n / delta^n)  ]
    
    where `r` is an arbitrary distance in the predefined grid of distances and
    `delta` is a shell width for the decay of the exponential around r. The 
    computation of neighbors is overriden to include all particles up to
    r_max + skin * delta, where `skin` and `delta` are tunable parameters.
    Any provided cutoff or nearest-neighbor method will thus be ignored.
    
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
        
    bounds : tuple of floats, default: (1, 2.5)
        Lower and upper bounds (r_min, r_max) to define the grid of distances.
        
    dr : float, default: 0.1
        Grid spacing to define the grid of distances [r_min, r_min+dr, ..., r_max].
        
    distance_grid : list or array, default: None
        Manually defined grid of distances. If different from `None`, it 
        overrides the linearly-spaced grid defined by `bounds` and `dr`.
        
    delta : float, default: 0.1
        Shell width to probe the local density at a distance r from the central
        particle using an exponential decay of the form:
        w(r_ij, r, \delta) = exp[-0.5 * (r_ij-r)^2 / (2 \delta^2) ]
        
    skin : float, default: 2.5
        Skin width (in units of `delta`) to consider neighbors further than the
        upper bound r_max of the grid of distances. Neighbors will then be 
        identified up to r_max + skin_width * delta.
        
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
    
    name = 'radial bond-orientational'
    symbol = 'rbo'
    
    
    def __init__(self, trajectory, lmin=1, lmax=8, orders=None, bounds=(1,2.5), dr=0.1, distance_grid=None, delta=0.1, skin=2.5, exponent=2):
        BondOrientationalDescriptor.__init__(self, trajectory, lmin=lmin, lmax=lmax, orders=orders)
        # Set the grid of distances
        self._set_bounds(dr, bounds, distance_grid)
        self.delta = delta
        self.skin = skin
        self.exponent = exponent
        
    @property
    def bounds(self):
        """
        Lower and upper bounds of the distance grid.
        """
        return self._bounds

    @bounds.setter
    def bounds(self, value):
        self._set_bounds(self._dr, value, self._distance_grid)

    @property
    def dr(self):
        """
        Grid spacing.
        """
        return self._dr

    @dr.setter
    def dr(self, value):
        self._set_bounds(value, self._bounds, self._distance_grid)
    
    @property
    def distance_grid(self):
        """
        Grid of distances over which to evaluate the descriptor.
        """
        return self._distance_grid
    
    @distance_grid.setter
    def distance_grid(self, value):
        self._set_bounds(self._dr, self._bounds, value)
    
    @property
    def n_features(self):
        """
        Number of features of the descriptor.
        """
        return len(self.grid) * len(self._distance_grid)
    
    def compute(self):
        # set up
        StructuralDescriptor._set_up(self, dtype=numpy.float64)
        n_frames = len(self.trajectory)
        row = 0
        # all relevant arrays
        pos_0 = self.dump('position', group=0)
        pos_1 = self.dump('position', group=1)
        idx_0 = self.dump('index', group=0)
        # compute the neighbors up to the largest distance in the grid + some skin distance
        R_cut = self.distance_grid[-1] + self.skin * self.delta
        if self.nearest_neighbors_method != 'FC':
            print("Warning: overriding the computation of nearest neighbors by using fixed cutoffs at a distance R={rcut:.6f}.".format(rcut=R_cut))
            self.nearest_neighbors_method = 'FC'
        else:
            # override any provided cutoff
            if any(self.cutoffs):
                print("Warning: overriding the provided cutoffs with R={rcut:.6f}.".format(rcut=R_cut))
        n_pairs = len(self.trajectory[0].pairs_of_species)
        self.cutoffs = numpy.array([R_cut for i in range(n_pairs)])
        self.nearest_neighbors(method=self.nearest_neighbors_method)
        
        for n in range(n_frames):
            box = self.trajectory[n].cell.side
            for i in range(len(idx_0[n])):
                # compute BO parameters for particle `i`
                neigh_i = self.neighbors[n][i]
                hist_n_i = numpy.empty_like(self.features[0], dtype=numpy.float64)
                feature_idx = 0
                for l in self.grid:
                    for r in self.distance_grid:
                        hist_n_i[feature_idx] = compute.radial_ql(l, r, 
                                                                  self.delta,
                                                                  self.exponent,
                                                                  neigh_i, 
                                                                  pos_0[n][i], 
                                                                  pos_1[n].T,
                                                                  box)
                        feature_idx += 1
                        
                self.features[row] = hist_n_i
                row += 1
        return self.features
        
    def _set_bounds(self, dr, bounds, distance_grid):
        # take the smallest side as maximal upper bound for the distance grid
        sides = numpy.array(self.trajectory.get_property('cell.side'))
        L = numpy.min(sides)
        
        if distance_grid is not None:
            self._distance_grid = numpy.asarray(distance_grid, dtype=numpy.float64)
            self._dr = None
            self._bounds = (self._distance_grid[0], self._distance_grid[-1])
        
        else:
            self._dr = dr
            if len(bounds) == 2 and bounds[0] >= 0 and bounds[1] >= 0 and bounds[0] < bounds[1] and bounds[1] <= L / 2:
                    rmin, rmax = bounds
                    r = numpy.arange(rmin + (self._dr / 2), rmax, self._dr, dtype=numpy.float64)
                    # set grid and bounds
                    self.distance_grid = r
                    self._bounds = (r[0], r[-1])
            else:
                raise ValueError('`bounds` is not correctly defined.')
