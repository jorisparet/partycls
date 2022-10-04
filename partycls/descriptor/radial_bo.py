import numpy
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

    verbose : int, default: 0
        Show progress information about the computation of the descriptor
        when verbose is 1, and remain silent when verbose is 0 (default).
        
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
        Grid of bond orientational orders `l` and distances `r` over which 
        the descriptor will be computed. It takes the form a of a list of
        tuples (l,r).
        
    features : ndarray
        Array of all the structural features for the particles in group=0 in
        accordance with the defined filters (if any). This attribute is 
        initialized when the method `compute` is called (default value is None).

    Examples:
    ---------
    
    >>> D = RadialBondOrientationalDescriptor('trajectory.xyz', distance_grid=[1.0, 1.1, 1.2], delta=0.2)
    >>> D.add_filter("species == 'A'", group=0)
    >>> D.compute()
    """
    
    name = 'radial bond-orientational'
    symbol = 'rbo'
    
    
    def __init__(self, trajectory, lmin=1, lmax=8, orders=None, bounds=(1,2.5), dr=0.1, distance_grid=None, delta=0.1, skin=2.5, exponent=2, verbose=0):
        BondOrientationalDescriptor.__init__(self, trajectory, lmin=lmin, lmax=lmax, orders=orders, verbose=verbose)
        # dummy values, to be set in the following lines
        self._distance_grid = [-1]
        self._orders = [-1]
        # set the grid of distances {r}
        self._set_bounds_distances(dr, bounds, distance_grid)
        # set the grid of orders {l} (overwrites method from BondOrientationalDescriptor)
        self._set_bounds_orders(lmin, lmax, orders)
        self.delta = delta
        self.skin = skin
        self.exponent = exponent
        
    @property
    def orders(self):
        """
        Grid of bond orders `l` used to compute the descriptor. 
        """
        return self._orders

    @orders.setter
    def orders(self, value):
        return self._set_bounds_orders(1, 8, value)

    @property
    def bounds(self):
        """
        Lower and upper bounds of the distance grid.
        """
        return self._bounds

    @bounds.setter
    def bounds(self, value):
        self._set_bounds_distances(self._dr, value, None)

    @property
    def dr(self):
        """
        Grid spacing for the distance grid.
        """
        return self._dr

    @dr.setter
    def dr(self, value):
        self._set_bounds_distances(value, self._bounds, None)
    
    @property
    def distance_grid(self):
        """
        Grid of distances used to compute the descriptor.
        """
        return self._distance_grid
    
    @distance_grid.setter
    def distance_grid(self, value):
        self._set_bounds_distances(self._dr, self._bounds, value)
    
    @property
    def mixed_grid(self):
        """
        Mixed grid of bond orientational orders `l` and
        distances `r` in the form of a list of tuples (l,r).
        """
        mixed_grid = []
        for l in self.grid:
            for r in self._distance_grid:
                mixed_grid.append((l,r))
        return mixed_grid

    def compute(self):
        # set up
        self._set_up(dtype=numpy.float64)
        self._manage_nearest_neighbors()
        n_frames = len(self.trajectory)
        row = 0
        # all relevant arrays
        pos_0 = self.dump('position', group=0)
        pos_all = self.trajectory.dump('position')
        box = self.trajectory.dump('cell.side')
        # compute extended neighbors with extended cutoffs
        # based on the largest distance in the distance grid
        R_cut = self.distance_grid[-1] + self.skin * self.delta
        n_pairs = len(self.trajectory[0].pairs_of_species)
        extended_cutoffs = numpy.array([R_cut for i in range(n_pairs)])
        self._compute_extended_neighbors(extended_cutoffs)
        # computation        
        for n in self._trange(n_frames):
            pos_all_n = pos_all[n].T
            for i in range(len(self.groups[0][n])):
                hist_n_i = numpy.empty_like(self.features[0], dtype=numpy.float64)
                feature_idx = 0
                for l in self._orders:
                    for r in self._distance_grid:
                        hist_n_i[feature_idx] = compute.radial_ql(l, r, 
                                                                  self.delta,
                                                                  self.exponent,
                                                                  self._extended_neighbors[n][i], 
                                                                  pos_0[n][i], 
                                                                  pos_all_n,
                                                                  box[n])
                        feature_idx += 1
                self.features[row] = hist_n_i
                row += 1
        return self.features
        
    def _set_bounds_distances(self, dr, bounds, distance_grid):
        # take the smallest side as maximal upper bound for the distance grid
        sides = numpy.array(self.trajectory.get_property('cell.side'))
        L = numpy.min(sides)
        
        if distance_grid is not None:
            self._distance_grid = numpy.array(distance_grid, dtype=numpy.float64)
            self._dr = None
            self._bounds = (self._distance_grid[0], self._distance_grid[-1])
        
        else:
            self._dr = dr
            if len(bounds) == 2 and bounds[0] >= 0 and bounds[1] >= 0 and bounds[0] < bounds[1] and bounds[1] <= L / 2:
                    rmin, rmax = bounds
                    r = numpy.arange(rmin + (self._dr / 2), rmax, self._dr, dtype=numpy.float64)
                    # set grid and bounds
                    self._distance_grid = r
                    self._bounds = (r[0], r[-1])
            else:
                raise ValueError('`bounds` is not correctly defined.')

        # update (l,r) grid
        self._set_grid()

    def _set_bounds_orders(self, lmin, lmax, orders):
        if orders is None:
            self._orders = numpy.array(range(lmin, lmax + 1))
        else:
            self._orders = numpy.array(orders)
        # update (l,r) grid
        self._set_grid()

    def _set_grid(self):
        self.grid = [(l,r) for l in self._orders for r in self._distance_grid]