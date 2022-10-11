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

    Attributes
    ----------
    trajectory : Trajectory
        Trajectory on which the structural descriptor will be computed.
        
    active_filters : list
        All the active filters on both groups prior to the computation of the
        descriptor.
        
    dimension : int
        Spatial dimension of the descriptor (2 or 3).
        
    grid : list
        Grid of bond orientational orders :math:`\{l_m\}` and distances 
        :math:`\{d_n\}` over which the descriptor will be computed. 
        It takes the form a of a list of tuples 
        ``[(l0,d0), (l0,d1), ..., (ln,dn), ...]``.
        
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
    """
    
    name = 'radial bond-orientational'
    symbol = 'rbo'
    
    
    def __init__(self, trajectory, lmin=1, lmax=8, orders=None,
                 bounds=(1,2.5), dr=0.1, distance_grid=None,
                 delta=0.1, skin=2.5, exponent=2,
                 accept_nans=True, verbose=False):
        """
        Parameters
        ----------
        trajectory : Trajectory
            Trajectory on which the structural descriptor will be computed.
            
        lmin : int, default: 1
            Minimum order :math:`l_\mathrm{min}`. This sets the lower bound of 
            the grid :math:`\{ l_m \}`.
            
        lmax : int, default: 8
            Maximum order :math:`l_\mathrm{max}`. This sets the upper bound of 
            the grid :math:`\{ l_m \}`.
            
        orders: list, default: None
            Sequence :math:`\{l_m\}` of specific orders to compute, *e.g.* 
            ``orders=[4,6]``. This has the priority over ``lmin`` and ``lmax``.
            
        bounds : tuple, default: (1, 2.5)
            Lower and upper bounds :math:`(d_{n_\mathrm{min}}, d_{n_\mathrm{max}})`
            to define the grid of distances :math:`\{ d_n \}`.
            
        dr : float, default: 0.1
            Grid spacing :math:`\Delta r` to define the grid of distances 
            :math:`\{ d_n \}` in the range 
            :math:`(d_{n_\mathrm{min}}, d_{n_\mathrm{max}})` by steps of size 
            :math:`\Delta r`.
            
        distance_grid : list, default: None
            Manually defined grid of distances :math:`\{ d_n \}`. If different 
            from ``None``, it  overrides the linearly-spaced grid defined by 
            ``bounds`` and ``dr``.
            
        delta : float, default: 0.1
            Shell width :math:`\delta` to probe the local density at a distances 
            :math:`\{d_n\}` from the central particle.
            
        skin : float, default: 2.5
            Skin width :math:`\sigma` (in units of ``delta``) to consider neighbors
            further than the upper bound :math:`d_{n_\mathrm{max}}` of the grid of 
            distances. Neighbors will then be  identified up to
            :math:`d_{n_\mathrm{max}} + \sigma \\times \delta`

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
        Grid of orders :math:`\{ l_m \}`.
        """
        return self._orders

    @orders.setter
    def orders(self, value):
        return self._set_bounds_orders(1, 8, value)

    @property
    def bounds(self):
        """
        Lower and upper bounds of the distance grid :math:`\{ d_n \}`.
        """
        return self._bounds

    @bounds.setter
    def bounds(self, value):
        self._set_bounds_distances(self._dr, value, None)

    @property
    def dr(self):
        """
        Grid spacing :math:`\Delta r` to define the distance grid :math:`\{ d_n \}`.
        """
        return self._dr

    @dr.setter
    def dr(self, value):
        self._set_bounds_distances(value, self._bounds, None)
    
    @property
    def distance_grid(self):
        """
        Grid of distances :math:`\{ d_n \}`.
        """
        return self._distance_grid
    
    @distance_grid.setter
    def distance_grid(self, value):
        self._set_bounds_distances(self._dr, self._bounds, value)
    
    @property
    def mixed_grid(self):
        """
        Mixed grids of bond orientational orders :math:`\{ l_m \}` and
        distances :math:`\{ d_n \}` in the form of a list of tuples ``(l,r)``.
        """
        mixed_grid = []
        for l in self.grid:
            for r in self._distance_grid:
                mixed_grid.append((l,r))
        return mixed_grid

    def compute(self):
        """
        Compute the radial bond-orientational correlations for the particles in group=0
        for the grid of orders :math:`\{ l_m \}` and grid of distances 
        :math:`\{ d_n \}`. Returns the data matrix and also
        overwrites the ``features`` attribute.

        Returns
        -------
        features : numpy.ndarray
            Data matrix with radial bond-orientational correlations.
        """
        # set up
        self._set_up(dtype=numpy.float64)
        self._manage_nearest_neighbors()
        n_frames = len(self.trajectory)
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
        start = 0
        for n in self._trange(n_frames):
            pos_0_n = pos_0[n].T
            pos_all_n = pos_all[n].T
            npart = len(self.groups[0][n])
            lr_idx = 0
            for l in self._orders:
                for r in self._distance_grid:
                    feat_lr_n = compute.radial_ql_all(l, r, 
                                                   self.delta,
                                                   self.exponent,
                                                   self._extended_neighbors[n],
                                                   self._extended_neighbors_number[n],
                                                   pos_0_n,
                                                   pos_all_n,
                                                   box[n])
                    self.features[start: start+npart, lr_idx] = feat_lr_n
                    lr_idx += 1
            start += npart
        self._handle_nans()
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