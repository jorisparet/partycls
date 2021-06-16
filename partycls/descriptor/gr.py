import numpy
from .descriptor import StructuralDescriptor
from .realspace_wrap import compute

class RadialDescriptor(StructuralDescriptor):
    """
    Structural descriptor based on radial correlations between particles.
    
    See the parent class for more details.
    
    Parameters
    ----------
    
    trajectory : str or an instance of `Trajectory`.
        Trajectory on which the structural descriptor will be computed.
        
    dr : float
        Bin width.
        
    n_shells : int, default: 3
        Number of coordination shells (based on the RDF of group=0). This sets
        the upper bound for the distance up to which correlations are computed.
        
    bounds : tuple, default: None
        Lower and upper bounds to describe the radial correlations. If set, 
        this has the priority over `n_shells`.
    
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
    
    >>> D = RadialDescriptor('trajectory.xyz', bounds=(0.0,3.0))
    >>> D.add_filter("species == 'A'")
    >>> D.compute()
    """

    name = 'radial'
    symbol = 'gr'
    
    def __init__(self, trajectory, dr=0.1, n_shells=3, bounds=None):
        StructuralDescriptor.__init__(self, trajectory)
        # set the grid automatically using coordination shells (`n_shells`)
        #  or user-defined limits (`bounds`) if provided
        self._set_bounds(dr, n_shells, bounds)
#        # default normalization (r**2*g(r))
#        self.normalize = self.squared_distance_RDF_normalization

    @property
    def n_shells(self):
        """
        Upper bound for correlation expressed in number of coordinations shells.
        """
        return self._n_shells
    
    @n_shells.setter
    def n_shells(self, value):
        return self._set_bounds(self._dr, value, None)

    @property
    def bounds(self):
        """
        Lower and upper bounds to describe the radial correlations.
        """
        return self._bounds
    
    @bounds.setter
    def bounds(self, value):
        self._set_bounds(self._dr, None, value)
        
    @property
    def dr(self):
        """
        Grid spacing.
        """
        return self._dr
    
    @dr.setter
    def dr(self, value):
        self._set_bounds(value, self._n_shells, self._bounds)
        
    def compute(self):
        """
        Compute the radial correlations for the particles in group=0 in the 
        range of distances given by `bounds`.

        Returns
        -------
        features : numpy.ndarray
            Radial correlations.

        """
        StructuralDescriptor._sanity_checks(self)
        n_frames = len(self.groups[0])
        pos_0 = self.dump('position', 0)
        pos_1 = self.dump('position', 1)
        idx_0 = self.dump('index', 0)
        idx_1 = self.dump('index', 1)
        features = numpy.empty((self.size, self.n_features), dtype=numpy.int64)
        row = 0
        for n in range(n_frames):
            box = self.trajectory[n].cell.side
            # pos_x arrays need to be transposed to be used with fortran
            hist_n = compute.radial_histogram(pos_0[n].T, pos_1[n].T, 
                                              idx_0[n], idx_1[n], box,
                                              self.grid, self.dr)
            # fill the array of features
            for hist_n_i in hist_n:
                features[row] = hist_n_i
                row += 1
        self.features = features
        return features

    def normalize(self, distribution, method="r2"):
        """
        Later.
        """
        if method == "r2":
            # TODO: this normalization is a bit inconsistent
            return (self.grid*self.grid) * self.normalize(distribution, method="gr")
        elif method == "gr":
            rho = self.trajectory[0].density
            x_1 = self.group_fraction(1)
            if self.dimension == 2:
                const = numpy.pi * self.dr**2
            if self.dimension == 3:
                const = 4.0 / 3.0 * numpy.pi * self.dr**3
            const = const * rho * x_1
            g_b = numpy.empty_like(self.grid)
            b_min = numpy.floor(self._bounds[0] / self.dr) # if r_min != 0
            for m in range(self.n_features):
                b = b_min + m + 1
                wb = (b**3 - (b-1)**3)
                g_b[m] = distribution[m] / wb
            return g_b / const
        else:
            raise ValueError("unknown value {}".format(methos))
            
    #TODO: do not compute the g(r) on the whole trajectory only for one cutoff...
    #TODO: duplicate code with `compute()`
    def _set_bounds(self, dr, n_shells, bounds):
        # take the smallest side as maximal upper bound for the grid
        sides = numpy.array(self.trajectory.get_property('cell.side'))
        L = numpy.min(sides)
        # use `n_shells`
        self._dr = dr
        if bounds is None:
            # first define full grid
            r = numpy.arange(self._dr/2, L/2, self._dr, dtype=numpy.float64)
            self._bounds = (r[0], r[-1]) # temporary
            self.grid = r # temporary
            # arrays
            pos_0 = self.dump('position', 0)
            pos_1 = self.dump('position', 1)
            idx_0 = self.dump('index', 0)
            idx_1 = self.dump('index', 1)
            box = self.trajectory[0].cell.side
            all_hist = numpy.empty((self.size, r.size), dtype=numpy.int64)
            n_frames = len(self.groups[0])
            row = 0
            for n in range(n_frames):
                # pos_x arrays need to be transposed to be used with fortran
                hist_n = compute.radial_histogram(pos_0[n].T, pos_1[n].T, 
                                                  idx_0[n], idx_1[n], box,
                                                  r, self._dr)
                # fill the array of features
                for hist_n_i in hist_n:
                    all_hist[row] = hist_n_i
                    row += 1
            # g(r)
            g = numpy.sum(all_hist, axis=0)
            g = self.normalize(g, method="gr")
            # find position of the n-th minimum in g(r)
            index = 0
            for shell in range(n_shells):
                g_tmp = g[index:]
                first_max = numpy.argmax(g_tmp)
                first_min = numpy.argmin(g_tmp[first_max:]) + first_max
                index += first_min   
            # set grid and bounds
            self.grid = r[0:index+1]
            self._n_shells = n_shells
            self._bounds = (r[0], r[index])
        
        # use user-defined limits if provided
        else:
            if len(bounds) == 2 and bounds[0] < bounds[1] and bounds[1] <= L/2:
                rmin, rmax = bounds
                r = numpy.arange(rmin+(self._dr/2), rmax, self._dr, dtype=numpy.float64)
                # set grid and bounds
                self.grid = r
                self._n_shells = None
                self._bounds = (r[0], r[-1])
            else:
                raise ValueError('`bounds` is not correctly defined.')
