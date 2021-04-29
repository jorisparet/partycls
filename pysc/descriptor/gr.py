from .descriptor import StructuralDescriptor
import numpy
from .realspace_wrap import compute

class RadialDescriptor(StructuralDescriptor):

    name = 'radial'
    symbol = 'gr'
    
    def __init__(self, trajectory, dr=0.1, n_shells=3, rlim=None):
        StructuralDescriptor.__init__(self, trajectory)
        self.n_shells = n_shells
        # set the grid automatically using coordination shells (`n_shells`)
        #  or user-defined limits (`rlim`) if provided
        self._bounds(dr, rlim)
#        # default normalization (r**2*g(r))
#        self.normalize = self.squared_distance_RDF_normalization

    @property
    def rlim(self):
        return self._rlim
    
    @rlim.setter
    def rlim(self, value):
        self._bounds(self._dr, value)
        
    @property
    def dr(self):
        return self._dr
    
    @dr.setter
    def dr(self, value):
        self._bounds(value, self._rlim)
            
    @property
    def n_features(self):
        return len(self.grid)
        
    def compute(self):
        StructuralDescriptor.sanity_checks(self)
        n_frames = len(self._groups[0])
        pos_0 = self.group_positions(0)
        pos_1 = self.group_positions(1)
        idx_0 = self.group_indices(0)
        idx_1 = self.group_indices(1)
        box = self.trajectory[0].cell.side
        features = numpy.empty((self.size, self.n_features), dtype=numpy.int64)
        row = 0
        for n in range(n_frames):
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
    
    def normalize_gr(self, dist):
        """
        Later.
        """
        rho = self.trajectory[0].density
        x_1 = self.group_fraction(1)
        if self.dimension == 2:
            const = numpy.pi * self.dr**2
        if self.dimension == 3:
            const = 4.0 / 3.0 * numpy.pi * self.dr**3
        const = const * rho * x_1
        g_b = numpy.empty_like(self.grid)
        b_min = numpy.floor(self._rlim[0]/self.dr) # if r_min != 0
        for m in range(self.n_features):
            b = b_min + m + 1
            wb = (b**3 - (b-1)**3)
            g_b[m] = dist[m] / wb
        return g_b / const

    def normalize(self, dist):
        """
        Later.
        """
        return (self.grid*self.grid) * self.normalize_gr(dist)
#        rho = self.trajectory[0].density
#        x_1 = self.group_fraction(1)
#        if self.dimension == 2:
#            const = numpy.pi * self.dr**2
#        if self.dimension == 3:
#            const = 4.0 / 3.0 * numpy.pi * self.dr**3
#        const = const * rho * x_1
#        g_b = numpy.empty_like(self.grid)
#        b_min = numpy.floor(self.rlim[0]/self.dr) # if r_min != 0
#        for m in range(self.n_features):
#            b = b_min + m + 1
#            wb = (b**3 - (b-1)**3)
#            g_b[m] = dist[m] / wb
#        return self.grid**2 * g_b / const
            
    #TODO: do not compute the g(r) on the whole trajectory only for one cutoff...
    #TODO: duplicate code with `compute()`
    def _bounds(self, dr, rlim):
        L = numpy.min(self.trajectory[0].cell.side)
        # use `n_shells`
        self._dr = dr
        if rlim is None:
            # first define full grid
            r = numpy.arange(self._dr/2, L/2, self._dr, dtype=numpy.float64)
            self._rlim = (r[0], r[-1]) # temporary
            self.grid = r # temporary
            # arrays
            pos_0 = self.group_positions(0)
            pos_1 = self.group_positions(1)
            idx_0 = self.group_indices(0)
            idx_1 = self.group_indices(1)
            box = self.trajectory[0].cell.side
            all_hist = numpy.empty((self.size, r.size), dtype=numpy.int64)
            n_frames = len(self._groups[0])
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
            g = self.normalize_gr(g)
            # find position of the n-th minimum in g(r)
            index = 0
            for shell in range(self.n_shells):
                g_tmp = g[index:]
                first_max = numpy.argmax(g_tmp)
                first_min = numpy.argmin(g_tmp[first_max:]) + first_max
                index += first_min   
            # set grid and bounds
            self.grid = r[0:index+1]
            self._rlim = (r[0], r[index])
        
        # use user-defined limits if provided
        else:
            if len(rlim) == 2 and rlim[0] < rlim[1] and rlim[1] <= L/2:
                rmin, rmax = rlim
                r = numpy.arange(rmin+(self._dr/2), rmax, self._dr, dtype=numpy.float64)
                # set grid and bounds
                self.grid = r
                self._rlim = (r[0], r[-1])
            else:
                raise ValueError('`rlim` is not correctly defined.')
