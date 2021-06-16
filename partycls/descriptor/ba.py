import numpy
from .descriptor import StructuralDescriptor, AngularStructuralDescriptor
from .realspace_wrap import compute

class BondAngleDescriptor(AngularStructuralDescriptor):
    """
    Structural descriptor based on bond angles between particles.
    
    See the parent class for more details.
    
    Parameters
    ----------
    
    trajectory : str or an instance of `Trajectory`.
        Trajectory on which the structural descriptor will be computed.
        
    dtheta : float
        Bin width in degrees.
    
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
        List of cutoff distances to identify the nearest neighbors using
        the fixed-cutoff ('FC') method.
        
    nearest_neighbors_method : str, default: 'FC'
        Nearest neighbor method, 'FC' or 'SANN'.
    
    Examples:
    ---------
    
    >>> D = BondAngleDescriptor('trajectory.xyz', dtheta=2.0)
    >>> D.nearest_neighbors_method = 'SANN'
    >>> D.add_filter("species == 'A'")
    >>> D.compute()    
    """
    
    name = 'bond-angle'
    symbol = 'ba'
        
    def __init__(self, trajectory, dtheta=3.0):
        AngularStructuralDescriptor.__init__(self, trajectory)
        self._dtheta = dtheta
        self.grid = numpy.arange(dtheta/2.0, 180.0, dtheta, dtype=numpy.float64)
        
    @property
    def dtheta(self):
        return self._dtheta
    
    @dtheta.setter
    def dtheta(self, value):
        self._dtheta = value
        self.grid = numpy.arange(value/2.0, 180.0, value, dtype=numpy.float64)
    
    def compute(self):      
        StructuralDescriptor._sanity_checks(self)
        # all relevant arrays
        n_frames = len(self.groups[0])        
        pos_0 = self.dump('position', 0)
        pos_1 = self.dump('position', 1)
        idx_0 = self.dump('index', 0)       
        features = numpy.empty((self.size, self.n_features), dtype=numpy.int64)
        row = 0
        # compute nearest neighbors
        self.nearest_neighbors(method=self.nearest_neighbors_method)   
        for n in range(n_frames):
            box = self.trajectory[n].cell.side
            for i in range(len(idx_0[n])):  
                # individual bond-angle distribution using nearest-neighbors
                neigh_i = self.neighbors[n][i]
                hist_n_i = compute.angular_histogram(idx_0[n][i],
                                                     pos_0[n][i], pos_1[n].T,
                                                     neigh_i, box, 
                                                     self.n_features, 
                                                     self.dtheta)
                features[row] = hist_n_i
                row += 1   
        self.features = features
        return features

    def normalize(self, distribution, method="sin"):
        """
        Later.
        """
        if method == "sin":
            return distribution / numpy.sum(distribution) / self.dtheta
        elif method == "pdf":
            norm_sin = numpy.sin(self.grid / 360.0 * 2 * numpy.pi)
            distribution = distribution / norm_sin
            return distribution / numpy.sum(distribution) / self.dtheta
        else:
            raise ValueError("unknown value {}".format(method))
