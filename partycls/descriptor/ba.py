import numpy
from .descriptor import AngularStructuralDescriptor
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
        Grid over which the structural features will be computed.
        
    features : ndarray
        Array of all the structural features for the particles in group=0 in
        accordance with the defined filters (if any). This attribute is 
        initialized when the method `compute` is called (default value is None).
    
    Examples:
    ---------
    
    >>> D = BondAngleDescriptor('trajectory.xyz', dtheta=2.0)
    >>> D.add_filter("species == 'A'", group=0)
    >>> D.compute()    
    """

    name = 'bond-angle'
    symbol = 'ba'

    def __init__(self, trajectory, dtheta=3.0, accept_nans=True, verbose=0):
        AngularStructuralDescriptor.__init__(self, trajectory,
                                             verbose=verbose,
                                             accept_nans=accept_nans)
        self._dtheta = dtheta
        self.grid = numpy.arange(dtheta / 2.0, 180.0, dtheta, dtype=numpy.float64)

    @property
    def dtheta(self):
        return self._dtheta

    @dtheta.setter
    def dtheta(self, value):
        self._dtheta = value
        self.grid = numpy.arange(value / 2.0, 180.0, value, dtype=numpy.float64)

    def compute(self):
        # set up
        self._set_up(dtype=numpy.int64)
        self._manage_nearest_neighbors()
        self._filter_neighbors()
        n_frames = len(self.trajectory)
        row = 0
        # all relevant arrays
        pos_0 = self.dump('position', group=0)
        pos_all = self.trajectory.dump('position')
        idx_0 = self.dump('_index', group=0)
        box = self.trajectory.dump('cell.side')
        # computation
        for n in self._trange(n_frames):
            pos_all_n = pos_all[n].T
            for i in range(len(self.groups[0][n])):
                hist_n_i = compute.angular_histogram(idx_0[n][i],
                                                     pos_0[n][i], pos_all_n,
                                                     self._neighbors[n][i], box[n],
                                                     self.n_features,
                                                     self.dtheta)
                self.features[row] = hist_n_i
                row += 1
        self._handle_nans()
        return self.features

    def normalize(self, distribution, method="sin"):
        """
        Normalize a bond angle distribution.

        Parameters
        ----------
        
        distribution : array
            Distribution to normalize.
            
        method : str, optional
            Normalization method:
            - method='sin': by construction, the probability density of
            has a sinusoidal enveloppe in 3D for uniformly distributed points 
            on a sphere (default) ;
            - method='pdf' : gives a flat probability density for uniformly 
            distributed points on a sphere ;

        Raises
        ------
        
        ValueError
            If `method` is invalid.

        Returns
        -------
        
        array
            Normalized distribution.

        """
        if method == "sin":
            return distribution / numpy.sum(distribution) / self.dtheta
        elif method == "pdf":
            norm_sin = numpy.sin(self.grid / 360.0 * 2 * numpy.pi)
            distribution = distribution / norm_sin
            return distribution / numpy.sum(distribution) / self.dtheta
        else:
            raise ValueError("unknown value {}".format(method))
