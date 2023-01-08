import numpy
from .descriptor import StructuralDescriptor

class PhysicsInformedDescriptor(StructuralDescriptor):
    """

    Physics-informed descriptor.
    
    ...
    ...
    ...

    .. math::
        X^\mathrm{C}(i) = (\: \Omega(i) \:) .
    
    See the tutorials for more details.

    Attributes
    ----------
    trajectory : Trajectory
        Trajectory on which the structural descriptor will be computed.
        
    active_filters : list
        All the active filters on both groups prior to the computation of the
        descriptor.
        
    dimension : int
        Spatial dimension of the descriptor (2 or 3).
        
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

    verbose : bool
        Show progress information and warnings about the computation of the 
        descriptor when verbose is ``True``, and remain silent when verbose is 
        ``False``.

    neighbors_boost : float, default: 1.5
        Scaling factor to estimate the number of neighbors relative to a
        an ideal gas with the same density. This is used internally to set
        the dimensions of lists of neighbors. A too small number creates a
        risk of overfilling the lists of neighbors, and a too large number
        increases memory usage. This only works if the associated ``Trajectory``
        has valid cutoffs in the ``Trajectory.nearest_neighbors_cutoffs`` list
        attribute. This sets the value of the ``max_num_neighbors`` attribute
        during the computation of the descriptor.

    max_num_neighbors : int, default: 100
        Maximum number of neighbors. This is used internally to set the dimensions
        of lists of neighbors. This number is automatically adjusted to limit
        memory usage if the associated ``Trajectory`` has valid cutoffs in the 
        ``Trajectory.nearest_neighbors_cutoffs`` list attribute. The
        default value ``100`` is used if no cutoffs can be used to estimate a
        better value. The default value is sufficient in most cases, otherwise 
        this number can manually be increased **before** computing the descriptor.

    """

    name = 'physics-informed'
    symbol = 'phinf'

    def __init__(self, trajectory, bounds=(1, 1.5), dr=0.1,
                 distance_grid=None, accept_nans=True, verbose=False):
        """
        Parameters
        ----------
        bounds : tuple, default: (1, 1.5)
            Lower and upper bounds :math:`(d_{n_\mathrm{min}}, d_{n_\mathrm{max}})`
            to define the grid of distances :math:`\{ d_n \}`, where consecutive 
            points in the the grid are separated by :math:`\Delta r`.
            
        dr : float, default: 0.1
            Grid spacing :math:`\Delta r` to define the grid of distances 
            :math:`\{ d_n \}` in the range 
            :math:`(d_{n_\mathrm{min}}, d_{n_\mathrm{max}})` by steps of size 
            :math:`\Delta r`.
            
        distance_grid : list, default: None
            Manually defined grid of distances :math:`\{ d_n \}`. If different 
            from ``None``, it  overrides the linearly-spaced grid defined by 
            ``bounds`` and ``dr``.        
        """
        StructuralDescriptor.__init__(self,
                                      trajectory,
                                      accept_nans=accept_nans,
                                      verbose=verbose)
        # dummy values, to be set in the following lines
        self._distance_grid = [-1]
        self._orders = [-1]
        # set the grid of distances {r}
        self._set_bounds(dr, bounds, distance_grid)

    @property
    def bounds(self):
        """
        Lower and upper bounds of the distance grid :math:`\{ d_n \}`.
        """
        return self._bounds

    @bounds.setter
    def bounds(self, value):
        self._set_bounds(self._dr, value, None)

    @property
    def dr(self):
        """
        Grid spacing :math:`\Delta r` to define the distance grid :math:`\{ d_n \}`.
        """
        return self._dr

    @dr.setter
    def dr(self, value):
        self._set_bounds(value, self._bounds, None)
    
    @property
    def distance_grid(self):
        """
        Grid of distances :math:`\{ d_n \}`.
        """
        return self._distance_grid
    
    @distance_grid.setter
    def distance_grid(self, value):
        self._set_bounds(self._dr, self._bounds, value)

    def compute(self):
        """
        Compute the physics-informed descriptor for the particles in group=0.
        Returns the data matrix and also overwrites the ``features`` attribute.

        Returns
        -------
        features : numpy.ndarray
            Data matrix.
        """
        # set up
        self._set_up(dtype=numpy.float64)
        n_frames = len(self.trajectory)
        # all relevant arrays
        ...
        # computation
        start = 0
        for n in self._trange(n_frames):
            npart = len(self.groups[0][n])
            for r, rn in enumerate(self._distance_grid):
                feat_r_n = ...
                ...
                self.features[start: start+npart, rn: (rn+1)*4] = feat_r_n
        self._handle_nans()
        return self.features

    def _set_bounds(self, dr, bounds, distance_grid):
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
                    r = numpy.arange(rmin, rmax + (self._dr / 2), self._dr, dtype=numpy.float64)
                    # set grid and bounds
                    self._distance_grid = r
                    self._bounds = (r[0], r[-1])
            else:
                raise ValueError('`bounds` is not correctly defined.')

        # update general grid
        self.grid = [[r for _ in range(4)] for r in self._distance_grid]
        
class JungDescriptor(PhysicsInformedDescriptor):
    """
    Alias for the class ``PhysicsInformedDescriptor``.
    """
    pass