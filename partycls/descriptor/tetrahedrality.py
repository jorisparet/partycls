import numpy
from .descriptor import AngularStructuralDescriptor
from .realspace_wrap import compute

class TetrahedralDescriptor(AngularStructuralDescriptor):
    """
    Tetrahedral descriptor.
    
    Computes the average deviation of the bond angles between a central 
    particle i and all the pairs of its nearest neighbors (j,k) from the ideal 
    angle in a tetrahedron (109.5°).
    
    This descriptor is scalar. Therefore, the `grid` attribute is not relevant.
    
    See the parent class for more details.
    
    Parameters
    ----------
    
    trajectory : str or an instance of `Trajectory`.
        Trajectory on which the structural descriptor will be computed.
    
    accept_nans: bool, default: True
        If False, discard any row from the array of features that contains a Nan
        element. If True, keep NaN elements in the array of features.

    verbose : bool, default: False
        Show progress information and warnings about the computation of the 
        descriptor when verbose is True, and remain silent when verbose is False.

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
    
    >>> D = TetrahedralDescriptor('trajectory.xyz')
    >>> D.add_filter("species == 'A'", group=0)
    >>> D.compute()
    """
    
    name = 'tetrahedral'
    symbol = 'tetra'
    
    def __init__(self, trajectory, accept_nans=True, verbose=False):
        AngularStructuralDescriptor.__init__(self, trajectory,
                                             accept_nans=accept_nans,
                                             verbose=verbose)
        self.grid = numpy.zeros(1, dtype=numpy.float64)
        
    def compute(self):
        self._set_up(dtype=numpy.float64)
        self._manage_nearest_neighbors()
        self._filter_neighbors()
        n_frames = len(self.trajectory)
        # all relevant arrays
        pos_0 = self.dump('position', group=0)
        pos_all = self.trajectory.dump('position')
        idx_0 = self.dump('_index', group=0)
        box = self.trajectory.dump('cell.side')
        # computation
        start = 0
        for n in self._trange(n_frames):
            pos_0_n = pos_0[n].T
            pos_all_n = pos_all[n].T
            npart = len(self.groups[0][n])
            feat_n = compute.tetrahedrality_all(idx_0[n],
                                                pos_0_n, pos_all_n,
                                                self._neighbors[n],
                                                self._neighbors_number[n],
                                                box[n])
            self.features[start: start+npart, 0] = feat_n
            start += npart
        self._handle_nans()
        return self.features