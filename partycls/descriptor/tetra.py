import numpy
from .descriptor import StructuralDescriptor, AngularStructuralDescriptor
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
    
    >>> D = TetrahedralDescriptor('trajectory.xyz')
    >>> D.nearest_neighbors_method = 'FC'
    >>> D.compute()
    """
    
    name = 'tetrahedral'
    symbol = 'tetra'
    
    def __init__(self, trajectory):
        AngularStructuralDescriptor.__init__(self, trajectory)
        self.grid = numpy.zeros(1, dtype=numpy.float64)
        
    def compute(self):
        StructuralDescriptor._set_up(self, dtype=numpy.float64)
        n_frames = len(self.trajectory)
        row = 0
        # all relevant arrays
        pos_0 = self.dump('position', group=0)
        pos_1 = self.dump('position', group=1)
        idx_0 = self.dump('index', group=0)
        # compute nearest neighbors
        self.nearest_neighbors(method=self.nearest_neighbors_method)
        for n in range(n_frames):
            box = self.trajectory[n].cell.side
            for i in range(len(idx_0[n])):
                # individual bond-angle distribution using nearest-neighbors
                neigh_i = self.neighbors[n][i]
                tetra_i = compute.tetrahedrality(idx_0[n][i],
                                                 pos_0[n][i], pos_1[n].T,
                                                 neigh_i, box)
                self.features[row] = tetra_i
                row += 1
        return self.features