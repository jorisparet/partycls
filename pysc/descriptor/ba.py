from .descriptor import AngularStructuralDescriptor
import numpy
from .realspace_wrap import compute

class BondAngleDescriptor(AngularStructuralDescriptor):
    
    name = 'bond-angle'
    symbol = 'ba'
        
    def __init__(self, trajectory, dtheta=3.0):
        AngularStructuralDescriptor.__init__(self, trajectory)
        self.dtheta = dtheta
        self.grid = numpy.arange(dtheta/2.0, 180.0, dtheta, dtype=numpy.float64)
        
    @property
    def n_features(self):
        return len(self.grid)
    
    def compute(self):         
        # all relevant arrays
        n_frames = len(self._groups[0])
        idx_0 = self.group_indices(0)
        pos_0, pos_1 = self.group_positions(0), self.group_positions(1)
        box = self.trajectory[0].cell.side
        features = numpy.empty((self.size, self.n_features), dtype=numpy.int64)
        row = 0
        # compute nearest neighbors
        self.nearest_neighbors(method=self.nearest_neighbors_method)   
        for n in range(n_frames):
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

    def normalize(self, dist):
        """
        Later.
        """
        return dist * (1.0 / numpy.sum(dist) / self.dtheta)
    
    #TODO: check normalization (directly include sinus?)
    def normalize_sin(self, dist):
        """
        Later.
        """
        norm_sin = numpy.sin(self.grid / 360.0 * 2 * numpy.pi)
        dist = dist / norm_sin
        return dist * (1.0 / numpy.sum(dist) / self.dtheta)
