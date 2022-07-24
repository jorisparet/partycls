import numpy
from .descriptor import StructuralDescriptor, AngularStructuralDescriptor
from .realspace_wrap import compute

class TetrahedralDescriptor(AngularStructuralDescriptor):

    name = 'bond-angle'
    symbol = 'tetra'
    
    def __init__(self, trajectory):
        AngularStructuralDescriptor.__init__(self, trajectory)
        self.grid = numpy.zeros(1, dtype=numpy.float64)
        
    def compute(self):
        StructuralDescriptor._sanity_checks(self)
        # all relevant arrays
        n_frames = len(self.groups[0])
        pos_0 = self.dump('position', 0)
        pos_1 = self.dump('position', 1)
        idx_0 = self.dump('index', 0)
        features = numpy.empty((self.size, self.n_features), dtype=numpy.float64)
        row = 0
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
                features[row] = tetra_i
                row += 1
        self.features = features
        return features    