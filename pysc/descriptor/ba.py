from .descriptor import AngularStructuralDescriptor
from .helpers import cartesian_to_spherical, pbc
import numpy
from .realspace_wrap import compute

class AngularDescriptor(AngularStructuralDescriptor):
    
    name = 'angular'
    symbol = 'ba'
        
    def __init__(self, trajectory, dtheta=3.0):
        AngularStructuralDescriptor.__init__(self, trajectory)
        self.dtheta = dtheta
        self.grid = numpy.arange(dtheta/2.0, 180.0, dtheta, dtype=numpy.float64)
        
    @property
    def n_features(self):
        return len(self.grid)
    
    def compute(self):
        # compute cutoffs if not provided
        if None in self.cutoffs:
            self._compute_cutoffs()
        # all relevant arrays
        n_frames = len(self._groups[0])
        idx_0 = self.group_indices(0)
        idx_1 = self.group_indices(1)
        idx_all = self.trajectory.dump('index')
        spe_0 = self.group_species_id(0)
        spe_1 = self.group_species_id(1)
        pos_0 = self.group_positions(0)
        pos_1 = self.group_positions(1)
        pos_all = self.trajectory.dump('position')
        pairs = numpy.asarray(self.trajectory[0].pairs_of_species_id)
        box = self.trajectory[0].cell.side
        cutoffs = numpy.asarray(self.cutoffs)
        rmax = 1.5*max(cutoffs)
        features = numpy.empty((self.size, self.n_features), dtype=numpy.int64)
        row = 0
        for n in range(n_frames):
            for i in range(len(idx_0[n])):
                
                # Find nearest neighbors
                #  - fixed cutoff (FC)
                if self.nearest_neighbors_method == 'FC':
                    neigh_i = compute.nearest_neighbors(idx_0[n][i], numpy.array(idx_1[n]),
                                                        numpy.array(pos_0[n][i]), numpy.array(pos_1[n]).T,
                                                        spe_0[n][i], numpy.array(spe_1[n]),
                                                        pairs, box, cutoffs)
                    
                #  - Solid Angle Nearest Neighbors (SANN)
                #  This will find all neighbors of `i`
                #  (including particles not in group=1)
                if self.nearest_neighbors_method == 'SANN':
                    neigh_i = compute.sann(numpy.array(pos_0[n][i]), numpy.array(pos_all[n]).T,
                                           idx_0[n][i], numpy.array(idx_all[n]), numpy.array(idx_1[n]),
                                           rmax, box)
                    
                # remove -1 spaces (non-neighbors) in the array
                # /!\ indices in the `neigh_i` array are not the same
                #      as self.group_indices()
                neigh_i = neigh_i[neigh_i >= 0]

                # individual bond-angle distribution using nearest-neighbors
                hist_n_i = compute.angular_histogram(idx_0[n][i],
                                                     numpy.array(pos_0[n][i]), numpy.array(pos_1[n]).T,
                                                     neigh_i, box, self.n_features, self.dtheta)
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
