import numpy
from .descriptor import StructuralDescriptor, AngularStructuralDescriptor
from .realspace_wrap import compute

class BondOrientationalDescriptor(AngularStructuralDescriptor):

    name = 'bond-orientational'
    symbol = 'bo'
    
    def __init__(self, trajectory, lmin=0, lmax=8, orders=None):
        AngularStructuralDescriptor.__init__(self, trajectory)
        if self.trajectory[0].number_of_dimensions == 2:
            raise ValueError('trajectory must be 3-dimensional to be used with a {} descriptor'.format(self.name))
        self._bounds(lmin, lmax, orders)

    @property    
    def orders(self):
        return self.grid
    
    @orders.setter
    def orders(self, values):
        self._bounds(0, 8, values)
        
    def compute(self):
        StructuralDescriptor.sanity_checks(self)
        # all relevant arrays
        n_frames = len(self._groups[0])
        idx_0 = self.group_indices(0)
        pos_0 = self.group_positions(0)
        pos_1 = self.group_positions(1)
        features = numpy.empty((self.size, self.n_features), dtype=numpy.float64)
        row = 0
        # compute nearest neighbors
        self.nearest_neighbors(method=self.nearest_neighbors_method)
        for n in range(n_frames):
            box = self.trajectory[n].cell.side
            for i in range(len(idx_0[n])):
                # compute BO parameters for particle `i`
                neigh_i = self.neighbors[n][i]
                hist_n_i = numpy.empty_like(self.grid, dtype=numpy.float64)
                for ln, l in enumerate(self.grid):
                    hist_n_i[ln] = compute.ql(l, neigh_i, pos_0[n][i], pos_1[n].T, box)
                features[row] = hist_n_i
                row += 1      
        self.features = features
        return features
    
    def normalize(self, dist):
        return dist

    def _bounds(self, lmin, lmax, orders):
        if orders is None:
            self.grid = numpy.array(range(lmin, lmax+1))
        else:
            self.grid = numpy.array(orders)
    
class LechnerDellagoDescriptor(BondOrientationalDescriptor):

    name = 'lechner-dellago'
    symbol = 'ld'
    
    def __init__(self, trajectory, lmin=0, lmax=8, orders=None):
        BondOrientationalDescriptor.__init__(self, trajectory, lmin=lmin, 
                                             lmax=lmax, orders=orders)

    def compute(self):
        StructuralDescriptor.sanity_checks(self)
        # all relevant arrays
        n_frames = len(self._groups[0])
        idx_0 = self.group_indices(0)
        idx_1 = self.group_indices(1)
        spe_1 = self.group_species_id(1)
        pos_0 = self.group_positions(0)
        pos_1 = self.group_positions(1)
        pairs = numpy.asarray(self.trajectory[0].pairs_of_species_id)
        cutoffs = numpy.asarray(self.cutoffs)
        features = numpy.empty((self.size, self.n_features), dtype=numpy.float64)
        row = 0
        # compute nearest neighbors
        self.nearest_neighbors(method=self.nearest_neighbors_method)
        for n in range(n_frames):
            box = self.trajectory[n].cell.side
            for i in range(len(idx_0[n])):
                # neighbors of neighbors of i
                neigh_i = self.neighbors[n][i]
                neigh_neigh_i = []
                for j in neigh_i:
                    neigh_j = compute.nearest_neighbors(idx_1[n][j], idx_1[n],
                                                        pos_1[n][j], pos_1[n].T,
                                                        spe_1[n][j], spe_1[n],
                                                        pairs, box, cutoffs)
                    neigh_j = numpy.asarray(neigh_j[neigh_j >= 0])
                    neigh_neigh_i.append(neigh_j)
                # compute BO parameters for particle `i`
                hist_n_i = numpy.empty_like(self.grid, dtype=numpy.float64)
                for ln, l in enumerate(self.grid):
                    hist_n_i[ln] = self._qbar_l(l, i, neigh_i, neigh_neigh_i, pos_0[n][i], pos_1[n], box)
                    #TODO: fix the shape of neigh_neigh_i to pass it to qbarl()
                    #hist_n_i[l] = compute.qbarl(l, numpy.array(neigh_i), numpy.array(neigh_neigh_i).T, pos_0[n][i], pos_1[n].T, box)
                features[row] = hist_n_i
                row += 1      
        self.features = features
        return features
    
    def _qbar_lm(self, l, i, neigh_i, neigh_neigh_i, pos_i, pos_j, box):
        Nbar_b = len(neigh_i) + 1
        q_lm_i = compute.qlm(l, neigh_i, pos_i, pos_j.T, box)
        q_lm_k = []
        for kn in range(len(neigh_i)):
            k = neigh_i[kn]
            q_lm_k.append(compute.qlm(l, neigh_neigh_i[kn], pos_j[k], pos_j.T, box))            
        qbar_lm = q_lm_i + numpy.sum(q_lm_k, axis=0)
        return qbar_lm / Nbar_b
    
    def _qbar_l(self, l, i, neigh_i, neigh_neigh_i, pos_i, pos_j, box):
        """
        Rotational invariant of order l for particle `i`.
        """
        qbar_lm = self._qbar_lm(l, i, neigh_i, neigh_neigh_i, pos_i, pos_j, box)
        return compute.rotational_invariant(l, qbar_lm)
