from .descriptor import AngularStructuralDescriptor
#from .helpers import cartesian_to_spherical, pbc
import numpy
#from scipy.special import sph_harm
from .realspace_wrap import compute

class BondOrientationalDescriptor(AngularStructuralDescriptor):

    name = 'bond-orientational'
    symbol = 'bo'
    
    def __init__(self, trajectory, lmax=10):
        AngularStructuralDescriptor.__init__(self, trajectory)
        if self.trajectory[0].number_of_dimensions == 2:
            raise ValueError('trajectory must be 3-dimensional to be used with a {} descriptor'.format(self.name))
        self.grid = numpy.array(range(0, lmax))

    @property
    def n_features(self):
        return len(self.grid)  
        
    def compute(self):
        # all relevant arrays
        n_frames = len(self._groups[0])
        idx_0 = self.group_indices(0)
        pos_0 = self.group_positions(0)
        pos_1 = self.group_positions(1)
        box = self.trajectory[0].cell.side
        features = numpy.empty((self.size, self.n_features), dtype=numpy.float64)
        row = 0
        # compute nearest neighbors
        self.nearest_neighbors(method=self.nearest_neighbors_method)
        for n in range(n_frames):
            for i in range(len(idx_0[n])):
                # compute BO parameters for particle `i`
                neigh_i = self.neighbors[n][i]
                hist_n_i = numpy.empty_like(self.grid, dtype=numpy.float64)
                for l in self.grid:
#                    hist_n_i[l] = self._q_l(l, i, neigh_i, pos_0[n][i], pos_1[n], box)
                    hist_n_i[l] = compute.ql(l, neigh_i, pos_0[n][i], pos_1[n].T, box)
                features[row] = hist_n_i
                row += 1      
        self.features = features
        return features

#    def _rotational_invariant(self, l, q_lm):
#        s = numpy.sum(q_lm*q_lm.conj(), axis=0)
#        s = numpy.real(s)
#        q_l = numpy.sqrt( 4.0*numpy.pi / (2*l+1) * s )
#        return q_l
    
#    def _q_l(self, l, i, neigh_i, pos_i, pos_j, box):
#        """
#        Rotational invariant of order l for particle `i`.
#        """
#        q_lm = self._q_lm(l, i, neigh_i, pos_i, pos_j, box)
#        return self._rotational_invariant(l, q_lm)
    
#    def _q_lm(self, l, i, neigh_i, pos_i, pos_j, box):
#        q_lm = numpy.zeros(2*l+1, dtype=complex)
#        # r_ij (cartesian)
#        r_xyz = pos_j[neigh_i] - pos_i
#        r_xyz = pbc(r_xyz, box)
#        # r_ij (spherical)
#        r_sph = cartesian_to_spherical(r_xyz)
#        for m in range(2*l+1):
#            Y_lm = sph_harm(m-l, l, r_sph[:,1], r_sph[:,2])
#            q_lm[m] = numpy.average(Y_lm)
#        return q_lm
    
#    def normalize(self, dist):
#        return dist * (1.0 / numpy.sum(dist))
    
    def normalize(self, dist):
        return dist
    
class LechnerDellagoDescriptor(BondOrientationalDescriptor):

    name = 'lechner-dellago'
    symbol = 'ld'
    
    def __init__(self, trajectory, lmax=10):
        BondOrientationalDescriptor.__init__(self, trajectory, lmax=lmax)

    def compute(self):
        # all relevant arrays
        n_frames = len(self._groups[0])
        idx_0 = self.group_indices(0)
        idx_1 = self.group_indices(1)
        spe_1 = self.group_species_id(1)
        pos_0 = self.group_positions(0)
        pos_1 = self.group_positions(1)
        pairs = numpy.asarray(self.trajectory[0].pairs_of_species_id)
        box = self.trajectory[0].cell.side
        cutoffs = numpy.asarray(self.cutoffs)
        features = numpy.empty((self.size, self.n_features), dtype=numpy.float64)
        row = 0
        # compute nearest neighbors
        self.nearest_neighbors(method=self.nearest_neighbors_method)
        for n in range(n_frames):
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
                for l in self.grid:
                    hist_n_i[l] = self._qbar_l(l, i, neigh_i, neigh_neigh_i, pos_0[n][i], pos_1[n], box)
                    #TODO: fix the shape of neigh_neigh_i to pass it to qbarl()
                    #hist_n_i[l] = compute.qbarl(l, numpy.array(neigh_i), numpy.array(neigh_neigh_i).T, pos_0[n][i], pos_1[n].T, box)
                features[row] = hist_n_i
                row += 1      
        self.features = features
        return features
    
    def _qbar_lm(self, l, i, neigh_i, neigh_neigh_i, pos_i, pos_j, box):
        Nbar_b = len(neigh_i) + 1
#        q_lm_i = self._q_lm(l, i, neigh_i, pos_i, pos_j, box)
        q_lm_i = compute.qlm(l, neigh_i, pos_i, pos_j.T, box)
        q_lm_k = []
        for kn in range(len(neigh_i)):
            k = neigh_i[kn]
#            q_lm_k.append(self._q_lm(l, k, neigh_neigh_i[kn], pos_j[k], pos_j, box))            
            q_lm_k.append(compute.qlm(l, neigh_neigh_i[kn], pos_j[k], pos_j.T, box))            
        qbar_lm = q_lm_i + numpy.sum(q_lm_k, axis=0)
        return qbar_lm / Nbar_b
    
    def _qbar_l(self, l, i, neigh_i, neigh_neigh_i, pos_i, pos_j, box):
        """
        Rotational invariant of order l for particle `i`.
        """
        qbar_lm = self._qbar_lm(l, i, neigh_i, neigh_neigh_i, pos_i, pos_j, box)
        return compute.rotational_invariant(l, qbar_lm)
#        return self._rotational_invariant(l, qbar_lm)
