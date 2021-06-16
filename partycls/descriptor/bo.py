import numpy
from .descriptor import StructuralDescriptor, AngularStructuralDescriptor
from .realspace_wrap import compute

class BondOrientationalDescriptor(AngularStructuralDescriptor):
    """
    Structural descriptor based on bond order parameters as defined by
    Steinhardt et al. (https://doi.org/10.1103%2FPhysRevB.28.784).
    
    See the parent class for more details.
    
    Parameters
    ----------
    
    trajectory : str or an instance of `Trajectory`.
        Trajectory on which the structural descriptor will be computed.
        
    lmin : int, default: 0
        Minimum degree. This set the lower bound of the grid.
        
    lmax : int, default: 8
        Minimum degree. This set the upper bound of the grid.
        
    orders: list, default: None
        Specific values of orders to compute, e.g. orders=[4,6]. This has
        the priority over `lmin` and `lmax`.
    
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
    
    >>> D = BondOrientationalDescriptor('trajectory.xyz', orders=[4,6])
    >>> D.nearest_neighbors_method = 'FC'
    >>> D.compute()
    """

    name = 'bond-orientational'
    symbol = 'bo'
    
    def __init__(self, trajectory, lmin=1, lmax=8, orders=None):
        AngularStructuralDescriptor.__init__(self, trajectory)
        if self.trajectory[0].n_dimensions == 2:
            raise ValueError('trajectory must be 3-dimensional to be used with a {} descriptor'.format(self.name))
        self._bounds(lmin, lmax, orders)

    @property    
    def orders(self):
        return self.grid
    
    @orders.setter
    def orders(self, values):
        self._bounds(1, 8, values)
        
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
                # compute BO parameters for particle `i`
                neigh_i = self.neighbors[n][i]
                hist_n_i = numpy.empty_like(self.grid, dtype=numpy.float64)
                for ln, l in enumerate(self.grid):
                    hist_n_i[ln] = compute.ql(l, neigh_i, pos_0[n][i], pos_1[n].T, box)
                features[row] = hist_n_i
                row += 1      
        self.features = features
        return features

    def _bounds(self, lmin, lmax, orders):
        if orders is None:
            self.grid = numpy.array(range(lmin, lmax+1))
        else:
            self.grid = numpy.array(orders)
    
class LechnerDellagoDescriptor(BondOrientationalDescriptor):
    """
    Structural descriptor based on locally averaged bond order parameters as
    defined by Lechner & Dellago (https://doi.org/10.1063/1.2977970).
    
    See the parent class for more details.
    
    Parameters
    ----------
    
    trajectory : str or an instance of `Trajectory`.
        Trajectory on which the structural descriptor will be computed.
        
    lmin : int, default: 0
        Minimum degree. This set the lower bound of the grid.
        
    lmax : int, default: 8
        Minimum degree. This set the upper bound of the grid.
        
    orders: list, default: None
        Specific values of orders to compute, e.g. orders=[4,6]. This has
        the priority over `lmin` and `lmax`.
    
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
    
    >>> D = LechnerDellagoDescriptor('trajectory.xyz', orders=[4,6])
    >>> D.nearest_neighbors_method = 'FC'
    >>> D.add_filter("species == 'A'")
    >>> D.compute()
    """

    name = 'lechner-dellago'
    symbol = 'ld'
    
    def __init__(self, trajectory, lmin=1, lmax=8, orders=None):
        BondOrientationalDescriptor.__init__(self, trajectory, lmin=lmin, 
                                             lmax=lmax, orders=orders)

    def compute(self):
        StructuralDescriptor._sanity_checks(self)
        # all relevant arrays
        n_frames = len(self.groups[0])        
        pos_0 = self.dump('position', 0)
        pos_1 = self.dump('position', 1)
        idx_0 = self.dump('index', 0)
        idx_1 = self.dump('index', 1)
        spe_1 = self.dump('species_id', 1)
        pairs = numpy.asarray(self.trajectory[0].pairs_of_species_id)
        features = numpy.empty((self.size, self.n_features), dtype=numpy.float64)
        row = 0
        # compute nearest neighbors
        self.nearest_neighbors(method=self.nearest_neighbors_method)
        cutoffs = numpy.array(self.cutoffs)
        for n in range(n_frames):
            box = self.trajectory[n].cell.side
            for i in range(len(idx_0[n])):
                # neighbors of neighbors of i
                neigh_i = self.neighbors[n][i]
                neigh_neigh_i = []
                # TODO: SANN for neighbors of j if set as nearest_neighbors_method
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
                    hist_n_i[ln] = self._qbar_l(l, neigh_i, neigh_neigh_i, pos_0[n][i], pos_1[n], box)
                    #TODO: improve Fortran calculation for Lechner-Dellago
                    # hist_n_i[ln] = compute.qbarl(l, numpy.array(neigh_i), 
                    #                               numpy.array(neigh_neigh_i).T, 
                    #                               pos_0[n][i], pos_1[n].T, box)
                features[row] = hist_n_i
                row += 1      
        self.features = features
        return features
    
    def _qbar_lm(self, l, neigh_i, neigh_neigh_i, pos_i, pos_j, box):
        Nbar_b = len(neigh_i) + 1
        q_lm_i = compute.qlm(l, neigh_i, pos_i, pos_j.T, box)
        q_lm_k = []
        for kn in range(len(neigh_i)):
            k = neigh_i[kn]
            q_lm_k.append(compute.qlm(l, neigh_neigh_i[kn], pos_j[k], pos_j.T, box))            
        qbar_lm = q_lm_i + numpy.sum(q_lm_k, axis=0)
        return qbar_lm / Nbar_b
    
    def _qbar_l(self, l, neigh_i, neigh_neigh_i, pos_i, pos_j, box):
        """
        Rotational invariant of order l for particle `i`.
        """
        qbar_lm = self._qbar_lm(l, neigh_i, neigh_neigh_i, pos_i, pos_j, box)
        return compute.rotational_invariant(l, qbar_lm)
