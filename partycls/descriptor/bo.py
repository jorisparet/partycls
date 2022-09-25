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
        
    lmin : int, default: 1
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
    
    Examples:
    ---------
    
    >>> D = BondOrientationalDescriptor('trajectory.xyz', orders=[4,6])
    >>> D.add_filter("species == 'A'", group=0)
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
        # set up
        StructuralDescriptor._set_up(self, dtype=numpy.float64)
        AngularStructuralDescriptor._manage_nearest_neighbors(self)
        AngularStructuralDescriptor._filter_neighbors(self)
        n_frames = len(self.trajectory)
        row = 0
        # all relevant arrays
        pos_0 = self.dump('position', group=0)
        pos_1 = self.dump('position', group=1)
        idx_0 = self.dump('index', group=0)
        # computation
        for n in range(n_frames):
            box = self.trajectory[n].cell.side
            for i in range(len(idx_0[n])):
                hist_n_i = numpy.empty_like(self.grid, dtype=numpy.float64)
                for ln, l in enumerate(self.grid):
                    hist_n_i[ln] = compute.ql(l, self._neighbors[n][i], pos_0[n][i], pos_1[n].T, box)
                self.features[row] = hist_n_i
                row += 1
        return self.features

    def _bounds(self, lmin, lmax, orders):
        if orders is None:
            self.grid = numpy.array(range(lmin, lmax + 1))
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
        
    lmin : int, default: 1
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
    
    Examples:
    ---------
    
    >>> D = LechnerDellagoDescriptor('trajectory.xyz', orders=[4,6])
    >>> D.add_filter("species == 'A'", group=0)
    >>> D.compute()
    """

    name = 'lechner-dellago'
    symbol = 'ld'

    def __init__(self, trajectory, lmin=1, lmax=8, orders=None):
        BondOrientationalDescriptor.__init__(self, trajectory, lmin=lmin,
                                             lmax=lmax, orders=orders)

    def compute(self):
        # set up
        StructuralDescriptor._set_up(self, dtype=numpy.float64)
        AngularStructuralDescriptor._manage_nearest_neighbors(self)
        AngularStructuralDescriptor._filter_neighbors(self)
        AngularStructuralDescriptor._filter_subsidiary_neighbors(self)
        n_frames = len(self.groups[0])
        row = 0
        # all relevant arrays
        pos_0 = self.dump('position', group=0)
        pos_1 = self.dump('position', group=1)
        idx_0 = self.dump('index', group=0)
        row = 0
        # computation
        for n in range(n_frames):
            box = self.trajectory[n].cell.side
            for i in range(len(idx_0[n])):
                hist_n_i = numpy.empty_like(self.grid, dtype=numpy.float64)
                for ln, l in enumerate(self.grid):
                    hist_n_i[ln] = self._qbar_l(l,
                                                self._neighbors[n][i],
                                                self._subsidiary_neighbors[n][i],
                                                pos_0[n][i], pos_1[n], box)
                    # TODO: improve Fortran calculation for Lechner-Dellago
                    # hist_n_i[ln] = compute.qbarl(l, numpy.array(neigh_i),
                    #                               numpy.array(neigh_neigh_i).T,
                    #                               pos_0[n][i], pos_1[n].T, box)
                self.features[row] = hist_n_i
                row += 1
        return self.features

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
