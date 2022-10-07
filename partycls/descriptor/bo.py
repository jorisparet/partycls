import numpy
from .descriptor import AngularStructuralDescriptor
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
    
    >>> D = BondOrientationalDescriptor('trajectory.xyz', orders=[4,6])
    >>> D.add_filter("species == 'A'", group=0)
    >>> D.compute()
    """

    name = 'bond-orientational'
    symbol = 'bo'

    def __init__(self, trajectory, lmin=1, lmax=8, orders=None, 
                 accept_nans=True, verbose=False):
        AngularStructuralDescriptor.__init__(self, trajectory,
                                             accept_nans=accept_nans,
                                             verbose=verbose)
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
        self._set_up(dtype=numpy.float64)
        self._manage_nearest_neighbors()
        self._filter_neighbors()
        n_frames = len(self.trajectory)
        # all relevant arrays
        pos_0 = self.dump('position', group=0)
        pos_all = self.trajectory.dump('position')
        box = self.trajectory.dump('cell.side')
        # computation
        start = 0
        for n in self._trange(n_frames):
            pos_0_n = pos_0[n].T
            pos_all_n = pos_all[n].T
            npart = len(self.groups[0][n])
            for ln, l in enumerate(self.grid):
                feat_n = compute.ql_all(l,
                                       self._neighbors[n],
                                       self._neighbors_number[n],
                                       pos_0_n, pos_all_n,
                                       box[n])
                self.features[start: start+npart, ln] = feat_n
            start += npart
        self._handle_nans()
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
    
    >>> D = LechnerDellagoDescriptor('trajectory.xyz', orders=[4,6])
    >>> D.add_filter("species == 'A'", group=0)
    >>> D.compute()
    """

    name = 'lechner-dellago'
    symbol = 'ld'

    def __init__(self, trajectory, lmin=1, lmax=8, orders=None, 
                 accept_nans=True, verbose=False):
        BondOrientationalDescriptor.__init__(self, trajectory, lmin=lmin,
                                             lmax=lmax, orders=orders, 
                                             accept_nans=accept_nans, verbose=verbose)

    def compute(self):
        # set up
        self._set_up(dtype=numpy.float64)
        self._manage_nearest_neighbors()
        self._filter_neighbors()
        self._filter_subsidiary_neighbors()
        n_frames = len(self.groups[0])
        row = 0
        # all relevant arrays
        pos_0 = self.dump('position', group=0)
        pos_all = self.trajectory.dump('position')
        box = self.trajectory.dump('cell.side')
        # computation
        for n in self._trange(n_frames):
            pos_all_n = pos_all[n].T
            for i in range(len(self.groups[0][n])):
                hist_n_i = numpy.empty_like(self.grid, dtype=numpy.float64)
                for ln, l in enumerate(self.grid):
                    hist_n_i[ln] = self._qbar_l(l,
                                                self._neighbors[n][i],
                                                self._subsidiary_neighbors[n][i],
                                                pos_0[n][i], pos_all_n, box[n])
                    # TODO: improve Fortran calculation for Lechner-Dellago
                    # hist_n_i[ln] = compute.qbarl(l, numpy.array(neigh_i),
                    #                               numpy.array(neigh_neigh_i).T,
                    #                               pos_0[n][i], pos_1[n].T, box)
                self.features[row] = hist_n_i
                row += 1
        self._handle_nans()
        return self.features

    def _qbar_lm(self, l, neigh_i, neigh_neigh_i, pos_i, pos_all, box):
        Nbar_b = len(neigh_i) + 1
        q_lm_i = compute.qlm(l, neigh_i, pos_i, pos_all, box)
        q_lm_k = []
        for kn in range(len(neigh_i)):
            k = neigh_i[kn]
            q_lm_k.append(compute.qlm(l, neigh_neigh_i[kn], pos_all[:,k], pos_all, box))
        qbar_lm = q_lm_i + numpy.sum(q_lm_k, axis=0)
        return qbar_lm / Nbar_b

    def _qbar_l(self, l, neigh_i, neigh_neigh_i, pos_i, pos_all, box):
        """
        Rotational invariant of order l for particle `i`.
        """
        qbar_lm = self._qbar_lm(l, neigh_i, neigh_neigh_i, pos_i, pos_all, box)
        return compute.rotational_invariant(l, qbar_lm)
