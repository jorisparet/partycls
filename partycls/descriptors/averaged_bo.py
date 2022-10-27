import numpy
from .bo import BondOrientationalDescriptor
from .realspace_wrap import compute

class LocallyAveragedBondOrientationalDescriptor(BondOrientationalDescriptor):
    """
    Locally averaged bond-orientational descriptor.

    The complex coefficients :math:`q_{lm}(i)` of particle :math:`i` can be 
    averaged over its :math:`N_b(i)` nearest neighbors, as suggested by Lechner 
    and Dellago :cite:`lechner_2008`,

    .. math::
        \\bar{q}_{lm}(i) = \\frac{1}{N_b(i)+1} \left[ q_{l m}(i) + \\sum_{j=1}^{N_b(i)} q_{l m}(j) \\right],

    and then made invariant,

    .. math::
        \\bar{Q}_{l}(i) = \\sqrt{ \\frac{4\pi}{2l + 1}\\sum_{m=-l}^l |\\bar{q}_{lm(i)}|^2 } ,

    to provide an improved descriptor for crystal structure detection.

    We then consider :math:`\\bar{Q}_l(i)` for a sequence of orders 
    :math:`\{ l_n \} = \{ l_\mathrm{min}, \dots, l_\mathrm{max} \}`. The resulting 
    feature vector for particle :math:`i` is given by

    .. math::
        X^\mathrm{LABO}(i) = (\: \\bar{Q}_{l_\mathrm{min}}(i) \;\; \dots \;\; \\bar{Q}_{l_\mathrm{max}}(i) \:) .
    
    See the tutorials for more details.
    
    Attributes
    ----------
    trajectory : Trajectory
        Trajectory on which the structural descriptor will be computed.
        
    active_filters : list
        All the active filters on both groups prior to the computation of the
        descriptor.
        
    dimension : int
        Spatial dimension of the descriptor (2 or 3).
        
    grid : numpy.ndarray
        Grid of orders :math:`\{ l_n \}`.
        
    features : numpy.ndarray
        Array of all the structural features for the particles in group=0 in
        accordance with the defined filters (if any). This attribute is 
        initialized when the method ``compute`` is called (default value is ``None``).

    groups : tuple
        Composition of the groups: ``groups[0]`` and ``groups[1]`` contain lists of all
        the ``Particle`` instances in groups 0 and 1 respectively. Each element of 
        the tuple is a list of ``Particle`` in ``trajectory``, *e.g.* ``groups[0][0]``
        is the list of all the particles in the first frame of ``trajectory`` that 
        belong to group=0.

    verbose : bool
        Show progress information and warnings about the computation of the 
        descriptor when verbose is ``True``, and remain silent when verbose is 
        ``False``.

    neighbors_boost : float, default: 1.5
        Scaling factor to estimate the number of neighbors relative to a
        an ideal gas with the same density. This is used internally to set
        the dimensions of lists of neighbors. A too small number creates a
        risk of overfilling the lists of neighbors, and a too large number
        increases memory usage. This only works if the associated ``Trajectory``
        has valid cutoffs in the ``Trajectory.nearest_neighbors_cutoffs`` list
        attribute. This sets the value of the ``max_num_neighbors`` attribute
        during the computation of the descriptor.

    max_num_neighbors : int, default: 100
        Maximum number of neighbors. This is used internally to set the dimensions
        of lists of neighbors. This number is automatically adjusted to limit
        memory usage if the associated ``Trajectory`` has valid cutoffs in the 
        ``Trajectory.nearest_neighbors_cutoffs`` list attribute. The
        default value ``100`` is used if no cutoffs can be used to estimate a
        better value. The default value is sufficient in most cases, otherwise 
        this number can manually be increased **before** computing the descriptor.
    """

    name = 'locally averaged bond-orientational'
    symbol = 'labo'

    def __init__(self, trajectory, lmin=1, lmax=8, orders=None, 
                 accept_nans=True, verbose=False):
        """
        Parameters
        ----------
        trajectory : Trajectory
            Trajectory on which the structural descriptor will be computed.
            
        lmin : int, default: 1
            Minimum order :math:`l_\mathrm{min}`. This sets the lower bound of 
            the grid :math:`\{ l_n \}`.
            
        lmax : int, default: 8
            Maximum order :math:`l_\mathrm{max}`. This sets the upper bound of 
            the grid :math:`\{ l_n \}`. For numerical reasons, 
            :math:`l_\mathrm{max}` cannot be larger than 16.
            
        orders: list, default: None
            Sequence :math:`\{l_n\}` of specific orders to compute, *e.g.* 
            ``orders=[4,6]``. This has the priority over ``lmin`` and ``lmax``.
        
        accept_nans: bool, default: True
            If ``False``, discard any row from the array of features that contains a 
            `NaN` element. If ``True``, keep `NaN` elements in the array of features.

        verbose : bool, default: False
            Show progress information and warnings about the computation of the 
            descriptor when verbose is ``True``, and remain silent when verbose 
            is ``False``.
        """
        BondOrientationalDescriptor.__init__(self, trajectory, lmin=lmin,
                                             lmax=lmax, orders=orders, 
                                             accept_nans=accept_nans, verbose=verbose)

    def compute(self):
        """
        Compute the locally averaged bond-orientational correlations for the particles
        in group=0 for the grid of orders :math:`\{ l_n \}`. Returns the data matrix 
        and also overwrites the ``features`` attribute.

        Returns
        -------
        features : numpy.ndarray
            Data matrix with bond-orientational correlations.
        """
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
                nn_i = self._neighbors_number[n][i]
                for ln, l in enumerate(self.grid):
                    hist_n_i[ln] = self._qbar_l(l,
                                                self._neighbors[n][i][0:nn_i],
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

class LechnerDellagoDescriptor(LocallyAveragedBondOrientationalDescriptor):
    """
    Alias for the class ``AveragedBondOrientationalDescriptor``.
    """
    pass