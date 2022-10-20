import numpy
from .descriptor import StructuralDescriptor


class CoordinationDescriptor(StructuralDescriptor):
    """
    Coordination descriptor.
    
    The coordination number :math:`n_\\alpha(i)` of a particle :math:`i` is given
    by the number of its nearest neighbors whose chemical species 
    is :math:`\\alpha`, where  :math:`\\alpha` is either one of the :math:`n`
    chemical species :math:`\{ \\alpha_i \}_{i=1 \dots n}` in the trajectory 
    (*i.e.* partial coordination number, :math:`\\alpha = \\alpha_i`) or all 
    species at once (*i.e.* total coordination number, 
    :math:`\\alpha = \mathrm{all}`).

    The resulting **full** feature vector for particle :math:`i` is given by

    .. math::
        X^\mathrm{N}(i) = (\: n_\mathrm{all}(i) \;\;  n_{\\alpha_1}(i) \;\; \dots \;\; n_{\\alpha_n}(i) \:) ,
    
    but its size depends on whether the user requests the total coordination 
    number, the partial ones, or both.

    .. note::

        By applying a filter on ``group=1``, the total/partial coordination 
        number(s) can also be computed by considering and **arbitrary** subset 
        of particles (*e.g.* particles whose radius is smaller than a certain 
        value).

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
        Grid of chemical species for which the coordination number is computed.
        
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

    name = 'coordination'
    symbol = 'coord'

    def __init__(self, trajectory, total=True, partial=False,
                 accept_nans=True, verbose=False):
        """
        Parameters
        ----------
        trajectory : Trajectory
            Trajectory on which the structural descriptor will be computed.

        total : bool, default: True
            Compute the total coordination number.

        partial : bool, default: False
            Compute the coordination number :math:`n_{\\alpha_i}` for each chemical 
            species :math:`\\alpha_i` separately.

        accept_nans: bool, default: True
            If ``False``, discard any row from the array of features that contains a 
            `NaN` element. If ``True``, keep `NaN` elements in the array of features.

        verbose : bool, default: False
            Show progress information and warnings about the computation of the 
            descriptor when verbose is ``True``, and remain silent when verbose 
            is ``False``.
        """
        StructuralDescriptor.__init__(self, trajectory,
                                      accept_nans=accept_nans,
                                      verbose=verbose)
        self._total = total
        self._partial = partial
        self._set_grid(total, partial)

    @property
    def total(self):
        return self._total

    @total.setter
    def total(self, value):
        self._total = value
        self._set_grid(value, self._partial)

    @property
    def partial(self):
        return self._partial

    @partial.setter
    def partial(self, value):
        self._partial = value
        self._set_grid(self._total, value)

    def compute(self):
        """
        Compute the coordination number(ss for the particles in group=0.
        Returns the data matrix and also overwrites the ``features`` attribute.

        Returns
        -------
        features : numpy.ndarray
            Data matrix with the coordination number(s).
        """
        self._set_up(dtype=numpy.int64)
        self._manage_nearest_neighbors()
        self._filter_neighbors()
        n_frames = len(self.trajectory)
        row = 0
        # all relevant arrays
        distinct_species = self.trajectory[0].distinct_species
        spe_all = self.trajectory.dump('species')
        nn_spe_i = numpy.empty_like(self.grid, dtype=numpy.int64)
        offset = 1 if self._total else 0
        # computation
        for n in self._trange(n_frames):
            for i in range(len(self.groups[0][n])):
                nn_i = self._neighbors_number[n][i]
                # total
                if self._total:
                    nn_spe_i[0] = nn_i
                # partial
                if self._partial:
                    neigh_i = self._neighbors[n][i,0:nn_i]
                    for j, spe_j in enumerate(distinct_species):
                        nn_spe_i[j+offset] = numpy.count_nonzero(spe_all[n][neigh_i] == spe_j)
                self.features[row] = nn_spe_i
                row += 1
        self._handle_nans()
        return self.features

    def _set_grid(self, total, partial):
        """
        Set the grid of species.
        """
        # safety check
        if not (total or partial):
            raise ValueError("at least one of `total` or `partial` must be True.")
        # grid
        species = self.trajectory[0].distinct_species
        self.grid = []
        #  total
        if total:
            self.grid.append('all')
        #  partial
        if partial:
            for spe in species:
                self.grid.append(spe)
        self.grid = numpy.array(self.grid)
        # reset data matrix since grid has changed
        self.features = None
