import numpy
from .descriptor import StructuralDescriptor
from .realspace_wrap import compute


class BondOrientationalDescriptor(StructuralDescriptor):
    """
    Bond-orientational descriptor.

    Bond-order parameters :cite:`steinhardt_1983` are standard measures of 
    structure in the first coordination shell. Let :math:`\mathbf{r}_i` be the 
    position of particle :math:`i` and define 
    :math:`\mathbf{r}_{ij} = \mathbf{r}_j - \mathbf{r}_i` and 
    :math:`r_{ij} = |\mathbf{r}_{ij}|`. Then consider the weighted microscopic 
    density around particle :math:`i`:

    .. math::
        \\rho(\mathbf{r}; i) = \\sum_{j=1}^{N_b(i)} w_j \delta(\mathbf{r} - \mathbf{r}_{ij})

    where :math:`w_j` is a particle-dependent weight and the sum involves a set of
    :math:`N_b(i)` particles, which defines the coordination shell of interest for 
    particle :math:`i`.

    We project the microscopic density on a unit-radius sphere, that is, 
    :math:`\hat{\\rho}(\hat{\mathbf{r}}; i) = \\sum_{j=1}^{N_b(i)} w_j \delta(\mathbf{r} - \hat{\mathbf{r}}_{ij})`,
    where :math:`\hat{\mathbf{r}} = \mathbf{r} / |\mathbf{r}|` and similarly 
    :math:`\hat{\mathbf{r}}_{ij} = \mathbf{r}_{ij}/|\mathbf{r}_{ij}|`. Expanding 
    in spherical harmonics yields

    .. math::
        \hat{\\rho}(\hat{\mathbf{r}}; i) = \\sum_{l=0}^\infty \sum_{m=-l}^l c_{l m}(i) Y_{l m}(\hat{\mathbf{r}}) ,

    with coefficients

    .. math::
        c_{l m}(i) =  \int d\mathbf{r} \\rho(\mathbf{r}; i) Y_{l m}(\hat{\mathbf{r}}) .

    In the conventional bond-order analysis, one sets the weights :math:`w_j` to 
    unity and considers the normalized complex coefficients,

    .. math::
        \\begin{align}
        q_{lm}(i) & = \\frac{1}{N_b(i)} \int d\mathbf{r} \\rho(\mathbf{r}; i) Y_{l m}(\hat{\mathbf{r}}) 
        \\nonumber \\\ & = \\frac{1}{N_b(i)} \\sum_{j=1}^{N_b(i)} Y_{l m}(\hat{\mathbf{r}}_{ij}) .
        \end{align}

    The rotational invariants,

    .. math::
        Q_{l}(i) = \\sqrt{ \\frac{4\pi}{2l + 1}\\sum_{m=-l}^l |q_{lm}(i)|^2 },

    provide a detailed structural description of the local environment around 
    particle :math:`i`.


    We then consider :math:`Q_l(i)` for a sequence of orders 
    :math:`\{ l_n \} = \{ l_\mathrm{min}, \dots, l_\mathrm{max} \}`. The resulting 
    feature vector for particle :math:`i` is given by

    .. math::
        X^\mathrm{BO}(i) = (\: Q_{l_\mathrm{min}}(i) \;\; \dots \;\; Q_{l_\mathrm{max}}(i) \:) .
    
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

    name = 'bond-orientational'
    symbol = 'bo'

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
        StructuralDescriptor.__init__(self, trajectory,
                                      accept_nans=accept_nans,
                                      verbose=verbose)
        self._dimension_check(dimension=3)
        self._bounds(lmin, lmax, orders)

    @property
    def orders(self):
        """
        Grid of orders :math:`\{ l_n \}`.
        """
        return self.grid

    @orders.setter
    def orders(self, values):
        self._bounds(1, 8, values)

    def compute(self):
        """
        Compute the bond-orientational correlations for the particles in group=0
        for the grid of orders :math:`\{ l_n \}`. Returns the data matrix and also
        overwrites the ``features`` attribute.

        Returns
        -------
        features : numpy.ndarray
            Data matrix with bond-orientational correlations.
        """
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
            self.grid = numpy.sort(orders)
        # check lmax
        self._check_lmax(self.grid)

    def _check_lmax(self, grid):
        if max(grid) > 16:
            raise ValueError("the largest possible value for an order l is 16.")
        

class SteinhardtDescriptor(BondOrientationalDescriptor):
    """
    Alias for the class ``BondOrientationalDescriptor``.
    """
    pass