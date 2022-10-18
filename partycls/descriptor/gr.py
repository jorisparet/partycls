import numpy
from .descriptor import StructuralDescriptor
from .realspace_wrap import compute


class RadialDescriptor(StructuralDescriptor):
    """
    Radial descriptor.

    Let :math:`\mathbf{r}_i` be the position of particle :math:`i` and define 
    :math:`r_{ij} = |\mathbf{r}_j - \mathbf{r}_i|` as the distance between 
    particle :math:`i` and its neighbors :math:`j`.

    We define :math:`n_i(r_m)` as the number of neighbors :math:`j` of particle 
    :math:`i` for which :math:`r_{ij}` is between 
    :math:`r_m = r_\mathrm{min} + m \\times \Delta r` and 
    :math:`r_{m+1} = r_\mathrm{min} + (m+1) \\times \Delta r` :cite:`paret_2020`. 
    Here, :math:`\Delta r` has the interpration of a bin width in a histogram.

    We then consider :math:`n_i(r_m)` for a set of distances :math:`\{ d_n \}` 
    separated by :math:`\Delta r`, 
    :math:`\{ d_n \} = \{ r_0, r_1, \dots, r_{n_\mathrm{max}} \}`. The resulting 
    feature vector for particle :math:`i` is given by

    .. math::
        X^\mathrm{R}(i) = (\: n_i(r_0) \;\; n_i(r_1) \;\; \dots \;\; n_i(r_{n_\mathrm{max}}) \:) .
    
    See tutorials for more details.

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
        Grid of distances :math:`\{ d_n \}`.
        
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
    """

    name = 'radial'
    symbol = 'gr'

    def __init__(self, trajectory, dr=0.1, n_shells=3, bounds=None,
                 accept_nans=True, verbose=False):
        """
        Parameters
        ----------
        trajectory : Trajectory
            Trajectory on which the structural descriptor will be computed.
            
        dr : float
            Bin width :math:`\Delta r`.
            
        n_shells : int, default: 3
            Number of coordination shells (based on the :math:`g(r)` of group=0). 
            This sets the upper bound :math:`r_\mathrm{max}` for the distance grid 
            :math:`\{d_n\}` up to which correlations are computed.
            
        bounds : tuple, default: None
            Lower and upper bounds :math:`(r_\mathrm{min},r_\mathrm{max})` to describe
            the radial correlations. If set, this has the priority over ``n_shells``.

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
        # set the grid automatically using coordination shells (`n_shells`)
        #  or user-defined limits (`bounds`) if provided
        self._set_bounds(dr, n_shells, bounds)
#        # default normalization (r**2*g(r))
#        self.normalize = self.squared_distance_RDF_normalization

    @property
    def n_shells(self):
        """
        Upper bound in the grid of distances :math:`\{d_n\}` expressed in number of 
        coordinations shells.
        """
        return self._n_shells

    @n_shells.setter
    def n_shells(self, value):
        return self._set_bounds(self._dr, value, None)

    @property
    def bounds(self):
        """
        Lower and upper bounds :math:`(r_\mathrm{min},r_\mathrm{max})` to describe
        the radial correlations.
        """
        return self._bounds

    @bounds.setter
    def bounds(self, value):
        self._set_bounds(self._dr, None, value)

    @property
    def dr(self):
        """
        Grid spacing :math:`\Delta r` for the grid of distances :math:`\{ d_n \}`.
        """
        return self._dr

    @dr.setter
    def dr(self, value):
        self._set_bounds(value, self._n_shells, self._bounds)

    def compute(self):
        """
        Compute the radial correlations for the particles in group=0
        for the grid of distances :math:`\{ d_n \}`. Returns the data matrix and also
        overwrites the ``features`` attribute.

        Returns
        -------
        features : numpy.ndarray
            Data matrix with radial correlations.
        """
        # set up
        self._set_up(dtype=numpy.int64)
        n_frames = len(self.trajectory)
        row = 0
        # all relevant arrays
        pos_0 = self.dump('position', group=0)
        pos_1 = self.dump('position', group=1)
        idx_0 = self.dump('_index', group=0)
        idx_1 = self.dump('_index', group=1)
        # computation
        for n in self._trange(n_frames):
            pos_0_n = pos_0[n].T
            pos_1_n = pos_1[n].T
            box = self.trajectory[n].cell.side
            hist_n = compute.radial_histogram(pos_0_n, pos_1_n,
                                              idx_0[n], idx_1[n], box,
                                              self.grid, self.dr)
            for hist_n_i in hist_n:
                self.features[row] = hist_n_i
                row += 1
        self._handle_nans()
        return self.features

    def normalize(self, distribution, method="r2"):
        """
        Normalize a radial distribution.

        Parameters
        ----------
        distribution : numpy.ndarray
            Distribution to normalize.
            
        method : str, default: "r2"
            Normalization method:

            - ``method='r2'``: returns math:`r^2 \\times g(r)` (default);
            - ``method='gr'`` : returns the standard :math:`g(r)` ;

        Raises
        ------
        ValueError
            If ``method`` is invalid.

        Returns
        -------
        numpy.ndarray
            Normalized distribution.
        """
        if method == "r2":
            # TODO: this normalization is a bit inconsistent
            return (self.grid * self.grid) * self.normalize(distribution, method="gr")
        elif method == "gr":
            rho = self.trajectory[0].density
            x_1 = self.group_fraction(1)
            if self.dimension == 2:
                const = numpy.pi * self.dr**2
            if self.dimension == 3:
                const = 4.0 / 3.0 * numpy.pi * self.dr**3
            const = const * rho * x_1
            g_b = numpy.empty_like(self.grid)
            b_min = numpy.floor(self._bounds[0] / self.dr)  # if r_min != 0
            for m in range(self.n_features):
                b = b_min + m + 1
                wb = (b**3 - (b - 1)**3)
                g_b[m] = distribution[m] / wb
            return g_b / const
        else:
            raise ValueError("unknown value {}".format(method))

    # TODO: do not compute the g(r) on the whole trajectory only for one cutoff...
    # TODO: duplicate code with `compute()`
    def _set_bounds(self, dr, n_shells, bounds):
        # take the smallest side as maximal upper bound for the grid
        sides = numpy.array(self.trajectory.get_property('cell.side'))
        L = numpy.min(sides)
        # use `n_shells`
        self._dr = dr
        if bounds is None:
            # first define full grid
            r = numpy.arange(self._dr / 2, L / 2, self._dr, dtype=numpy.float64)
            self._bounds = (r[0], r[-1])  # temporary
            self.grid = r  # temporary
            # arrays
            pos_0 = self.dump('position', 0)
            pos_1 = self.dump('position', 1)
            idx_0 = self.dump('_index', 0)
            idx_1 = self.dump('_index', 1)
            box = self.trajectory[0].cell.side
            all_hist = numpy.empty((self.n_samples, r.size), dtype=numpy.int64)
            n_frames = len(self.groups[0])
            row = 0
            for n in range(n_frames):
                # pos_x arrays need to be transposed to be used with fortran
                hist_n = compute.radial_histogram(pos_0[n].T, pos_1[n].T,
                                                  idx_0[n], idx_1[n], box,
                                                  r, self._dr)
                # fill the array of features
                for hist_n_i in hist_n:
                    all_hist[row] = hist_n_i
                    row += 1
            # g(r)
            g = numpy.sum(all_hist, axis=0)
            g = self.normalize(g, method="gr")
            # find position of the n-th minimum in g(r)
            index = 0
            for shell in range(n_shells):
                g_tmp = g[index:]
                first_max = numpy.argmax(g_tmp)
                first_min = numpy.argmin(g_tmp[first_max:]) + first_max
                index += first_min
            # set grid and bounds
            self.grid = r[0:index + 1]
            self._n_shells = n_shells
            self._bounds = (r[0], r[index])

        # use user-defined limits if provided
        else:
            if len(bounds) == 2 and bounds[0] < bounds[1] and bounds[1] <= L / 2:
                rmin, rmax = bounds
                r = numpy.arange(rmin + (self._dr / 2), rmax, self._dr, dtype=numpy.float64)
                # set grid and bounds
                self.grid = r
                self._n_shells = None
                self._bounds = (r[0], r[-1])
            else:
                raise ValueError('`bounds` is not correctly defined.')
