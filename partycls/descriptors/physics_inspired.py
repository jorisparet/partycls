import numpy
from .descriptor import StructuralDescriptor
from .realspace_wrap import compute

class PhysicsInspiredDescriptor(StructuralDescriptor):
    """

    Physics-informed descriptor.
    
    ...
    ...
    ...

    .. math::
        X^\mathrm{C}(i) = (\: \:) .
    
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

    name = 'physics-inspired'
    symbol = 'phins'

    def __init__(self, trajectory, bounds=(1, 1.5), dr=0.1,
                 distance_grid=None,
                 energy_field='potential_energy',
                 cell_field='cell_perimeter',
                 accept_nans=True, verbose=False):
        """
        Parameters
        ----------
        bounds : tuple, default: (1, 1.5)
            Lower and upper bounds :math:`(d_{n_\mathrm{min}}, d_{n_\mathrm{max}})`
            to define the grid of distances :math:`\{ d_n \}`, where consecutive 
            points in the the grid are separated by :math:`\Delta r`.
            
        dr : float, default: 0.1
            Grid spacing :math:`\Delta r` to define the grid of distances 
            :math:`\{ d_n \}` in the range 
            :math:`(d_{n_\mathrm{min}}, d_{n_\mathrm{max}})` by steps of size 
            :math:`\Delta r`.
            
        distance_grid : list, default: None
            Manually defined grid of distances :math:`\{ L_n \}`. If different 
            from ``None``, it  overrides the linearly-spaced grid defined by 
            ``bounds`` and ``dr``.

        energy_field : str, default: "potential_energy"
            Name of the particle property used to store the potential energy.

        cell_field : str, default: "cell_perimeter"
            Name of the particle property used to store the appropriate Voronoi
            cell metric (*i.e.* the cell perimeter in 2D and the cell surface in
            3D).

        accept_nans: bool, default: True
            If ``False``, discard any row from the array of features that contains a 
            `NaN` element. If ``True``, keep `NaN` elements in the array of features.

        verbose : bool, default: False
            Show progress information and warnings about the computation of the 
            descriptor when verbose is ``True``, and remain silent when verbose 
            is ``False``.
        """
        StructuralDescriptor.__init__(self,
                                      trajectory,
                                      accept_nans=accept_nans,
                                      verbose=verbose)
        # dummy values, to be set in the following lines
        #self._distance_grid = [-1]
        # set the grid of distances {r}
        self._set_bounds(dr, bounds, distance_grid)
        # fields
        self.energy_field = energy_field
        self.cell_field = cell_field

    @property
    def bounds(self):
        """
        Lower and upper bounds of the distance grid :math:`\{ L_n \}`.
        """
        return self._bounds

    @bounds.setter
    def bounds(self, value):
        self._set_bounds(self._dr, value, None)

    @property
    def dr(self):
        """
        Grid spacing :math:`\Delta r` to define the distance grid :math:`\{ L_n \}`.
        """
        return self._dr

    @dr.setter
    def dr(self, value):
        self._set_bounds(value, self._bounds, None)
    
    @property
    def distance_grid(self):
        """
        Grid of distances :math:`\{ L_n \}`.
        """
        return self._distance_grid
    
    @distance_grid.setter
    def distance_grid(self, value):
        self._set_bounds(self._dr, self._bounds, value)

    def compute(self):
        """
        Compute the physics-informed descriptor for the particles in group=0.
        Returns the data matrix and also overwrites the ``features`` attribute.

        Returns
        -------
        features : numpy.ndarray
            Data matrix.
        """
        # set up
        self._set_up(dtype=numpy.float64)
        n_frames = len(self.trajectory)
        # all relevant arrays
        pos_0 = self.dump('position', group=0)
        pos_1 = self.dump('position', group=1)
        spe_1 = self.dump('species_id', group=1)
        Epot_1 = self.dump('{}'.format(self.energy_field), group=1)
        cell_1 = self.dump('{}'.format(self.cell_field), group=1)
        box = self.trajectory.dump('cell.side')
        # the descriptor is computed for each species ID separately (1, 2, ...)
        # and also for all species at once, which is expressed as -1.
        distinct_species = list(self.trajectory[0].distinct_species_id)
        distinct_species += [-1]
        # computation
        start = 0
        n_L = len(self._distance_grid)
        for n in self._trange(n_frames):
            pos_0_n = pos_0[n].T
            pos_1_n = pos_1[n].T
            spe_1_n = spe_1[n]
            Epot_1_n = Epot_1[n]
            cell_1_n = cell_1[n]
            npart = len(self.groups[0][n])
            for spe_n, beta in enumerate(distinct_species):
                feat_spe_n = compute.physics_inspired_all(self._distance_grid,
                                                           pos_0_n, pos_1_n,
                                                           spe_1_n, beta,
                                                           Epot_1_n, cell_1_n,
                                                           box[n])
                self.features[start: start+npart, 4*n_L*spe_n: 4*n_L*(spe_n+1) ] = feat_spe_n
            start += npart
        self._handle_nans()
        return self.features


    # def compute(self):
    #     """
    #     Compute the physics-informed descriptor for the particles in group=0.
    #     Returns the data matrix and also overwrites the ``features`` attribute.

    #     Returns
    #     -------
    #     features : numpy.ndarray
    #         Data matrix.
    #     """
    #     # set up
    #     self._set_up(dtype=numpy.float64)
    #     n_frames = len(self.trajectory)
    #     # all relevant arrays
    #     pos_0 = self.dump('position', group=0)
    #     pos_all = self.trajectory.dump('position')
    #     spe_0 = self.dump('species_id', group=0)
    #     spe_all = self.trajectory.dump('species_id')
    #     Epot_0 = self.dump('particle.{}'.format(self.potential_energy_field), group=0)
    #     Epot_all = self.trajectory.dump('particle.{}'.format(self.potential_energy_field))
    #     cell_perimeter_0 = self.dump('particle.{}'.format(self.cell_perimeter_field), group=0)
    #     cell_perimeter_all = self.trajectory.dump('particle.{}'.format(self.cell_perimeter_field))
    #     box = self.trajectory.dump('cell.side')
    #     # compute extended neighbors with extended cutoffs
    #     # based on the largest distance in the distance grid
    #     R_cut = self.distance_grid[-1] + self.skin
    #     n_pairs = len(self.trajectory[0].pairs_of_species)
    #     extended_cutoffs = numpy.array([R_cut for _ in range(n_pairs)])
    #     self._compute_extended_neighbors(extended_cutoffs)
    #     # the descriptor is computed for each species separately and
    #     # also for all species at once.
    #     # we use two separate lists to keep track of this:
    #     # (1) a list of species_id [1, 2,, 3, ..., -1] (-1 = all).
    #     # (2) a list of boolean [False, False, ..., True] to know
    #     # if all species must be considered or not.
    #     distinct_species = list(self.trajectory[0].distinct_species_id)
    #     use_all_species = [False for _ in range(len(distinct_species))]
    #     distinct_species += [-1]
    #     use_all_species += [True]
    #     # computation
    #     start = 0
    #     for n in self._trange(n_frames):
    #         pos_0_n = pos_0[n].T
    #         pos_all_n = pos_all[n].T
    #         spe_0_n = spe_0[n]
    #         spe_all_n = spe_all[n]
    #         Epot_0_n = Epot_0[n]
    #         Epot_all_n = Epot_all[n]
    #         cell_perimeter_0_n = cell_perimeter_0[n]
    #         cell_perimeter_all_n = cell_perimeter_all[n]
    #         npart = len(self.groups[0][n])
    #         count_L_spe = 0
    #         M = numpy.max(self._extended_neighbors_number[n])
    #         for L in self._distance_grid:
    #             for spe_target, use_all in zip(distinct_species, use_all_species):
    #                 feat_L_spe_n = compute.physics_inspired_all(L,
    #                                                         self._extended_neighbors[n],
    #                                                         self._extended_neighbors_number[n],
    #                                                         pos_0_n,
    #                                                         pos_all_n,
    #                                                         spe_0_n,
    #                                                         spe_all_n,
    #                                                         spe_target,
    #                                                         use_all,
    #                                                         Epot_0_n,
    #                                                         Epot_all_n,
    #                                                         cell_perimeter_0_n,
    #                                                         cell_perimeter_all_n,
    #                                                         M,
    #                                                         box[n])
    #                 self.features[start: start+npart, 4*count_L_spe: 4*(count_L_spe+1)] = feat_L_spe_n
    #                 count_L_spe += 1
    #         start += npart
    #     self._handle_nans()
    #     return self.features

    def _set_bounds(self, dr, bounds, distance_grid):
        # take the smallest side as maximal upper bound for the distance grid
        sides = numpy.array(self.trajectory.get_property('cell.side'))
        L = numpy.min(sides)
        
        if distance_grid is not None:
            self._distance_grid = numpy.array(distance_grid, dtype=numpy.float64)
            self._dr = None
            self._bounds = (self._distance_grid[0], self._distance_grid[-1])
        
        else:
            self._dr = dr
            if len(bounds) == 2 and bounds[0] >= 0 and bounds[1] >= 0 and bounds[0] < bounds[1] and bounds[1] <= L / 2:
                    rmin, rmax = bounds
                    r = numpy.arange(rmin, rmax + (self._dr / 2), self._dr, dtype=numpy.float64)
                    # set grid and bounds
                    self._distance_grid = r
                    self._bounds = (r[0], r[-1])
            else:
                raise ValueError('`bounds` is not correctly defined.')

        # update general grid
        self._set_grid()

    def _set_grid(self):
        observables = ['DEN', 'EPOT', 'PERI', 'EPOT_VAR2']
        species = list(self.trajectory[0].distinct_species) + ['ALL']
        grid = []
        for L in self._distance_grid:
            for spe in species:
                for obs in observables:
                    grid.append( ('length={}'.format(L), 'type={}'.format(spe), 'obs={}'.format(obs)) )
        self.grid = grid
        
class JungDescriptor(PhysicsInspiredDescriptor):
    """
    Alias for the class ``PhysicsInformedDescriptor``.
    """
    pass