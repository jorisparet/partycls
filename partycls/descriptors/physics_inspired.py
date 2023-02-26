import numpy
from .descriptor import StructuralDescriptor
from .realspace_wrap import compute

class PhysicsInspiredDescriptor(StructuralDescriptor):
    """

    Physics-inspired descriptor.
    
    The physics-inspired descriptor is a collection of established structural
    order parameters (observables), first introduced by *Jung et al.* 
    :cite:`jung_2023` to train *GlassMLP*, a deep neural network used to predict
    the long-time dynamics in deeply supercooled liquids.

    For each particle :math:`i`, four distinct observables are computed and 
    coarse-grained over a set of lengthscales :math:`\{L_n\}`:

    1. Local density: 

    .. math::
        \\bar{\\rho}_i(L_n;\\beta_m) = \sum_{j \in N_{\\beta_m}^i} \exp(-\\frac{r_{ij}}{L_n}) ,

    where the sum runs over all :math:`N_{\\beta_m}^i` particles of type
    :math:`\\beta_m` within distance :math:`r_{ij} < R_c` from particle :math:`i`
    (including :math:`i` itself if its type is :math:`\\beta_m`). Particles :math:`j`
    beyond the cutoff distance :math:`R_c` will be ignored for efficiency reasons.

    2. Potential energy:

    .. math::
        \\bar{E}_i(L_n;\\beta_m) = \\frac{1}{\\bar{\\rho}_i(L_n;\\beta_m)} \sum_{j \in N_{\\beta_m}^i} E_j \exp(-\\frac{r_{ij}}{L_n}) ,
        
    where :math:`E_j` is extracted from the pair potential :math:`V(r)`,
    such that :math:`E_j = \sum_{j \\neq k} V(r_{jk})`.

    3. Perimeter or surface of the Voronoi cell centered on particle :math:`i`
    (2- and 3-dimensional cases, respectively):

    .. math::
        \\bar{\mathcal{C}}_i(L_n;\\beta_m) = \\frac{1}{\\bar{\\rho}_i(L_n;\\beta_m)} \sum_{j \in N_{\\beta_m}^i} \mathcal{C}_j \exp(- \\frac{r_{ij}}{L_n}) .

    4. Local variance of the potential energy:

    .. math::
        \mathrm{Var}( \\bar{E}_i(L_n;\\beta_m) ) = \\frac{1}{\\bar{\\rho}_i(L_n;\\beta_m)} \sum_{j \in N_{\\beta_m}^i} \left( E_j - \\bar{E}_i(L_n;\\beta_m) \\right)^2 \exp(-\\frac{r_{ij}}{L_n}) .

    The prefactor :math:`\\frac{1}{\\bar{\\rho}_i(L_n;\\beta_m)}` is used for
    normalization purposes. These four observables are then arranged in

    .. math::
        \\bar{S}_i(L_n;\\beta_m) = \{ \\bar{\\rho}_i(L_n;\\beta_m), \\bar{E}_i(L_n;\\beta_m), \\bar{\mathcal{C}}_i(L_n;\\beta_m),  \mathrm{Var}( \\bar{E}_i(L_n;\\beta_m) )\}
        
    and computed for each distance in :math:`\{ L_n \}` and for each particle
    type in :math:`\{ \\beta_m \}`, including one additional
    computation where **all** particles are considered at once regardless of
    their type. The total number :math:`M` of features produced is given by
    :math:`M = 4 \\times N_L \\times (N_\\beta + 1)`, where :math:`N_L` is the total
    number of coarse-graining lengths and :math:`N_\\beta` is the total number of
    distinct particle types in the system.

    The resulting feature vector for particle :math:`i` is given by

    .. math::
        X^\mathrm{P}(i) = (\: \\bar{S}_i(L_1;\\beta_1) \;\; \dots \;\; \\bar{S}_i(L_n;\\beta_m) \;\; \dots \;\;  \\bar{S}_i(L_{N_L};\\beta_{N_{\\beta}+1} ) \:) .
    
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

    def __init__(self, trajectory, bounds=(1.0, 2.0), dr=0.5,
                 distance_grid=None,
                 cutoff=None,
                 energy_field='potential_energy',
                 cell_field='cell_surface',
                 accept_nans=True, verbose=False):
        """
        Parameters
        ----------
        bounds : tuple, default: (1, 2)
            Lower and upper bounds :math:`(L_\mathrm{min}, L_\mathrm{max})`
            to define the grid of distances :math:`\{ L_n \}`, where consecutive 
            points in the the grid are separated by :math:`\Delta r`.
            
        dr : float, default: 0.5
            Grid spacing :math:`\Delta r` to define the grid of distances 
            :math:`\{ L_n \}` in the range 
            :math:`(L_\mathrm{min}, L_\mathrm{max})` by steps of size 
            :math:`\Delta r`.
            
        distance_grid : list, default: None
            Manually defined grid of distances :math:`\{ L_n \}`. If different 
            from ``None``, it  overrides the linearly-spaced grid defined by 
            ``bounds`` and ``dr``.

        cutoff : float, default: None
            Distance :math:`R_c` up to which the particles contribute to the coarse-grained
            average of the central particle. Particles beyond this cutoff distance
            will be ignored, for efficiency reasons. This distance should be larger
            than the largest distance in the grid :math:`\{ L_n \}` (large enough not to
            discard significant contributions). If ``None``, the value will automatically
            be chosen to be half the side of the simulation cell.

        energy_field : str, default: "potential_energy"
            Name of the particle property used to store the potential energy :math:`E_i`.

        cell_field : str, default: "cell_surface"
            Name of the particle property used to store the appropriate Voronoi
            cell metric :math:`\mathcal{C}_i` (*i.e.* the cell perimeter in 2D and the cell
            surface in 3D).

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
        # set the grid of distances {r}
        self._set_bounds(dr, bounds, distance_grid)
        # handle cutoff
        sides = self.trajectory.dump('cell.side')
        sides = numpy.asarray(sides)
        if cutoff is None:
            self.cutoff = numpy.min(sides) / 2
        else:
            self.cutoff = cutoff
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
        pos_all = self.trajectory.dump('position')
        spe_0 = self.dump('species_id', group=0)
        spe_all = self.trajectory.dump('species_id')
        Epot_0 = self.dump('particle.{}'.format(self.energy_field), group=0)
        Epot_all = self.trajectory.dump('particle.{}'.format(self.energy_field))
        cell_0 = self.dump('particle.{}'.format(self.cell_field), group=0)
        cell_all = self.trajectory.dump('particle.{}'.format(self.cell_field))
        box = self.trajectory.dump('cell.side')
        # compute extended neighbors with an extended cutoff
        n_pairs = len(self.trajectory[0].pairs_of_species)
        extended_cutoffs = numpy.array([self.cutoff for _ in range(n_pairs)])
        self._compute_extended_neighbors(extended_cutoffs)
        # the descriptor is computed for each species ID separately (1, 2, ...)
        # and also for all species at once, which is expressed as -1.
        distinct_species = list(self.trajectory[0].distinct_species_id)
        distinct_species += [-1]
        # computation
        start = 0
        n_L = len(self._distance_grid)
        for n in self._trange(n_frames):
            pos_0_n = pos_0[n].T
            pos_all_n = pos_all[n].T
            spe_0_n = spe_0[n],
            spe_all_n = spe_all[n]
            Epot_0_n = Epot_0[n]
            Epot_all_n = Epot_all[n]
            cell_0_n = cell_0[n]
            cell_all_n = cell_all[n]
            npart = len(self.groups[0][n])
            max_num_neigh_n = numpy.max(self._extended_neighbors_number[n])
            for spe_n, beta in enumerate(distinct_species):
                feat_spe_n = compute.physics_inspired_all(self._distance_grid,
                                                          self._extended_neighbors[n],
                                                          self._extended_neighbors_number[n],
                                                          max_num_neigh_n,
                                                          pos_0_n, pos_all_n,
                                                          spe_0_n, spe_all_n, beta,
                                                          Epot_0_n, Epot_all_n, 
                                                          cell_0_n, cell_all_n,
                                                          box[n])
                self.features[start: start+npart, 4*n_L*spe_n: 4*n_L*(spe_n+1) ] = feat_spe_n
            start += npart
        self._handle_nans()
        return self.features

    def _set_bounds(self, dr, bounds, distance_grid):
        if distance_grid is not None:
            self._distance_grid = numpy.array(distance_grid, dtype=numpy.float64)
            self._dr = None
            self._bounds = (self._distance_grid[0], self._distance_grid[-1])
        
        else:
            self._dr = dr
            if len(bounds) == 2 and bounds[0] >= 0 and bounds[1] >= 0 and bounds[0] < bounds[1]:
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
        observables = ['DEN', 'EPOT', 'CELL', 'VAR_EPOT']
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