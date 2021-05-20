import numpy
from pysc.trajectory import Trajectory
from pysc.core.utils import _standardize_condition
from .realspace_wrap import compute

class StructuralDescriptor:    
    """
    Base class for structural descriptors.
    
    The descriptor is calculated for the provided trajectory `trajectory`. This can be:
    - an object implementing the `Trajectory` interface ;
    - the path to a trajectory file in a format recognized by pyscl ;
    
    A structural descriptor S(x) is a collection of N individual empirical correlation 
    functions {s_i(x)} at the particle level, defined over a grid {x_j} of M features.
    These are stored in the `features` array as a matrix usually refered to as the "data set":
    
    s_0(x_0) s_0(x_1) ... s_0(x_M)
    s_1(x_0) s_1(x_1) ... s_1(x_M)
    ...      ...          ...
    s_N(x_0) s_N(x_1) ... s_N(x_M)
    
    The `features` array is None by default and is computed only when the `compute()` method is called.
    
    The correlations can be calculated between two arbitrary subsets of particles called "groups":
    - group 0 is the main group, i.e. particles for which the correlations are being calculated ;
    - group 1 is the secondary group, i.e. particles that are being considered when calculating the correlations ;
    These groups are formed by adding filters on particles' properties (species, radius, position, etc.).
    
    Parameters
    ----------
    
    trajectory : str or an instance of `Trajectory`.
        Trajectory on which the structural descriptor will be computed.
    
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
    
    >>> D = StructuralDescriptor('trajectory.xyz')
    >>> D.add_filter("species == 'A'", group=0)
    >>> D.add_filter("species == 'B'", group=1)
    >>> D.active_filters
    [("particle.species == 'A'", 0), ("particle.species == 'B'", 1)]
    >>> D.clear_filters(0)
    >>> D.active_filters
    [("particle.species == 'B'", 1)]
    """
    
    def __init__(self, trajectory):
        # Trajectory
        # TODO: we can't change format or backend when passing a string
        if isinstance(trajectory, str):
            self.trajectory = Trajectory(trajectory)
        else:
            self.trajectory = trajectory
        # Default: consider all particles for the correlation
        self._groups = ([], [])
        self._group_init(0)
        self._group_init(1)
        # Active filters (none by default)
        self.active_filters = []
        # Dimension is guessed from the first frame of the trajectory
        self.dimension = self.trajectory[0].number_of_dimensions
        # Features
        #  default is None (when the object is created)
        #  correctly assigned when is method compute() is called
        self.grid = None
        self.features = None

    def _group_init(self, group):
        """
        Initialize the group `group` with all the particles by default.
        """
        self._groups[group].clear()
        for system in self.trajectory:
            frame = []
            for particle in system.particle:
                frame.append(particle)
            self._groups[group].append(frame.copy())
    
    def add_filter(self, condition, group=0):
        """
        Add a filter on the group (0 or 1) to select the subset of particles
        that respects the provided condition.
        
        `condition` should have the following format:
            
        <attribute> _operator_ <value>
        
        where:
        - <attribute> is a particle property (accepts aliases) ;
        - _operator_ is a logical operator (<, <=, ==, >=, >) ;
        - <value> is the corresponding value of <attribute> with the proper type ;
        
        Examples:
        ---------
        - "particle.radius < 0.5", "radius < 0.5" and "rad < 0.5" are all valid conditions ;
        - "particle.species == 'A'", "species == 'A'" and "spe == 'A'" are all valid conditions ;        
        """
        condition = _standardize_condition(condition)
        self.active_filters.append((condition, group))
        # Iterate over frames
        for frame in self._groups[group]:
            to_remove = []
            # First find particles to remove from the current frame
            for particle in frame:
                if not(eval(condition)):
                    to_remove.append(particle)
            # Actually remove them
            for p_to_rem in to_remove:
                frame.remove(p_to_rem)
                
    def clear_filters(self, group=0):
        """
        Clear all active filters on `group`.
        All particles are included again in `group`.
        """
        # Reset `group` with all the particles
        self._group_init(group)
        # Update the `active_filters` list
        for fltr in self.active_filters:
            if fltr[1] == group:
                self.active_filters.remove(fltr)
                
    def clear_all_filters(self):
        """
        Clear all active filters in both groups.
        All particles are included in both groups again.
        """
        self._group_init(0)
        self._group_init(1)
        self.active_filters = []
        
    def group_size(self, group):
        """
        Return the number of particles in `group`.
        """
        N = 0
        for frame in self._groups[group]:
            N += len(frame)
        return N
            
    def group_indices(self, group):
        """
        Return the indices of the particles included in `group`.
        Keeps the trajectory format (frames).
        """
        group_idx = []
        for frame in self._groups[group]:
            frame_indices = []
            for particle in frame:
                frame_indices.append(particle.index)
            group_idx.append(frame_indices)
        return group_idx
    
    def group_positions(self, group):
        """
        Return the positions of the particles in `group`.
        Keeps the trajectory format (i.e. list structure).
        """
        _pos = []
        for frame in self._groups[group]:
            _pos_frame = numpy.empty((len(frame), self.dimension), dtype=numpy.float64)
            for n, particle in enumerate(frame):
                _pos_frame[n] = particle.position
            _pos.append(_pos_frame)
        return _pos

    def group_species(self, group):
        """
        Return the species of the particles in `group`.
        Keeps the trajectory format (i.e. list structure).
        """
        _species = []
        dtype_species = type(self.trajectory[0].particle[0].species)
        for frame in self._groups[group]:
            _species_frame = numpy.empty(len(frame), dtype=dtype_species)
            for n, particle in enumerate(frame):
                _species_frame[n] = particle.species
            _species.append(_species_frame)
        return _species      
    
    def group_species_id(self, group):
        """
        Return the species' ID of the particles in `group`.
        Keeps the trajectory format (i.e. list structure).
        """
        _species_id = []
        for frame in self._groups[group]:
            _species_frame = numpy.empty(len(frame), dtype=numpy.int64)
            for n, particle in enumerate(frame):
                _species_frame[n] = particle.species_id
            _species_id.append(_species_frame)
        return _species_id
    
    def group_fraction(self, group):
        """
        Return the fraction of particles inside `group` over the whole trajectory.
        """
        N_group = self.group_size(group)
        N_tot = numpy.sum([len(sys.particle) for sys in self.trajectory])
        return N_group / N_tot
    
    @property
    def size(self):
        """
        Total number of particles in the descriptor (i.e. in group=0).
        """
        return sum([len(frame) for frame in self._groups[0]])
    
    @property
    def n_features(self):
        """
        Number of features of the descriptor.
        """
        pass
    
    @property
    def average(self):
        """
        Average feature vector of the descriptor.
        """
        return numpy.mean(self.features, axis=0)
    
    def compute(self):
        pass
        
    def sanity_checks(self):
        assert (self.group_size(0) > 0 and self.group_size(1) > 0), 'groups cannot be empty.'

class AngularStructuralDescriptor(StructuralDescriptor):
    """
    Base class for angular structural descriptors.
    
    Descriptors that exploit angular correlations and require nearest-neighbors information 
    will inherit of this class. Two methods to identify nearest-neighbors are available:
    - "Fixed cutoff" (symbol: 'FC'): uses the partial radial distribution functions to compute
      the cutoffs between each possible pair of species (s1, s2) ;
    - "Solid-Angle based Nearest Neighbors" (symbol: 'SANN'): see van Meel et al. (https://doi.org/10.1063/1.4729313)
    
    The nearest-neighbors method can be changed by modifying the attribute `nearest_neighbors_method`
    to 'FC' (default) or 'SANN'.
    
    When using the 'FC' method, it is also possible to specify the cutoffs manually
    for a pair of species (s1, s2) by using the method `set_cutoff`. The cutoffs
    that were not set manually will be computed automatically.
    
    Parameters
    ----------
    
    trajectory : str or an instance of `Trajectory`.
        Trajectory on which the structural descriptor will be computed.
    
    Attributes
    ----------
    
    cutoffs : list of float
        List of cutoff distances to identify the nearest neighbors using
        the fixed-cutoff ('FC') method.
        
    nearest_neighbors_method : str, default: 'FC'
        Nearest neighbor method, 'FC' or 'SANN'.
    
    Examples:
    ---------
    
    """
    
    def __init__(self, trajectory):
        StructuralDescriptor.__init__(self, trajectory)
        self.cutoffs = [None for n in range(len(self.trajectory[0].pairs_of_species))]
        # 'FC' = Fixed Cutoff (default)
        # 'SANN' = Solid Angle Nearest Neighbors
        self.nearest_neighbors_method = 'FC'
        
    def set_cutoff(self, s1, s2, rcut, mirror=True):
        """
        Set the nearest-neighbor cutoff for the pair of species (s1, s2).
        The cutoff of the mirror pair (s2, s1) is set automatically if the `mirror` 
        parameter is True (default).
        """
        pairs = self.trajectory[0].pairs_of_species
        idx_12 = pairs.index((s1, s2))
        self.cutoffs[idx_12] = rcut
        if mirror:
            idx_21 = pairs.index((s2, s1))
            self.cutoffs[idx_21] = rcut    

    #TODO: define self.neighbors as an attribute for the class
    def nearest_neighbors(self, method='FC'):
        """
        Compute the nearest neighbors of particles in group=0 using one of the
        following methods:
        - "Fixed cutoff" (method='FC'): uses the partial radial distribution functions 
          to compute the cutoffs between each possible pair of species (s1, s2) ;
        - "Solid-Angle based Nearest Neighbors" (method='SANN'): see  
           van Meel et al. (https://doi.org/10.1063/1.4729313) ;
        """
        # indices
        idx_0, idx_1 = self.group_indices(0), self.group_indices(1)
        idx_all = self.trajectory.dump('index')
        # species
        spe_0, spe_1 = self.group_species_id(0), self.group_species_id(1)
        pairs = numpy.asarray(self.trajectory[0].pairs_of_species_id)
        # positions
        pos_0, pos_1 = self.group_positions(0), self.group_positions(1)
        pos_all = self.trajectory.dump('position')
        # compute all/missing cutoffs
        if None in self.cutoffs: self._compute_cutoffs()
        cutoffs = numpy.array(self.cutoffs)
        # boundaries
        n_frames = len(self._groups[0])
        box = self.trajectory[0].cell.side
        # list of neighbors
        self.neighbors = [[] for n in range(n_frames)]
        
        # Fixed cutoff
        if method == 'FC':
            for n in range(n_frames):
                for i in range(len(idx_0[n])):
                        neigh_i = compute.nearest_neighbors(idx_0[n][i], idx_1[n],
                                                            pos_0[n][i], pos_1[n].T,
                                                            spe_0[n][i], spe_1[n],
                                                            pairs, box, cutoffs)
                        neigh_i = neigh_i[neigh_i >= 0]
                        self.neighbors[n].append(neigh_i)

        #  Solid Angle Nearest Neighbors (SANN)
        #   This will find all neighbors of `i`
        #   (including particles not in group=1)            
        if method == 'SANN':
            # scaling factor for first guess as trying neighbors
            rmax = 1.5 * numpy.max(cutoffs)
            for n in range(n_frames):
                for i in range(len(idx_0[n])):
                    neigh_i = compute.sann(pos_0[n][i], pos_all[n].T,
                                           idx_0[n][i], idx_all[n], idx_1[n],
                                           rmax, box)
                    neigh_i = neigh_i[neigh_i >= 0]
                    self.neighbors[n].append(neigh_i)

    #TODO: if fixed-cutoff method, let the user choose `dr`
    def _compute_cutoffs(self):
        from .gr import RadialDescriptor
        pairs = self.trajectory[0].pairs_of_species
        for pair in pairs:
            if self.cutoffs[pairs.index(pair)] is None:
                s1, s2 = pair
                rlim = (0.0, numpy.min(self.trajectory[0].cell.side/2))
                descriptor = RadialDescriptor(self.trajectory, dr=0.1, rlim=rlim)
                descriptor.add_filter("species == '{}'".format(s1), group=0)
                descriptor.add_filter("species == '{}'".format(s2), group=1)
                descriptor.compute()
                # grid and average descriptor
                r = descriptor.grid
                h_12 = descriptor.average
                # normalized g(r)
                g_12 = descriptor.normalize_gr(h_12)
                # find the first minimum of g_12(r)
                first_max = numpy.argmax(g_12)
                first_min = numpy.argmin(g_12[first_max:]) + first_max
                rcut = r[first_min]
                # set the cutoff
                self.set_cutoff(s1, s2, rcut)
    
class DummyDescriptor(StructuralDescriptor):
    
    name = 'dummy'
    symbol = 'dm'
    
    def __init__(self):
        self.grid = [0, 1]
        self.features = None

    @property
    def size(self):
        return self.features.shape[0]
        
    @property
    def n_features(self):
        return self.features.shape[1]
    
    def normalize(self, dist):
        return dist * (1.0 / numpy.sum(dist))
