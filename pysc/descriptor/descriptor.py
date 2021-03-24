from pysc.trajectory import Trajectory
from pysc.trajectory.particle import aliases
import numpy
import re

def _standardize_condition(condition):
    """
    Check that the condition is correctly formated (i.e <attr> _operator_ <val>).
    """
    regexp = re.search('(\w+)\s?(<|<=|==|>=|>)\s?([\'|\"]?\w+[\'|\"])', condition)
    if regexp:
        attr = regexp.group(1)
        operator = regexp.group(2)
        value = regexp.group(3)
        # if attribute is an alias, replace it with full attribute
        if not(attr.startswith('particle')):
            if attr in aliases:
                attr = aliases[attr]
            else:
                raise ValueError('attribute "{}" is not recognized'.format(attr))
        # reformat condition
        condition = '{} {} {}'.format(attr, operator, value)
        return condition
    else:
        raise ValueError('"{}" is not a valid condition'.format(condition))

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
    Examples:
    ---------
    - self.add_filter("species == 'A'", group=0)
    - self.add_filter("species == 'B'", group=1)
    It is then possible to access to the particles' properties of each group separately.
    The active filters are stored in the `active_filters` attribute.
    Active filters on a given group can be cleared using the method `clear_filters`.
    Active filters on both groups can be cleared using the method `clear_all_filters`. 
    """
    
    def __init__(self, trajectory):
        # Trajectory
        if isinstance(trajectory, Trajectory):
            self.trajectory = trajectory
        elif isinstance(trajectory, str):
            self.trajectory = Trajectory(trajectory)
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
        Returns the indices of the particles included in `group`.
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
        Returns the positions of the particles in `group`.
        Keeps the trajectory format (frames).
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
        Returns the species of the particles in `group`.
        Keeps the trajectory format (frames).
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
        Returns the species' ID of the particles in `group`.
        Keeps the trajectory format (frames).
        """
        _species_id = []
        for frame in self._groups[group]:
            _species_frame = numpy.empty(len(frame), dtype=numpy.int64)
            for n, particle in enumerate(frame):
                _species_frame[n] = particle._species_id
            _species_id.append(_species_frame)
        return _species_id
    
    def group_fraction(self, group):
        """
        Fraction of particles inside `group` over the whole trajectory.
        """
        N_group = self.group_size(group)
        N_tot = numpy.sum([sys.number_of_particles for sys in self.trajectory])
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

    #TODO: let the user choose the nearest-neighbor method
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
