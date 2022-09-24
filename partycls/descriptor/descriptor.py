import numpy
from partycls.trajectory import Trajectory
from partycls.core.utils import standardize_condition
from .realspace_wrap import compute
from partycls.particle import aliases
from partycls.neighbors_wrap import nearest_neighbors as nearest_neighbors_f90

class StructuralDescriptor:
    """
    Base class for structural descriptors.
    
    The descriptor is calculated for the provided trajectory `trajectory`. This can be:
    - an object implementing the `Trajectory` interface ;
    - the path to a trajectory file in a format recognized by partyclsl ;
    
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
        
    groups : tuple
        Composition of the groups: groups[0] and groups[1] contain lists of all
        the `Particle` objects in groups 0 and 1 respectively. Each element of 
        the tuple is a list of `Particle` in `trajectory`, e.g. groups[0][0] is 
        the list of all the particles in the first frame of `trajectory` that 
        belong to group=0.
    
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
        self.groups = ([], [])
        self._group_init(0)
        self._group_init(1)
        # Active filters (none by default)
        self.active_filters = []
        # Dimension is guessed from the first frame of the trajectory
        self.dimension = self.trajectory[0].n_dimensions
        # Features
        #  default is None (when the object is created)
        #  correctly assigned when is method compute() is called
        self.grid = None
        self.features = None

    def __str__(self):
        rep = 'Descriptor(name="{}", dimension={}, filters={})'
        return rep.format(self.name, self.dimension, self.active_filters)

    def __repr__(self):
        return self.__str__()

    def _group_init(self, group):
        """
        Initialize the group `group` with all the particles by default.
        """
        self.groups[group].clear()
        for system in self.trajectory:
            frame = []
            for particle in system.particle:
                frame.append(particle)
            self.groups[group].append(frame.copy())

    def add_filter(self, condition, group=0):
        """
        Add a filter on the group (0 or 1) to select the subset of particles
        that respects the provided condition.

        Parameters
        ----------
        condition : str
            The condition should have the following format:
    
            <attribute> _operator_ <value>
            
            where:
            - <attribute> is a particle property (accepts aliases) ;
            - _operator_ is a logical operator (<, <=, ==, >=, >) ;
            - <value> is the corresponding value of <attribute> with the proper type ;
        
        group : int, optional
            Index of the group to which the filter must be applied.
            The default is 0.

        Returns
        -------
        None.
        
        Examples:
        ---------
        >>> S = StructuralDescriptor('trajectory.xyz')
        >>> S.add_filter("particle.radius < 0.5")
        >>> S.add_filter("species == 'A'", group=1)
        >>> S.add_filter("x < 0", group=0) # particles on the left side of the box
      
        """
        condition = standardize_condition(condition)
        self.active_filters.append((condition, group))
        # Iterate over frames
        for frame in self.groups[group]:
            to_remove = []
            # First find particles to remove from the current frame
            for particle in frame:
                if not(eval(condition)):
                    to_remove.append(particle)
            # Actually remove them
            for p_to_rem in to_remove:
                frame.remove(p_to_rem)
        self._group_check()

    def clear_filters(self, group=0):
        """
        Clear all active filters on `group`.
        All particles are included again in `group`.        

        Parameters
        ----------
        group : int, optional
            Index of the group on which to clear the filters. The default is 0.

        Returns
        -------
        None.

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

        Returns
        -------
        None.

        """
        self._group_init(0)
        self._group_init(1)
        self.active_filters = []

    def group_size(self, group):
        """
        Return the number of particles in `group`.
        """
        N = 0
        for frame in self.groups[group]:
            N += len(frame)
        return N

    def get_group_property(self, what, group):
        """
        Return a list of numpy arrays with the properties of the particles in
        group `group`. The list size is the number of systems in the 
        trajectory.

        Parameters
        ----------
        what : str
            Requested particle property. 
            
            `what` must be a particle property or an alias.
            
            The following particle aliases are accepted:
            - 'position': 'particle.position'
            - 'pos': 'particle.position'
            - 'position[0]': 'particle.position[0]', 
            - 'pos[0]': 'particle.position[0]'
            - 'x': 'particle.position[0]'
            - 'position[1]': 'particle.position[1]',
            - 'pos[1]': 'particle.position[1]'
            - 'y': 'particle.position[1]'
            - 'position[2]': 'particle.position[2]'
            - 'pos[2]': 'particle.position[2]'
            - 'z': 'particle.position[2]'
            - 'species': 'particle.species'
            - 'spe': 'particle.species'
            - 'label': 'particle.label'
            - 'index': 'particle.index'
            - 'mass': 'particle.mass'
            - 'radius': 'particle.radius'

        Returns
        -------
        to_dump : list
            List of the requested particle property with length equal to the 
            number of frames in the trajectory. Each element of the list is a
            numpy.ndarray of the requested particle property.
            
        Examples
        --------
        >>> traj = Trajectory('trajectory.xyz')
        >>> D = StructuralDescriptor(traj)
        >>> D.get_group_property('position', 0)
        >>> D.get_group_property('x', 1)
        >>> D.get_group_property('energy', 0)
        
        """
        if what in aliases:
            what = aliases[what]
        if what.startswith('particle'):
            what = what.split('.')[-1]
        to_dump = []
        for frame in self.groups[group]:
            to_dump_frame = []
            for particle in frame:
                to_dump_frame.append(eval('particle.{}'.format(what)))
            to_dump.append(numpy.array(to_dump_frame))
        return to_dump

    def dump(self, what, group):
        """
        Alias for the method get_group_property.
        """
        return self.get_group_property(what, group)

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
        return sum([len(frame) for frame in self.groups[0]])

    @property
    def n_features(self):
        """
        Number of features of the descriptor.
        """
        return len(self.grid)

    @property
    def average(self):
        """
        Average feature vector of the descriptor.
        """
        return numpy.mean(self.features, axis=0)

    def compute(self):
        pass

    def normalize(self, dist):
        """
        Generic normalization function for child classes. Returns the input
        distribution unchanged.
        """
        return dist

    def _group_check(self):
        # check that groups are not empty
        for gn in range(2):
            if self.group_size(gn) == 0:
                raise AssertionError("group {} is empty. Check the filters on your descriptor.".format(gn))        

    def _set_up(self, dtype=numpy.int64):
        # initialize the data matrix
        self.features = numpy.empty((self.size, self.n_features), dtype=dtype)


class AngularStructuralDescriptor(StructuralDescriptor):
    """
    Base class for angular structural descriptors.
    
    See the parent class for more details.
    
    Descriptors that exploit angular correlations and require nearest-neighbors information 
    will inherit of this class. Two methods to identify nearest-neighbors are available:
    - "Fixed cutoff" (symbol: 'FC'): uses the partial radial distribution functions to compute
      the cutoffs between each possible pair of species (s1, s2) ;
    - "Solid-Angle based Nearest Neighbors" (symbol: 'SANN'): see van Meel et al. (https://doi.org/10.1063/1.4729313)
    
    Neighbors can also be read from the trajectory file directly. Otherwise, 
    they can be computed by changing the attribute `nearest_neighbors_method`
    to 'FC' or 'SANN'.
    
    When using the 'FC' method, it is also possible to specify the cutoffs 
    manually for a pair of species (s1, s2) by using the method `set_cutoff`. 
    The cutoffs that were not set manually will be computed automatically.
    
    Parameters
    ----------
    
    trajectory : str or an instance of `Trajectory`.
        Trajectory on which the structural descriptor will be computed.
    
    Attributes
    ----------
    
    cutoffs : list of float
        List of cutoff distances to identify the nearest neighbors using
        the fixed-cutoff ('FC') method.
        
    nearest_neighbors_method : str, default: 'auto'
        Nearest neighbor method, 'FC' or 'SANN'. If method is 'auto', neighbors
        are read directly from the trajectory (if provided). If no neighbors 
        are found, it uses method='FC' instead.
        
    neighbors : list
        Lists of nearest neighbors for all the particles in group=0. None by
        default and filled when calling the method `nearest_neighbors`.
    """

    def __init__(self, trajectory):
        StructuralDescriptor.__init__(self, trajectory)
        self.cutoffs = [None for n in range(len(self.trajectory[0].pairs_of_species))]
        # 'FC' = Fixed Cutoff (default)
        # 'SANN' = Solid Angle Nearest Neighbors
        self.nearest_neighbors_method = 'auto'
        self.neighbors = None

    def set_cutoff(self, s1, s2, rcut, mirror=True):
        """
        Set the nearest-neighbor cutoff for the pair of species (s1, s2).
        The cutoff of the mirror pair (s2, s1) is set automatically if the `mirror` 
        parameter is True (default).        

        Parameters
        ----------
        s1 : str
            Symbol of the first species.
        s2 : str
            Symbol of the second species.
        rcut : float
            Value of the cutoff for the pair (s1,s2).
        mirror : bool, optional
            Set the cutoff for the mirror pair (s2,s1). The default is True.

        Returns
        -------
        None.
        """
        pairs = self.trajectory[0].pairs_of_species
        idx_12 = pairs.index((s1, s2))
        self.cutoffs[idx_12] = rcut
        if mirror:
            idx_21 = pairs.index((s2, s1))
            self.cutoffs[idx_21] = rcut

    def _filter_neighbors(self):
        """
        Create a separate list of neighbors from the neighbors in the trajectory
        by removing those that are not in group=1 for partial correlations.
        """
        n_frames = len(self.trajectory)
        self._neighbors = [[] for n in range(n_frames)]
        for frame, system in enumerate(self.trajectory):
            for pi in self.groups[0][frame]:
                neigh_pi = pi.nearest_neighbors
                selected_neigh_pi = []
                for j in neigh_pi:
                    pj = system.particle[j]
                    if pj in self.groups[1][frame]:
                        selected_neigh_pi.append(j)    
                self._neighbors[frame].append(selected_neigh_pi)

    def _filter_subsidiary_neighbors(self):
        n_frames = len(self.trajectory)
        self._subsidiary_neighbors = [[] for n in range(n_frames)]
        for frame, system in enumerate(self.trajectory):
            for neigh_pi in self._neighbors[frame]:
                selected_neigh_neigh_pi = []
                for j in neigh_pi:
                    pj = system.particle[j]
                    neigh_pj = pj.nearest_neighbors
                    selected_neigh_pj = []
                    for k in neigh_pj:
                        pk = system.particle[k]
                        if pk in self.groups[1][frame]:
                            selected_neigh_pj.append(k)
                    selected_neigh_neigh_pi.append(selected_neigh_pj)
                self._subsidiary_neighbors[frame].append(selected_neigh_neigh_pi)

    def _manage_nearest_neighbors(self):
        # check if nearest neighbors were already computed
        neighbors_None = False
        for system in self.trajectory:
            nearest_neighbors = system.dump('nearest_neighbors')
            if None in nearest_neighbors:
                neighbors_None = True
                break
        # if not, compute them using the method specified in the trajectory
        if neighbors_None:
            self.trajectory.compute_nearest_neighbors()

    def _compute_neighbors_fixed_cutoffs(self, cutoffs):
        """
        Compute the neighbors of particles in group=0 at given distances.
        This is different from nearest neighbors in the trajectory.
        """
        n_frames = len(self.trajectory)
        self._generalized_neighbors = [[] for n in range(n_frames)]
        #  indices
        idx_0 = self.dump('index', group=0)
        idx_1 = self.dump('index', group=1)
        #  species
        spe_0_id = self.dump('species_id', group=0)
        spe_1_id = self.dump('species_id', group=1)
        pairs_of_species_id = numpy.asarray(self.trajectory[0].pairs_of_species_id)
        #  positions
        pos_0 = self.dump('position', group=0)
        pos_1 = self.dump('position', group=1)
        #  box
        box = self.trajectory.dump('cell.side')
        for frame in range(n_frames):
            for i in range(len(idx_0[frame])):
                neigh_i = nearest_neighbors_f90.fixed_cutoffs(idx_0[frame][i], idx_1[frame],
                                                              pos_0[frame][i], pos_1[frame].T,
                                                              spe_0_id[frame][i], spe_1_id[frame],
                                                              pairs_of_species_id, box[frame],
                                                              cutoffs)
                neigh_i = neigh_i[neigh_i >= 0]
                self._generalized_neighbors[frame].append(neigh_i)

    def nearest_neighbors(self, method='auto'):
        """
        Compute the nearest neighbors of particles in group=0 using one of the
        following methods:
        - "Fixed cutoff" (method='FC'): uses the partial radial distribution functions 
          to compute the cutoffs between each possible pair of species (s1, s2) ;
        - "Solid-Angle based Nearest Neighbors" (method='SANN'): see  
           van Meel et al. (https://doi.org/10.1063/1.4729313) ;

        Parameters
        ----------
        method : str, default: 'auto'
            Method to identify nearest neighbors. Must be 'FC' or 'SANN'. If
            method is 'auto', neighbors are read directly from the trajectory
            (if provided). If no neighbors are found, it uses method='FC' 
            instead.

        Returns
        -------
        None.
        """
        #  number of frames
        n_frames = len(self.groups[0])
        #  list of neighbors
        self.neighbors = [[] for n in range(n_frames)]
        
        # Read neighbors from the trajectory
        if method == 'auto':
            for n in range(n_frames):
                self.neighbors[n] = numpy.array(self.dump('neighbors', group=0)[n])
            # if no neighbors in the trajectory, fallback to method='FC'
            if None in self.neighbors[0]:
                self.nearest_neighbors_method = 'FC'
                self.nearest_neighbors(method=self.nearest_neighbors_method)
            return
                
        # Compute neighbors
        #  indices
        idx_0 = self.dump('index', 0)
        idx_1 = self.dump('index', 1)
        idx_all = self.trajectory.get_property('index')
        #  species
        spe_0 = self.dump('species_id', 0)
        spe_1 = self.dump('species_id', 1)
        pairs = numpy.asarray(self.trajectory[0].pairs_of_species_id)
        #  positions
        pos_0 = self.dump('position', 0)
        pos_1 = self.dump('position', 1)
        pos_all = self.trajectory.get_property('position')
        # compute all/missing cutoffs
        if None in self.cutoffs:
            self._compute_cutoffs()
        cutoffs = numpy.array(self.cutoffs)
        #  boundaries
        box = self.trajectory[0].cell.side

        # Fixed cutoff ('FC')
        if method == 'FC':
            for n in range(n_frames):
                for i in range(len(idx_0[n])):
                    neigh_i = compute.nearest_neighbors(idx_0[n][i], idx_1[n],
                                                        pos_0[n][i], pos_1[n].T,
                                                        spe_0[n][i], spe_1[n],
                                                        pairs, box, cutoffs)
                    neigh_i = neigh_i[neigh_i >= 0]
                    self.neighbors[n].append(neigh_i)
            return

        #  Solid Angle Nearest Neighbors ('SANN')
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
            return

    # TODO: if fixed-cutoff method, let the user choose `dr`
    def _compute_cutoffs(self):
        from .gr import RadialDescriptor
        pairs = self.trajectory[0].pairs_of_species
        for pair in pairs:
            if self.cutoffs[pairs.index(pair)] is None:
                s1, s2 = pair
                # use the smallest side of the smallest box in case of
                #  non-constant volume trajectory
                sides = numpy.array(self.trajectory.get_property('cell.side'))
                L = numpy.min(sides)
                bounds = (0.0, L / 2)
                descriptor = RadialDescriptor(self.trajectory, dr=0.1, bounds=bounds)
                descriptor.add_filter("species == '{}'".format(s1), group=0)
                descriptor.add_filter("species == '{}'".format(s2), group=1)
                descriptor.compute()
                # grid and average descriptor
                r = descriptor.grid
                h_12 = descriptor.average
                # normalized g(r)
                g_12 = descriptor.normalize(h_12, method="gr")
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

    def normalize(self, dist):
        return dist * (1.0 / numpy.sum(dist))
