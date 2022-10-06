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

    verbose : int, default: 0
        Show progress information about the computation of the descriptor
        when verbose is 1, and remain silent when verbose is 0 (default).
    
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

    def __init__(self, trajectory, accept_nans=True, verbose=0):
        # Trajectory
        # TODO: we can't change format or backend when passing a string
        if isinstance(trajectory, str):
            self.trajectory = Trajectory(trajectory)
        else:
            self.trajectory = trajectory
        # dismiss rows of self.features who contain nans
        self._accept_nans = accept_nans
        # verbose for computation
        self.verbose = verbose
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

    @property
    def accept_nans(self):
        """
        Boolean property.
        
        If False, discard any row from the array of features that contains a
        NaN element. Warning: this will delete the selected rows from the
        array of featyres. Use the method dismiss_nans instead to return a 
        filtered copy of the array of features.
        
        If True, keep NaN elements in the array of features.
        """
        return self._accept_nans

    @accept_nans.setter
    def accept_nans(self, accept):
        self._accept_nans = accept
        if not accept and self.features is not None:
            self._handle_nans()

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
            - 'mass': 'particle.mass'
            - 'radius': 'particle.radius'
            - 'nearest_neighbors': 'particle.nearest_neighbors'
            - 'neighbors': 'particle.nearest_neighbors'
            - 'neighbours': 'particle.nearest_neighbors'
            - 'voronoi_signature': 'particle.voronoi_signature'
            - 'signature': 'particle.voronoi_signature'

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
        >>> D.get_group_property('position', group=0)
        >>> D.get_group_property('x', group=1)
        >>> D.get_group_property('energy', group=0)
        
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
            to_dump.append(numpy.array(to_dump_frame, dtype=object))
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
        if self._accept_nans and self.verbose:
            isfinite = numpy.isfinite(self.features)
            collapsed_rows = numpy.product(isfinite, axis=1, dtype=bool)
            num_nans = self.size - numpy.sum(collapsed_rows)
            if num_nans > 0:
                print("Warning: {} NaN sample(s) in the array of features. This will compromise the computation of the average.".format(num_nans))
        return numpy.mean(self.features, axis=0)

    def _trange(self, bound):
        if self.verbose == 1:
            try:
                from tqdm import trange
                return trange(bound, desc=self.symbol)
            except ImportError:
                print('install tqdm to show the progress bar')
                return range(bound)
        return range(bound)

    def compute(self):
        pass

    def normalize(self, dist):
        """
        Generic normalization function for child classes. Returns the input
        distribution unchanged.
        """
        return dist

    def discard_nans(self):
        """
        Return the array of features where each row that contains
        NaN is filtered out.
        """
        isfinite = numpy.isfinite(self.features)
        collapsed_rows = numpy.product(isfinite, axis=1, dtype=bool)
        num_nans = self.size - numpy.sum(collapsed_rows)
        if num_nans > 0 and self.verbose:
            print('Warning: discarding {} NaN samples from the analysis.'.format(num_nans))
        return self.features[collapsed_rows]

    def _group_check(self):
        # check that groups are not empty
        for gn in range(2):
            if self.group_size(gn) == 0:
                raise AssertionError("group {} is empty. Check the filters on your descriptor.".format(gn))        

    def _set_up(self, dtype=numpy.int64):
        # initialize the data matrix
        self.features = numpy.empty((self.size, self.n_features), dtype=dtype)

    def _handle_nans(self):
        if not self._accept_nans:
            self.features = self.discard_nans()

class AngularStructuralDescriptor(StructuralDescriptor):
    """
    Base class for angular structural descriptors.
    
    See the parent class for more details.
    
    Descriptors that exploit angular correlations and require 
    neighbors information will inherit from this class.
    
    Parameters
    ----------
    
    trajectory : str or an instance of `Trajectory`.
        Trajectory on which the structural descriptor will be computed.
    """

    def __init__(self, trajectory, accept_nans=True, verbose=0):
        StructuralDescriptor.__init__(self, trajectory,
                                      accept_nans=accept_nans,
                                      verbose=verbose)

    def _manage_nearest_neighbors(self):
        """
        Check if nearest neighbors were already computed. If not,
        compute them using the method specified in the trajectory.
        """
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

    def _manage_nearest_neighbors_cutoffs(self):
        """
        Check if nearest neighbors cutoffs were already computed. If
        not, compute them from the trajectory.
        """
        if None in self.trajectory.nearest_neighbors_cutoffs:
            self.trajectory.compute_nearest_neighbors_cutoffs()

    def _filter_neighbors(self):
        """
        Create a list of neighbors separate from the trajectory for particles
        in group=0 by removing the neighbors that are not in group=1
        (e.g. for partial correlations).
        """
        # is there an active filter on group=1?
        needs_filtering = False
        for _, group in self.active_filters:
            if group == 1:
                needs_filtering = True
                break
        # if no filter, use unfiltered neighbors from the trajectory
        if not needs_filtering:
            self._neighbors = self.dump('nearest_neighbors', group=0)
        # else filter out nearest neighbors not in group=1
        else:
            n_frames = len(self.trajectory)
            idx_1 = self.dump('_index', group=1)
            self._neighbors = [[] for _ in range(n_frames)]
            for n in range(n_frames):
                idx_1_n = set(idx_1[n])
                for pi in self.groups[0][n]:
                    neigh_pi = set(pi.nearest_neighbors)
                    selected_neigh_pi = list(neigh_pi & idx_1_n)
                    self._neighbors[n].append(selected_neigh_pi)

    def _filter_subsidiary_neighbors(self):
        """
        Create a subsidiary list of neighbors for the neighbors of
        particles in group=0 by removing the neighbors that are not
        in group=1 (e.g. for partial correlations).
        """
        n_frames = len(self.trajectory)
        idx_1 = self.dump('_index', group=1)
        self._subsidiary_neighbors = [[] for n in range(n_frames)]
        for n, system in enumerate(self.trajectory):
            idx_1_n = set(idx_1[n])
            for neigh_pi in self._neighbors[n]:
                selected_neigh_neigh_pi = []
                for j in neigh_pi:
                    neigh_pj = set(system.particle[j].nearest_neighbors)
                    selected_neigh_pj = list(neigh_pj & idx_1_n)
                    selected_neigh_neigh_pi.append(selected_neigh_pj)
                self._subsidiary_neighbors[n].append(selected_neigh_neigh_pi)

    def _compute_extended_neighbors(self, cutoffs):
        """
        Compute the extended neighbors of particles in group=0 at given distances.
        This is different from the **nearest** neighbors in the trajectory for which
        the nearest neighbors are clearly defined on the basis of various methods.
        """
        n_frames = len(self.trajectory)
        # self._extended_neighbors = [[] for _ in range(n_frames)]
        # TODO: directly initialize lists with arrays at the right size
        self._extended_neighbors = []
        self._extended_neighbors_number = []
        #  indices
        idx_0 = self.dump('_index', group=0)
        idx_1 = self.dump('_index', group=1)
        #  species
        spe_0_id = self.dump('species_id', group=0)
        spe_1_id = self.dump('species_id', group=1)
        n_species = len(self.trajectory[0].distinct_species)
        #  positions
        pos_0 = self.dump('position', group=0)
        pos_1 = self.dump('position', group=1)
        #  box
        box = self.trajectory.dump('cell.side')
        #  cutoffs squared
        cutoffs_sq = numpy.array(cutoffs, dtype=numpy.float64)**2
        # TODO: why T?
        cutoffs_sq = cutoffs_sq.reshape(n_species, n_species).T
        for n in range(n_frames):
            n_0 = len(idx_0[n])
            # TODO: now hard coded 100?
            neighbors = numpy.zeros((n_0, 100), dtype=numpy.int64, order='F')
            num_neighbors = numpy.zeros(n_0, dtype=numpy.int64)
            pos_0_n = pos_0[n].T
            pos_1_n = pos_1[n].T
            nearest_neighbors_f90.fixed_cutoffs_distinct(idx_0[n], idx_1[n],
                                                         pos_0_n, pos_1_n,
                                                         spe_0_id[n], spe_1_id[n],
                                                         box[n], cutoffs_sq,
                                                         num_neighbors, neighbors)
            # TODO: add some assertion test here (unless there are checks in f90)
            self._extended_neighbors.append(neighbors)
            self._extended_neighbors_number.append(num_neighbors)
            # TODO: provide function accessor to neighbors in list form as below, though slow
            #for i, neigh_i in enumerate(neighbors):
            #    self._extended_neighbors[n].append(neigh_i[0:num_neighbors[i]])


class DummyDescriptor(StructuralDescriptor):

    name = 'dummy'
    symbol = 'dm'

    def __init__(self):
        self.grid = [0, 1]
        self.features = None

    def normalize(self, dist):
        return dist * (1.0 / numpy.sum(dist))
