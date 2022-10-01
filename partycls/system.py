"""
The physical system at hand.

The system of interest in a classical atomistic simulations is
composed of interacting point particles, usually enclosed in a
simulation cell.

This class is inspired by the `atooms` framework authored by Daniele Coslovich
See https://framagit.org/atooms/atooms 
"""

import re
import numpy
from .particle import aliases
from .core.utils import standardize_condition, NearestNeighborsMethod, _nearest_neighbors_methods_
from .neighbors_wrap import nearest_neighbors as nearest_neighbors_f90


class System:
    """
    A system is composed of a collection of particles that lie within an
    orthorhombic cell.
    
    Parameters
    ----------
    
    particle : list of `Particle`, optional, default: None
        A list of instances of `Particle`.
    
    cell : Cell, optional, default: None
        The cell (simulation box).
    
    Attributes
    ----------
        
    particle : list of `Particle`
        All the particles in the system.
        
    cell : Cell
        The cell where all the particles lie.
    
    Examples
    --------
    
    >>> p = [Particle(position=[0.0, 0.0, 0.0], species='A'),
             Particle(position=[1.0, 1.0, 1.0], species='B')]
    >>> c = Cell([5.0, 5.0, 5.0])
    >>> sys = System(particle=p, cell=c)
    """

    def __init__(self, particle=None, cell=None):
        if particle is None:
            particle = []
        self.particle = particle
        self.cell = cell
        # nearest neighbors
        self._nearest_neighbors_method = NearestNeighborsMethod.Auto
        self.nearest_neighbors_cutoffs = [None for pair in self.pairs_of_species]

    @property 
    def nearest_neighbors_method(self):
        return self._nearest_neighbors_method.value

    @nearest_neighbors_method.setter
    def nearest_neighbors_method(self, value):
        self._nearest_neighbors_method = NearestNeighborsMethod(value.lower())

    @property
    def n_dimensions(self):
        """
        Number of spatial dimensions, guessed from the length of
        `particle[0].position`.
        """
        if len(self.particle) > 0:
            return len(self.particle[0].position)
        else:
            return 0

    @property
    def density(self):
        """
        Number density of the system.

        It will raise a ValueException if `cell` is None.
        """
        if self.cell is None:
            return ValueError('cannot compute density without a cell')
        return len(self.particle) / self.cell.volume

    @property
    def distinct_species(self):
        """
        Sorted numpy array of all the distinct species in the system.
        """
        return numpy.array(sorted(set(self.get_property('species'))))

    @property
    def pairs_of_species(self):
        """
        List of all the possible pairs of species.
        """
        pairs = []
        for s1 in self.distinct_species:
            for s2 in self.distinct_species:
                pairs.append((s1, s2))
        return pairs

    @property
    def pairs_of_species_id(self):
        """
        List of all the possible pairs of species ID.
        """
        pairs = []
        for i in range(len(self.distinct_species)):
            for j in range(len(self.distinct_species)):
                pairs.append((i + 1, j + 1))
        return pairs

    @property
    def chemical_fractions(self):
        """
        Numpy array with the chemical fractions of each species in the system.
        """
        species = self.get_property('species')
        fractions = numpy.empty(len(self.distinct_species))
        for i, species_i in enumerate(self.distinct_species):
            fractions[i] = numpy.sum(species == species_i) / len(self.particle)
        return fractions

    def get_property(self, what, subset=None):
        """
        Return a numpy array with the system property specified by `what`.
        If `what` is a particle property, return the property for all particles
        in the system, or for a given subset of particles specified by `subset`.
        
        Parameters
        ----------
        what : str
            Requested system property.
        
            `what` must be of the form 
            "particle.<attribute>" or "cell.<attribute>". 
            
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
            
        subset : str, optional
            Subset of particles for which the property must be dumped. Must be 
            of the form "particle.<attribute>" unless "<attribute>" is an 
            alias. The default is None (all particles will be included).
            This is ignored if `what` is cell property.
            
        Returns
        -------
        to_dump : numpy.ndarray
            Array of the requested system property.
        
        Examples
        --------
        >>> traj = Trajectory('trajectory.xyz')
        >>> sys = traj[0]
        >>> pos_0 = sys.get_property('position')
        >>> spe_0 = sys.get_property('species')
        >>> sides = sys.get_property('cell.side')
        """
        if what in aliases:
            what = aliases[what]

        # Set the property to a given subset?
        if subset is not None:
            condition = standardize_condition(subset)
        else:
            condition = 'True'

        # Make array of the attribute
        attr = what.split('.')[-1]
        if what.startswith('particle'):
            data = []
            for particle in self.particle:
                if eval(condition):
                    data.append(eval('particle.{}'.format(attr)))
            data = numpy.array(data, dtype=object)
        elif what.startswith('cell'):
            what = what.split('.')[-1]
            regexp = re.search(r'(\w+)\[(\w+)\]', what)
            # cell iterable property
            if regexp:
                what = regexp.group(1)
                idx = int(regexp.group(2))
                data = getattr(self.cell, what)[idx]
            else:
                data = getattr(self.cell, attr)
        else:
            raise ValueError('Unknown attribute %s' % what)
        return data

    def dump(self, what, subset=None):
        """
        Alias for the method get_property.
        """
        return self.get_property(what, subset=subset)

    def set_property(self, what, value, subset=None):
        """
        Set a system property `what` to `value`. If `what` is a particle 
        property, set the property for all the particles in the system or for a 
        given subset of particles specified by `subset`.

        Parameters
        ----------
        what : str
            Name of the property to set. This is considered to be a particle
            property by default, unless it starts with "cell", e.g. 
            "cell.side".
            
        value : int, float, list, or numpy.ndarray
            Value(s) of the property to set. An instance of `int` or `float`
            will set the same value for all concerned particles. An instance
            of `list` or `numpy.ndarray` will assign a specific value to each
            particle. In this case, the size of `value` should respect the
            number of concerned particles.
            
        subset : str, optional
            Particles to which the property must be set. The default is None.
            This is ignored if `what` is cell property.

        Returns
        -------
        None.

        Examples
        --------
        >>> sys.set_property('mass', 1.0)
        >>> sys.set_property('radius', 0.5, "species == 'A'")
        >>> labels = [0, 1, 0] # 3 particles in the subset
        >>> sys.set_property('label', labels, "species == 'B'")
        >>> sys.set_property('cell.side[0]', 2.0)

        """

        # Set the property to a given subset?
        if subset is not None:
            condition = standardize_condition(subset)
        else:
            condition = 'True'

        # Set the same scalar value to each selected particle/cell
        if not isinstance(value, (list, numpy.ndarray)):
            if what.startswith('cell'):
                what = what.split('.')[-1]
                regexp = re.search('(\w+)\[(\w+)\]', what)
                # cell iterable property
                if regexp:
                    what = regexp.group(1)
                    idx = int(regexp.group(2))
                    getattr(self.cell, what)[idx] = value
                else:
                    setattr(self.cell, what, value)
            else:
                if what.startswith('particle'):
                    what = what.split('.')[-1]
                for particle in self.particle:
                    if eval(condition):
                        setattr(particle, what, value)

        # Set a specific value to each particle/cell with a list/array
        else:
            if what.startswith('cell'):
                what = what.split('.')[-1]
                setattr(self.cell, what, value)
            else:
                if what.startswith('particle'):
                    what = what.split('.')[-1]
                c = 0
                for particle in self.particle:
                    if eval(condition):
                        setattr(particle, what, value[c])
                        c += 1

    def compute_nearest_neighbors(self, method, cutoffs):
        """
        Compute the nearest neighbors for all the particles in the trajectory using
        the provided method.

        Parameters
        ----------
        method : str, default: None
            Method to identify the nearest neighbors. Must be one of 
            ['fixed', 'sann', 'voronoi].

        cutoffs : list
            List containing the cutoffs distances for each pair of species
            in the system (for method 'fixed' and 'sann').

        Returns
        -------
        None.
        """

        # Set up
        if isinstance(method, str):
            method = NearestNeighborsMethod(method.lower())
            self._nearest_neighbors_method = method
        else:
            self._nearest_neighbors_method = method
        positions = self.dump('position')
        species_id = self.dump('species_id')
        pairs_of_species_id = numpy.asarray(self.pairs_of_species_id)
        indices = self.dump('particle._index')
        box = self.dump('cell.side')
        self.nearest_neighbors_cutoffs = cutoffs
        
        # Computation
        #  Fixed-cutoffs ('fixed')
        if method is NearestNeighborsMethod.Fixed:
            positions = positions.T
            cutoffs_sq = numpy.array(cutoffs)**2
            for p in self.particle:
                neigh_i = nearest_neighbors_f90.fixed_cutoffs(p._index, indices,
                                                              p.position, positions,
                                                              p.species_id, species_id,
                                                              pairs_of_species_id, box,
                                                              cutoffs_sq)
                neigh_i = neigh_i[neigh_i >= 0]
                p.nearest_neighbors = list(neigh_i)
            return
        
        #  Solid-Angle Nearest Neighbors ('sann')
        if method is NearestNeighborsMethod.SANN:
            positions = positions.T
            rmax = 1.5 * numpy.max(cutoffs)
            for p in self.particle:
                neigh_i = nearest_neighbors_f90.sann(p.position, positions,
                                                     p._index, indices,
                                                     rmax, box)
                neigh_i = neigh_i[neigh_i >= 0]
                p.nearest_neighbors = list(neigh_i)
            return

        #  Voronoi neighbors ('voronoi')
        if method is NearestNeighborsMethod.Voronoi:
            try:
                import pyvoro

                if self.n_dimensions == 2:
                    raise NotImplementedError("The computation of Voronoi neighbors is currently not possible in dimension 2.")

                # parameters
                limits = [[-L/2, L/2] for L in self.cell.side]
                radii = self.dump('radius')
                periodic = self.cell.periodic
                # For efficiency, voro++ divides the box into a
                #  grid of cubic blocks. In order to achieve maximum
                #  performance, a block should contain 3-8 particles.
                #  To do so, one can compute the side of a cube that
                #  would contain 5.5 particles based on the number
                #  density to set the block size (i.e. dispersion).
                dispersion = (5.5 / self.density)**(1.0/3.0)
                # computation
                voronoi = pyvoro.compute_voronoi(positions,
                                                 limits, 
                                                 dispersion,
                                                 radii=radii,
                                                 periodic=periodic)
                # attribution
                for i, pi in enumerate(self.particle):
                    neigh_i = []
                    for face in voronoi[i]['faces']:
                        neigh_i.append(face['adjacent_cell'])
                    pi.nearest_neighbors = neigh_i

            except ModuleNotFoundError:
                raise ModuleNotFoundError('No `pyvoro` module found.')
            return

    def compute_voronoi_signatures(self):
        """
        Compute the Voronoi signatures of all the particles in the system
        using the radical Voronoi tessellation method.
        
        Particle radii must be set using the `set_property` method if the 
        original trajectory file does not contain such information.

        Creates a `voronoi_signature` property for the particles.
        """
        try:
            import pyvoro

            if self.n_dimensions == 2:
                raise NotImplementedError("The computation of Voronoi signatures is currently not possible in dimension 2.")

            # parameters
            positions = self.dump('position')
            limits = [[-L/2, L/2] for L in self.cell.side]
            radii = self.dump('radius')
            periodic = self.cell.periodic
            dispersion = (5.5 / self.density)**(1.0/3.0)
            # computation
            voronoi = pyvoro.compute_voronoi(positions,
                                             limits, 
                                             dispersion,
                                             radii=radii,
                                             periodic=periodic)
            # attribution
            for i, pi in enumerate(self.particle):
                faces = []
                for face in voronoi[i]['faces']:
                    faces.append(len(face['vertices']))
                    signature = [faces.count(i) for i in range(3, max(faces)+1)]
                    signature = '_'.join(map(str, signature))
                    pi.voronoi_signature = signature

        except ModuleNotFoundError:
            raise ModuleNotFoundError('No `pyvoro` module found.')


    def show(self, backend='matplotlib', color='species', **kwargs):
        """
        Show a snapshot of the system and color particles
        according to an arbitrary property, such as species, cluster label, 
        etc. Current visualization backends are 'matplotlib', 'ovito' and 
        '3dmol'.

        Parameters
        ----------
        backend : str, optional
            Name of the backend to use for visualization. 
            The default is 'matplotlib'.
        color : str, optional
            Name of the particle property to use as basis for coloring the 
            particles. This property must be defined for all the particles in the system.
            The default is 'species'.
        **kwargs : additional keyworded arguments (backend-dependent).

        Raises
        ------
        ValueError
            In case of unknown `backend`.

        Returns
        -------
        Figure or View (backend-dependent)
        
        Examples
        --------
        >>> sys.show(frame=0, color='label', backend='3dmol')
        >>> sys.show(frame=1, color='energy', backend='matplotlib', cmap='viridis')

        """
        from .helpers import show_matplotlib, show_ovito, show_3dmol
        if backend == 'matplotlib':
            _show = show_matplotlib
        elif backend == 'ovito':
            _show = show_ovito
        elif backend == '3dmol':
            _show = show_3dmol
        else:
            raise ValueError('unknown backend for visualization')
        return _show(self, color=color, **kwargs)

    def fold(self):
        """
        Fold the particle positions into the central cell.

        Returns
        -------
        None.

        """
        for p in self.particle:
            p.fold(self.cell)

    def __str__(self):
        rep = 'System(number_of_particles={}, species={}, chemical_fractions={}, cell={})'
        return rep.format(len(self.particle),
                          self.distinct_species,
                          self.chemical_fractions,
                          self.cell.side)

    def __repr__(self):
        return self.__str__()
