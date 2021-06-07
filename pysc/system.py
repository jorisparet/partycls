"""
This class is inspired by the `atooms` framework authored by Daniele Coslovich
See https://framagit.org/atooms/atooms 
"""

import numpy
from .particle import Particle, aliases
from .cell import aliases as cell_aliases
from pysc.core.utils import standardize_condition

# combine aliases
aliases.update(cell_aliases)

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
    >>> c = Cell(side = [5.0, 5.0, 5.0])
    >>> sys = System(particle=p, cell=c)
    """
    
    def __init__(self, particle=None, cell=None):
        if particle is None:
            particle = []
        self.particle = particle
        self.cell = cell
    
    @property
    def number_of_dimensions(self):
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
        return numpy.array(sorted(set(self.dump('species'))))
    
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
                pairs.append((i+1, j+1))
        return pairs
    
    @property
    def chemical_fractions(self):
        """
        Numpy array with the chemical fractions of each species in the system.
        """
        species = self.dump('species')
        fractions = numpy.empty(len(self.distinct_species))
        for i, species_i in enumerate(self.distinct_species):
            fractions[i] = numpy.sum(species == species_i) / len(self.particle)
        return fractions

    def dump(self, what):
        """
        Return a numpy array with the system property specified by `what`.
        
        Parameters
        ----------
        what : str
            Requested system property.
        
            `what` must be of the form 
            "particle.<attribute>" or "cell.<attribute>". 
            
            The following aliases are allowed:
            - "pos" ("particle.position")
            - "position" ("particle.position")
            - "position_x" ("particle.position[0]")
            - "x" ("particle.position[0]")
            - "position_y" ("particle.position[1]")
            - "y" ("particle.position[1]")
            - "position_z" ("particle.position[2]")
            - "z" ("particle.position[2]")
            - "spe" ("particle.species")
            - "species" ("particle.species")
            - "species_id" ("particle.species_id")
            - "radius" ("particle.radius")
            - "label" ("particle.label")
            - "index" ("particle.index")
            - "mass" ("particle.mass")
            - "box" ("cell.side")
        
        Returns
        -------
        to_dump : numpy.ndarray
            Array of the requested system property.
        
        Examples
        --------
        >>> traj = Trajectory('trajectory.xyz')
        >>> sys = traj[0]
        >>> pos_0 = sys.dump('position')
        >>> spe_0 = sys.dump('species')
        """
        if what in aliases:
            what = aliases[what]
        
        # Make array of the attribute
        attr = what.split('.')[-1]
        if what.startswith('particle'):
            data = numpy.array([eval('p.{}'.format(attr)) for p in self.particle])
        elif what.startswith('cell'):
             data = numpy.array(getattr(self.cell, attr))
        else:
            raise ValueError('Unknown attribute %s' % what)
        return data
    
    def set_property(self, what, value, subset=None):
        """
        Set a property `what` to `value` for all the particles in the 
        system or for a given subset of particles specified by `subset`.
        

        Parameters
        ----------
        what : str
            Name of the property to set.
        value : int, float, list, or numpy.ndarray
            Value(s) of the property to set. An instance of `int` or `float`
            will set the same value for all concerned particles. An instance
            of `list` or `numpy.ndarray` will assign a specific value to each
            particle. In this case, the size of `value` should respect the
            number of concerned particles.
        subset : str, optional
            Particles to which the property must be set. The default is None.

        Returns
        -------
        None.

        Examples
        --------
        >>> sys.set_property('mass', 1.0)
        >>> sys.set_property('radius', 0.5, "species == 'A'")
        >>> labels = [0, 1, 0] # 3 particles in the subset
        >>> sys.set_property('label', labels, "species == 'B'")

        """
        
        # Set the property to a given subset?
        if subset is not None:
            condition = standardize_condition(subset)
        else:
            condition = 'True'
        # Set the same scalar value to each selected particle
        if not isinstance(value, (list, numpy.ndarray)):
            for particle in self.particle:
                if eval(condition):
                    setattr(particle, what, value)
        # Set a specific value to each particle with a list/array
        else:
            c = 0
            for particle in self.particle:
                if eval(condition):
                    setattr(particle, what, value[c])
                    c += 1

    def show(self, backend='matplotlib', color='species', *args, **kwargs):
        """
        Show a snapshot of the system and color particles
        according to an arbitrary property, such as species, cluster label, 
        etc. Current visualization backends are 'matplotlib' and '3dmol'.

        Parameters
        ----------
        backend : str, optional
            Name of the backend to use for visualization. 
            The default is 'matplotlib'.
        color : str, optional
            Name of the particle property to use as basis for coloring the 
            particles. This property must be defined for all the particles in the system.
            The default is 'species'.
        *args : additional non-keyworded arguments (backend-dependent).
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
        from .helpers import show_matplotlib, show_3dmol
        if backend == 'matplotlib':
            _show = show_matplotlib
        elif backend == '3dmol':
            _show = show_3dmol
        else:
            raise ValueError('unknown backend for visualization')
        return _show(self, color=color, *args, **kwargs)   

    def __str__(self):
        rep = 'System(number_of_particles={}, species={}, chemical_fractions={}, cell={})'
        return rep.format(len(self.particle),
                          self.distinct_species,
                          self.chemical_fractions,
                          self.cell.side)
    
    def __repr__(self):
        return self.__str__()