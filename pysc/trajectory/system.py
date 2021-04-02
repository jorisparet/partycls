# Reference to atooms

from .particle import Particle, aliases
import numpy

class System:
    """
    A system is composed of a collection of particles that lie within a
    rectangular cell.
    
    Parameters
    ----------
    
    particle : list of `Particle`, optional, default: None
        A list of instances of `Particle`.
    
    cell : Cell, optional, default: None
        The rectangular cell (simulation box).
    
    Attributes
    ----------
        
    particle : list of `Particle`
        All the particles in the system.
        
    cell : Cell
        The rectangular cell where all the particles lie.
    
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
        
    def add_particle(self, particle):
        """
        Add a `Particle` to the `self.particle` list.
        """
        if isinstance(particle, Particle):
            self.particle.append(particle)
        else:
            raise TypeError('can only add an instance of `Particle`')
    
    @property
    def number_of_dimensions(self):
        """
        Return the number of spatial dimensions, guessed from the length of
        `particle[0].position`.
        """
        if len(self.particle) > 0:
            return len(self.particle[0].position)
        else:
            return 0
        
    @property
    def number_of_particles(self):
        """
        Return the total number of particles in the system.
        """
        return len(self.particle)
        
    @property
    def density(self):
        """
        Return the number density of the system.

        It will raise a ValueException if `cell` is None.
        """
        if self.cell is None:
            return ValueError('cannot compute density without a cell')
        return len(self.particle) / self.cell.volume
    
    @property
    def distinct_species(self):
        """
        Return a sorted numpy array of all the distinct species in the system.
        """
        return numpy.array(sorted(set(self.dump('species'))))
    
    @property
    def number_of_species(self):
        """
        Return the number of distinct species in the system.
        """
        return len(self.distinct_species)
    
    @property
    def pairs_of_species(self):
        """
        Return a list of all the possible pairs of species.
        """
        pairs = []
        for s1 in self.distinct_species:
            for s2 in self.distinct_species:
                pairs.append((s1, s2))
        return pairs

    @property
    def pairs_of_species_id(self):
        """
        Return a list of all the possible pairs of species ID.
        """
        pairs = []
        for i in range(self.number_of_species):
            for j in range(self.number_of_species):
                pairs.append((i+1, j+1))
        return pairs
    
    @property
    def chemical_fractions(self):
        """
        Return a numpy array of the chemical fractions of each species 
        in the system.
        """
        species = self.dump('species')
        fractions = numpy.empty(self.number_of_species)
        for i, species_i in enumerate(self.distinct_species):
            fractions[i] = numpy.sum(species == species_i) / self.number_of_particles
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
            - "x" ("particle.position_x")
            - "y" ("particle.position_y")
            - "z" ("particle.position_z")
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
        
        to_dump : ndarray
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
            data = numpy.array([getattr(p, attr) for p in self.particle])
        elif what.startswith('cell'):
             data = numpy.array(getattr(self.cell, attr))
        else:
            raise ValueError('Unknown attribute %s' % what)
        return data

    def __str__(self):
        rep = 'System(number_of_particles={}, species={}, chemical_fractions={}, cell={})'
        return rep.format(self.number_of_particles,
                          self.distinct_species,
                          self.chemical_fractions,
                          self.cell.side)
    
    def __repr__(self):
        return self.__str__()