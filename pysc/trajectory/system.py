# Reference to atooms

from .particle import Particle, aliases
import numpy

class System:
    
    """System class."""
    
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
            raise TypeError
    
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
    def number_of_particles(self):
        return len(self.particle)
        
    @property
    def density(self):
        """
        Density of the system.

        It will raise a ValueException if `cell` is None.
        """
        if self.cell is None:
            return ValueError('cannot compute density without a cell')
        return len(self.particle) / self.cell.volume
    
    @property
    def distinct_species(self):
        """
        Sorted array of all the distinct species in the `System`.
        """
        return numpy.array(sorted(set(self.dump('species'))))
    
    @property
    def number_of_species(self):
        """
        Number of distinct species in the system.
        """
        return len(self.distinct_species)
    
    @property
    def pairs_of_species(self):
        """
        Array of all the possible pairs of species.
        """
        pairs = []
        for s1 in self.distinct_species:
            for s2 in self.distinct_species:
                pairs.append((s1, s2))
        return pairs

    @property
    def pairs_of_species_id(self):
        pairs = []
        for i in range(self.number_of_species):
            for j in range(self.number_of_species):
                pairs.append((i+1, j+1))
        return pairs
    
    @property
    def chemical_fractions(self):
        species = self.dump('species')
        fractions = numpy.empty(self.number_of_species)
        for i, species_i in enumerate(self.distinct_species):
            fractions[i] = numpy.sum(species == species_i) / self.number_of_particles
        return fractions

    def dump(self, what):
        """
        Return a numpy array with the system propety specified by `what`.
        
        `what` must be of the form 
        `particle.<attribute>` or `cell.<attribute>`. 
        
        The following aliases are allowed:
        - `pos` (`particle.position`)
        - `position` (particle.position)
        - `spe` (`particle.species`)
        - `species` (`particle.species`)
        - `species_id` (`particle._species_id`)
        - `rad` (`particle.radius`)
        - `radius` (`particle.radius`)
        - `box` (`cell.side`)
        - `id` (`particle.index`)
        - `index` (`particle.index`)
        - `label` (`particle.label`)
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
