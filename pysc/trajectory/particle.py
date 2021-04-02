# Reference to atooms

import numpy

# Aliases for particles' properties
aliases = {'position': 'particle.position',
           'pos': 'particle.position',
           'position_x': 'particle.position_x',
           'x': 'particle.position_x',
           'position_y': 'particle.position_y',
           'y': 'particle.position_y',
           'position_z': 'particle.position_z',
           'z': 'particle.position_z',
           'species': 'particle.species',
           'spe': 'particle.species',
           'species_id': 'particle.species_id',
           'label': 'particle.label',
           'index': 'particle.index',
           'mass': 'particle.mass',
           'radius': 'particle.radius'}

class Particle:
    """
    A particle is defined by its position, its type, and an optional cluster
    label (default is -1).
    
    Parameters
    ----------
    
    position : list of float or float array, optional, default: None
        The position of the particle. 
        If not given, it will be set to [0.0, 0.0, 0.0].
        
    species : str, optional, default: "A"
        Particle type / species.
    
    label : int, optional, default: -1
        Cluster label of the particle. 
        Default is -1 (i.e. not belonging to any cluster).
    
    Attributes
    ----------
    
    position : float array
        The position of the particle.
        
    species : str
        Particle type / species.
        
    species_id : int
        A numeral ID for the species. Automatically given in the context
        of a trajectory.
        
    label : int
        Cluster label of the particle. 
        
    index : int
        A unique index to identify the particle.
    
    Examples
    --------
    
    
    """
    
    index = 0
    
    def __init__(self, position=None, species='A', label=-1):
        if position is None:
            self.position = numpy.zeros(3)
        else:
            self.position = numpy.asarray(position)
        self.species = species
        # Cluster label
        self.label = label
        # Index of the particle
        self.index = Particle.index
        Particle.index += 1
        # Numeral id of the species (given in a trajectory)
        self.species_id = None
        
    def __str__(self):
        rep = 'Particle('
        for attr, value in self.__dict__.items():
            if not attr.startswith('_'):
                rep += '{}={}, '.format(attr, value)
        rep = rep[:-2]+')'
        return rep
    
    def __repr__(self):
        return self.__str__()
    
    @property
    def position_x(self):
        return self.position[0]

    @property
    def position_y(self):
        return self.position[1]
    
    @property
    def position_z(self):
        return self.position[2]