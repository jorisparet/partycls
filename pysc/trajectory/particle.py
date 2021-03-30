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