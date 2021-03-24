# Reference to atooms

import numpy

# Aliases for particles' properties
aliases = {'pos': 'particle.position',
           'x': 'particle.position_x',
           'y': 'particle.position_y',
           'z': 'particle.position_z',
           'spe': 'particle.species',
           'rad': 'particle.radius',
           'id': 'particle.index',
           'position': 'particle.position',
           'position_x': 'particle.position_x',
           'position_y': 'particle.position_y',
           'position_z': 'particle.position_z',
           'species': 'particle.species',
           'species_id': 'particle._species_id',
           'radius': 'particle.radius',
           'index': 'particle.index',
           'label': 'particle.label'}

class Particle:
    
    index = 0
    
    def __init__(self, position=None, species='A', label=-1, radius=0.5):
        if position is None:
            self.position = numpy.zeros(3)
        else:
            self.position = numpy.asarray(position)
        self.species = species
        # Cluster label
        self.label = label
        self.radius = radius
        # Index of the particle
        self.index = Particle.index
        Particle.index += 1
        # Numeral id of the species (given in a trajectory)
        self._species_id = None
        
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