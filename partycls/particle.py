"""
This class is inspired by the `atooms` framework authored by Daniele Coslovich
See https://framagit.org/atooms/atooms 
"""

import numpy

# Aliases for particles' properties
aliases = {'position': 'particle.position',
           'pos': 'particle.position',
           'position[0]': 'particle.position[0]',  
           'pos[0]': 'particle.position[0]',
           'x': 'particle.position[0]',
           'position[1]': 'particle.position[1]', 
           'pos[1]': 'particle.position[1]',
           'y': 'particle.position[1]',
           'position[2]': 'particle.position[2]',
           'pos[2]': 'particle.position[2]',
           'z': 'particle.position[2]',
           'species': 'particle.species',
           'spe': 'particle.species',
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
        
    radius : float, optional, defaut: 0.5
        Particle radius.
    
    Attributes
    ----------
    
    position : float array
        The position of the particle.
        
    species : str
        Particle type / species.
        
    label : int
        Cluster label of the particle. 
        
    radius : float
        Particle radius.
        
    index : int
        A unique index to identify the particle.
    
    Examples
    --------
    
    >>> p = Particle([0.0, 0.0, 0.0], species='A')
    >>> p = Particle([0.0, 0.0], species='B')
    """
    
    def __init__(self, position=None, species='A', label=-1, radius=0.5):
        if position is None:
            self.position = numpy.zeros(3)
        else:
            self.position = numpy.asarray(position)
        self.species = species
        # Cluster label
        self.label = label
        # Particle radois
        self.radius = radius
        # Index of the particle
        self.index = id(self)
        
    def __str__(self):
        rep = 'Particle('
        for attr, value in self.__dict__.items():
            if not attr.startswith('_'):
                rep += '{}={}, '.format(attr, value)
        rep = rep[:-2]+')'
        return rep
    
    def __repr__(self):
        return self.__str__()
    
    def fold(self, cell):
        """
        Fold the particle position into the central cell.

        Parameters
        ----------
        cell : Cell
            Simulation cell.

        Returns
        -------
        None

        """
        def _periodic_vector_unfolded(vec, box):
            return vec - numpy.rint(vec / box) * box
        self.position[:] = _periodic_vector_unfolded(self.position, cell.side)