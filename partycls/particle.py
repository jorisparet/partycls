"""
Point particles in a cartesian reference frame.

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
           'species_id': 'particle.species_id',
           'spe_id': 'particle.species_id',
           'label': 'particle.label',
           'mass': 'particle.mass',
           'radius': 'particle.radius',
           'nearest_neighbors': 'particle.nearest_neighbors',
           'neighbors': 'particle.nearest_neighbors',
           'neighbours': 'particle.nearest_neighbors',
           'voronoi_signature': 'particle.voronoi_signature',
           'signature': 'particle.voronoi_signature'}


class Particle:
    """
    A particle is defined by its position, its type, and additional attributes
    like a radius, a cluster label, a list of neighbors, etc.
    
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
        
    nearest_neighbors : list of int, default: None
        Indices of the particle's nearest neighbors in the System.
    
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
        
    nearest_neighbors : list of int
        Zero-based indices of the particle's nearest neighbors in the System.
    
    Examples
    --------
    
    >>> p = Particle([0.0, 0.0, 0.0], species='A')
    >>> p = Particle([0.0, 0.0], species='B')
    """

    def __init__(self, position=None, species='A', label=-1, radius=0.5, nearest_neighbors=None):
        if position is None:
            self.position = numpy.zeros(3)
        else:
            self.position = numpy.asarray(position)
        self.species = species
        # Cluster label
        self.label = label
        # Particle radius
        self.radius = radius
        # Neighbors
        self.nearest_neighbors = nearest_neighbors

    def __str__(self):
        rep = 'Particle('
        for attr, value in self.__dict__.items():
            if not attr.startswith('_'):
                rep += '{}={}, '.format(attr, value)
        rep = rep[:-2] + ')'
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
