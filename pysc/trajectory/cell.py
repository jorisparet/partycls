"""
This class is inspired by the `atooms` framework authored by Daniele Coslovich
See https://framagit.org/atooms/atooms 
"""

import numpy

class Cell:
    """
    Rectangular bidimensional or tridimensional cell.
    
    Parameters
    ----------
    
    side : list of float or float array, optional, default: None
        List of lengths for the sides of the cell.
    
    Attributes
    ----------
    
    side : float array
        List of lengths for the sides of the cell.
    
    Examples
    --------
    
    >>> c = Cell(side=[2.0, 2.0, 2.0])
    >>> c.volume
    8.0
    """
    
    def __init__(self, side=None):
        if side is None:
            self.side = numpy.zeros(3)
        else:
            self.side = numpy.asarray(side, dtype=numpy.float64)
            
    @property
    def volume(self):
        """
        Return the volume of the cell.
        """
        return numpy.prod(self.side)