"""
This class is inspired by the `atooms` framework authored by Daniele Coslovich
See https://framagit.org/atooms/atooms 
"""

import numpy

class Cell:
    """
    Orthorhombic cell.
    
    Parameters
    ----------
    
    side : list of float or float array
        List of lengths for the sides of the cell.
    
    Attributes
    ----------
    
    side : float array
        List of lengths for the sides of the cell. The default is None (will
        set `True` in each direction).
        
    periodic : bool array
        Periodicity of the cell on each axis.
    
    Examples
    --------
    
    >>> c = Cell([2.0, 2.0, 2.0])
    >>> c.volume
    8.0
    """
    
    def __init__(self, side, periodic=None):
        self.side = numpy.asarray(side, dtype=numpy.float64)

        # Periodic boundary conditions apply separately on each axis
        if periodic is None:
            self.periodic = numpy.empty_like(self.side, dtype=bool)
            self.periodic[:] = True
        else:
            self.periodic = numpy.asarray(periodic, dtype=bool)
            
    @property
    def volume(self):
        """
        Volume of the cell.
        """
        return numpy.prod(self.side)
    
    def __str__(self):
        return 'Cell(side={}, periodic={}, volume={})'.format(self.side,
                                                              self.periodic,
                                                              self.volume)
    
    def __repr__(self):
        return self.__str__()
