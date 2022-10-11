"""
Simulation cell.

This class is inspired by the framework `atooms <https://framagit.org/atooms/atooms>`_
authored by `Daniele Coslovich <https://www2.units.it/daniele.coslovich/>`_.
"""

import numpy


class Cell:
    """
    Orthorhombic cell.
    
    Attributes
    ----------
    side : numpy.ndarray
        List of lengths for the sides of the cell.
        
    periodic : numpy.ndarray
        Periodicity of the cell on each axis.
    """

    def __init__(self, side, periodic=None):
        """
        Parameters
        ----------
        side : list
            List of lengths for the sides of the cell.

        periodic : list, default: None
            Periodicity of the cell on each axis. Default is ``None`` (sets ``True``)
            in each direction.

        Example
        -------
        >>> c = Cell([2.0, 2.0, 2.0], periodic=[True, True, True ])
        """
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
