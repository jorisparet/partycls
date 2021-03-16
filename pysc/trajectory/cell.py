# Reference to atooms

import numpy

class Cell:
    
    def __init__(self, side=None):
        if side is None:
            self.side = numpy.zeros(3)
        else:
            self.side = numpy.asarray(side, dtype=numpy.float64)
            
    @property
    def volume(self):
        return numpy.prod(self.side)