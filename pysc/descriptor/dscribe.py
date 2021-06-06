import numpy
from .descriptor import StructuralDescriptor

__all__ = ['DscribeDescriptor', 'DscribeChemicalDescriptor']

def _system_to_ase_atoms(system, chemistry, pbc):
    from ase import Atoms
    if chemistry:
        atoms = Atoms(system.dump('particle.species'),
                      system.dump('particle.position'),
                      cell=system.cell.side,
                      pbc=pbc)
    else:
        atoms = Atoms(['H'] * len(system.particle),
                      system.dump('particle.position'),
                      cell=system.cell.side,
                      pbc=pbc)
    return atoms


class DscribeDescriptor(StructuralDescriptor):
    """
    Adapter for generic DScribe descriptors, without chemical species 
    information. Essentially, all the particles are considered as hydrogen
    atoms.
    """

    # Class-level switch to use chemical information
    _chemistry = False

    def __init__(self, trajectory, backend, *args, **kwargs):

        StructuralDescriptor.__init__(self, trajectory)
        self.name = backend.__name__
        self.symbol = backend.__name__.lower()

        # Use chemical species
        if self._chemistry:
            kwargs['species'] = self.trajectory[0].distinct_species
        else:
            kwargs['species'] = ['H']

        # Periodic boundary conditions
        cell = self.trajectory[0].cell
        if cell is not None:
            kwargs['periodic'] = cell.periodic.all()
            self._periodic = cell.periodic.all()
        else:
            kwargs['periodic'] = False
            self._periodic = False

        # DScribe backend setup
        self.backend = backend(*args, **kwargs)
        self.grid = range(self.backend.get_number_of_features())

    def compute(self):
        self.features = numpy.empty((self.size, self.n_features))
        row = 0
        for system in self.trajectory:
            system = _system_to_ase_atoms(system,
                                          chemistry=self._chemistry,
                                          pbc=self._periodic)
            features = self.backend.create(system)
            self.features[row: row+features.shape[0], :] = features
            row += features.shape[0]

        return self.features

    def normalize(self, dist):
        return dist

class DscribeChemicalDescriptor(DscribeDescriptor):
    """
    Adapter for generic DScribe descriptors, with chemical species information.
    """

    _chemistry = True
    
