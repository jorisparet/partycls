import numpy
from .descriptor import StructuralDescriptor

__all__ = ['DscribeDescriptor', 'DscribeChemicalDescriptor']

def _system_to_ase_atoms(system, chemistry, pbc):
    from ase import Atoms
    if chemistry:
        atoms = Atoms(system.get_property('particle.species'),
                      system.get_property('particle.position'),
                      cell=system.cell.side,
                      pbc=pbc)
    else:
        atoms = Atoms(['H'] * len(system.particle),
                      system.get_property('particle.position'),
                      cell=system.cell.side,
                      pbc=pbc)
    return atoms

def _arrays_to_ase_atoms(positions, species, side, pbc):
    from ase import Atoms
    atoms = Atoms(species,
                  positions,
                  cell=side,
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
        StructuralDescriptor._sanity_checks(self)
        self.features = numpy.empty((self.size, self.n_features))
        row = 0
        for i, system in enumerate(self.trajectory):
            positions = self.dump('position', 1)[i]
            if self._chemistry:
                species = self.dump('species', 1)[i]
            else:
                species = ['H'] * len(self.dump('species', 1)[i])
            side = system.cell.side
            system = _arrays_to_ase_atoms(positions, species, side,
                                          pbc=self._periodic)
            other_positions = self.dump('position', 0)[i]
            features = self.backend.create(system, positions=other_positions)
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
    
