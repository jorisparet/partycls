import numpy
from .descriptor import StructuralDescriptor

__all__ = ['DscribeDescriptor', 'DscribeChemicalDescriptor']

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
        # NaNs
        accept_nans_key = 'accept_nans'
        accept_nans = True
        if accept_nans_key in kwargs.keys():
            accept_nans = kwargs[accept_nans_key]
            kwargs.pop(accept_nans_key)
        # verbose
        verbose_key = 'verbose'
        verbose = False
        if verbose_key in kwargs.keys():
            verbose = kwargs[verbose_key]
            kwargs.pop(verbose_key)
        StructuralDescriptor.__init__(self, trajectory,
                                      accept_nans=accept_nans,
                                      verbose=verbose)
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
        # set up
        self._set_up(dtype=numpy.float64)
        row = 0
        # computation
        for i, system in enumerate(self._tqdm(self.trajectory)):
            positions = self.dump('position', group=1)[i]
            if self._chemistry:
                species = self.dump('species', group=1)[i]
            else:
                species = ['H'] * len(self.dump('species', group=1)[i])
            side = system.cell.side
            system = _arrays_to_ase_atoms(positions, species, side,
                                          pbc=self._periodic)
            other_positions = self.dump('position', group=0)[i]
            features = self.backend.create(system, positions=other_positions)
            self.features[row: row + features.shape[0], :] = features
            row += features.shape[0]
        self._handle_nans()
        return self.features

    def normalize(self, dist):
        return dist

    def _tqdm(self, iterable):
        if self.verbose == 1:
            try:
                from tqdm import tqdm
                return tqdm(iterable,
                            desc='Computing {} descriptor'.format(self.symbol))
            except ImportError:
                print('Warning: install tqdm to show the progress bar.')
                return iterable
        return iterable     


class DscribeChemicalDescriptor(DscribeDescriptor):
    """
    Adapter for generic DScribe descriptors, with chemical species information.
    """

    _chemistry = True
