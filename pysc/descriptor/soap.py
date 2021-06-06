import numpy
from .descriptor import StructuralDescriptor


class SOAPDescriptor(StructuralDescriptor):
    """
    """

    name = 'SOAP'
    symbol = 'soap'

    def __init__(self, trajectory, **kwargs):        
        from dscribe.descriptors import SOAP
        
        StructuralDescriptor.__init__(self, trajectory)
        #self._species = s.distinct_species()
        self._species = ['H', 'C']
        # TODO: descriptor parameters
        self._kwargs = kwargs
        self.soap = SOAP(
            species=self._species,
            rcut=4.0,
            nmax=10,
            lmax=7,
            periodic=True,
            sparse=False
        )
        self.grid = range(self.soap.get_number_of_features())
        
    def compute(self):
        from ase import Atoms
        import atooms.trajectory as trj

        def system_to_ase_atoms(system):
            from ase import Atoms
            atoms = Atoms(system.dump('particle.species'),
                          system.dump('particle.position'),
                          cell=system.cell.side,
                          pbc=True)
            return atoms

        self.features = numpy.empty((self.size, self.n_features))
        row = 0
        for s in self.trajectory:
            # should be a decorator        
            for p in s.particle:
                if p.species == 'A':
                    p.species = 'H'
                else:
                    p.species = 'C'
            system = system_to_ase_atoms(s)
            features = self.soap.create(system)
            self.features[row: row+features.shape[0], :] = features
            row += features.shape[0]
            
        self.features = numpy.array(self.features)
        return self.features

    def normalize(self, dist):
        return dist
