# Reference to atooms

from .system import System
from .particle import Particle
from .cell import Cell
import numpy

class Trajectory:
    
    def __init__(self, filename, fmt='xyz', first=0, last=None, step=1):
        self.filename = filename
        self.fmt = fmt
        self._systems = []
        self._read(first, last, step)

    def __getitem__(self, item):
        return self._systems[item]
    
    def __len__(self):
        return len(self._systems)

    def __str__(self):
        rep = 'Trajectory(filename="{}", number_of_frames={})'
        return rep.format(self.filename, self.__len__())
    
    def __repr__(self):
        return self.__str__()

    def add_system(self, system):
        self._systems.append(system)
        
    def remove(self, item):
        self._systems.pop(item)
    
    def dump(self, what):
        to_dump = []
        for system in self._systems:
            to_dump.append(system.dump(what))
        return to_dump
    
    def _read(self, first, last, step):
        
        # Select the correct parser to read the file
        if self.fmt == 'xyz':
            self._parser_xyz()
        elif self.fmt == 'rumd':
            self._parser_rumd()
        else:
            pass
        
        # Add numeral ID to each species (1 to N_species)
        self._make_species_numeral()
        
        # Sanity checks
        #  constant number of particles
        n_particles = set([sys.number_of_particles for sys in self._systems])
        assert(len(n_particles) == 1), 'the number of particles should be kept constant in the trajectory.'
        #  constant volume
        volumes = set([tuple(sys.cell.side) for sys in self._systems])
        assert(len(volumes) == 1), 'the volume of the cell should be kept constant in the trajectory.'
        
        # Slice the trajectory
        self._slice(first, last, step)

    def _write(self, output_path, fmt='xyz', fields=[], precision=6):
        # Select the method matching `fmt` to write the file
        if fmt == 'xyz':
            self._write_xyz(output_path, fields, precision)
        elif fmt == 'rumd':
            self._write_rumd(output_path, fields, precision)
        else:
            raise ValueError('Unknown trajectory format "{}"'.format(fmt))
        
    def _parser_xyz(self):
        """
        Read the trajectory from a XYZ file and put 
        the different frames in a list of `System`.
        """

        def _system_info(info):
            """Information on the system (dimension, fields, etc.)"""
            import re
            # loop over all properties
            for p in info:
                # search columns
                re_cols = re.search('^[C|c]olumns:(.+)$', p)
                # keep fields/columns for a later use
                if re_cols:
                    fields = re_cols.group(1).split(',')
                # search cell
                re_cell = re.search('^[C|c]ell:(.+)$', p)
                if re_cell:
                    cell = re_cell.group(1).split(',')
                    cell = [float(L) for L in cell]
            # look for community/cluster field
            cluster_field = 'community' in fields or 'cluster' in fields
            return cell, cluster_field
        
        # Read whole file
        with open(self.filename) as trajectory:
            while True:
                firstline = trajectory.readline()
                # Stop if EOF
                if not(firstline):
                    break
                # Current frame
                n_particles = int(firstline)
                frame_info = trajectory.readline().split()
                sides, cluster_field = _system_info(frame_info)
                cell = Cell(sides)
                system = System(cell=cell)
                dimension = len(sides)
                
                # Loop over particles
                for p in range(n_particles):
                    line = trajectory.readline().split()                    
                    # particle type
                    p_type = line[0]
                    # position (2D or 3D)
                    if dimension == 2:
                        p_pos = numpy.array(line[1:3], dtype=numpy.float64)
                    if dimension == 3:
                        p_pos = numpy.array(line[1:4], dtype=numpy.float64)
                    # community/cluster
                    if cluster_field:
                        if dimension == 2:
                            p_label = int(line[3])
                        if dimension == 3:
                            p_label = int(line[4])
                    else:
                        p_label = -1    
                    # create the Particle object
                    particle = Particle(position=p_pos, species=p_type, label=p_label)
                    system.add_particle(particle)
                    
                # Add system to trajectory
                self.add_system(system)

    def _write_xyz(self, output_path, fields, precision):
        if not output_path.endswith('.xyz'):
            output_path += '.xyz'
        with open(output_path, 'w') as file:
            for system in self._systems:
                file.write('{}\n'.format(system.number_of_particles))
                columns = 'columns:id,pos,cluster'
                if 'radius' in fields: columns += ',radius'
                header = columns + ' cell:{}\n'
                file.write(header.format(','.join('{:.{}f}'.format(L, precision) for L in system.cell.side)))
                for particle in system.particle:
                    line = '{} '.format(particle.species)
                    line += '{} '.format(' '.join('{:.{}f}'.format(p_i, precision) for p_i in particle.position))
                    line += '{} '.format(particle.label)
                    if 'radius' in fields:
                        line += '{} '.format(particle.radius)
                    line += '\n'
                    file.write(line)

    def _parser_rumd(self):
        """
        Read the trajectory from a RUMD file and put 
        the different frames in a list of `System`.
        """
        
        import gzip
        
        def _system_info(info):
            """Information on the system (dimension, fields, etc.)"""
            import re
            # loop over all properties
            for p in info:
                # search columns
                re_cols = re.search('^columns=(.+)$', p)
                # keep fields/columns for a later use
                if re_cols:
                    fields = re_cols.group(1).split(',')
                # search cell
                re_cell = re.search('^[sim_box|boxLengths]=(\w+),(.+)$', p)
                if re_cell:
                    assert('Rectangular' in re_cell.group(1)), 'simulation box must be rectangular.'
                    cell = re_cell.group(2).split(',')
                    cell = [float(L) for L in cell]
            # look for community/cluster field
            cluster_field = 'community' in fields or 'cluster' in fields
            return cell, cluster_field
        
        with gzip.open(self.filename, mode='rt') as trajectory:
            while True:
                firstline = trajectory.readline()
                # Stop if EOF
                if not(firstline):
                    break
                # Current frame
                n_particles = int(firstline)
                frame_info = trajectory.readline().split()
                sides, cluster_field = _system_info(frame_info)
                cell = Cell(sides)
                system = System(cell=cell)
                dimension = len(sides)            

                # Loop over particles
                for p in range(n_particles):
                    line = trajectory.readline().split()                    
                    # particle type
                    p_type = line[0]
                    # position (2D or 3D)
                    if dimension == 2:
                        p_pos = numpy.array(line[1:3], dtype=numpy.float64)
                    if dimension == 3:
                        p_pos = numpy.array(line[1:4], dtype=numpy.float64)
                    # community/cluster
                    if cluster_field:
                        if dimension == 2:
                            p_label = int(line[3])
                        if dimension == 3:
                            p_label = int(line[4])
                    else:
                        p_label = -1    
                    # create the Particle object
                    particle = Particle(position=p_pos, species=p_type, label=p_label)                        
                    system.add_particle(particle)
                    
                # Add system to trajectory
                self.add_system(system)

    #TODO: check if actually readable by RUMD
    #TODO: allow to write radii
    def _write_rumd(self, output_path, fields, precision):
        import gzip
        if not output_path.endswith('.xyz.gz'):
            output_path += '.xyz.gz'
        with gzip.open(output_path, 'wt') as file:
            for system in self._systems:
                file.write('{}\n'.format(system.number_of_particles))
                header = 'ioformat=1 boxLengths={} '
                dimension = system.number_of_dimensions
                if dimension == 2:
                    columns = 'columns=type,x,y,cluster\n'
                if dimension == 3:
                    columns = 'columns=type,x,y,z,cluster\n'
                header = header + columns
                file.write(header.format(','.join('{:.{}f}'.format(L, precision) for L in system.cell.side)))
                for particle in system.particle:
                    line = '{} '.format(particle.species_id-1)
                    line += '{} '.format(' '.join('{:.{}f}'.format(p_i, precision) for p_i in particle.position))
                    line += '{} '.format(particle.label)
                    line += '\n'
                    file.write(line)
        
    def _make_species_numeral(self):
        """
        Standardize the names of the species to [1, ..., N_species] by
        changing the attribute `particle.species_id` of each particle in
        the trajectory.
        """
        for system in self._systems:
            distinct_species = list(system.distinct_species)
            for particle in system.particle:
                particle.species_id = distinct_species.index(particle.species) + 1

    # TODO: check if always working
    # TODO: handle fractions
    def _slice(self, first, last, step):
        nframes = len(self._systems)
        assert(first >= 0 and first < nframes), 'invalid first frame.'
        if last is not None:
            assert(last >= 0 and last < nframes), 'invalid last frame.'
            assert(first <= last), 'first frame must be inferior to last frame.'
        else:
            last = nframes
        # Remove unwanted frames
        all_frames = list(range(nframes))
        frames = all_frames[first:last+1:step]
        kept_systems = []
        for frame in frames:
            kept_systems.append(self._systems[frame])
        # Owerwrite non-sliced list
        self._systems = kept_systems