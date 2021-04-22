"""
This class is inspired by the `atooms` framework authored by Daniele Coslovich
See https://framagit.org/atooms/atooms 
"""

from .system import System
from .particle import Particle
from .cell import Cell
import numpy

def tipify(s):
    """
    Convert a string into the best matching type.
    """
    try:
        return int(s)
    except ValueError:
        try:
            return float(s)
        except ValueError:
            return s


class Trajectory:
    """
    A trajectory is composed by one or several frames, each frame being an 
    instance of `System`. Trajectory instances are iterable. By default, only
    the positions and particle types are being read from the trajectory file.
    Additional particle properties in the file can be read using the 
    `additional_fields` parameter.
    
    Parameters
    ----------
    
    filename : str
        Path to the trajectory file to read.
        
    fmt : str, optional, default: 'xyz'
        Format of the trajectory.
        
    additional_fields : list of str, optional, default: []
        Additional fields (i.e. particle properties) to read from the trajectory.   
        
    first : int, optional, default: 0
        Index of the first frame to consider in the trajectory. Starts at zero.
        
    last : int, optional, default: None
        Index of the last frame to consider in the trajectory. Default is the 
        last frame.
        
    step : int, optional, default: 1
        Step between each frame to consider in the trajectory. For example,
        if `step=2`, one every two frames is read.
        
    
    Attributes
    ----------
    
    filename : str
        Name of the original trajectory file.
        
    fmt : str, optional, default: "xyz"
        Format of the original trajectory file.
        
    additional_fields : list, optional, default: []
        List of additional particle properties that were extracted from the
        original trajectory file.
    
    Examples
    --------
    
    >>> from pysc.trajectory import Trajectory
    >>> traj = Trajectory('trajectory.xyz', additional_fields=['mass'])
    """
    
    def __init__(self, filename, fmt='xyz', backend=None, additional_fields=[], first=0, last=None, step=1):
        self.filename = filename
        self.fmt = fmt
        self.backend = backend
        self.additional_fields = additional_fields
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
        """
        Add an instance of `System` to the trajectory.
        """
        self._systems.append(system)
        
    def remove(self, frame):
        """
        Remove the system at position `frame` from the trajectory.
        """
        self._systems.pop(frame)
    
    def dump(self, what):
        """
        Return a list of numpy arrays with the system property specified by 
        `what`. The list size is the number of systems in the trajectory.
        
        Parameters
        ----------
        
        what : str
            Requested system property.
        
            `what` must be of the form 
            "particle.<attribute>" or "cell.<attribute>". 
            
            The following aliases are allowed:
            - "pos" ("particle.position")
            - "position" ("particle.position")
            - "x" ("particle.position_x")
            - "y" ("particle.position_y")
            - "z" ("particle.position_z")
            - "spe" ("particle.species")
            - "species" ("particle.species")
            - "species_id" ("particle.species_id")
            - "radius" ("particle.radius")
            - "label" ("particle.label")
            - "index" ("particle.index")
            - "mass" ("particle.mass")
            - "box" ("cell.side")
        
        Returns
        -------
        
        to_dump : list of ndarrays with a size equal to the number of systems.
            List of the requested system property. Each element of the list
            is a ndarray from a frame (system) in the trajectory.
        
        Examples
        --------
        
        >>> traj = Trajectory('trajectory.xyz')
        >>> pos = traj.dump('position')
        >>> spe = traj.dump('species')
        """
        to_dump = []
        for system in self._systems:
            to_dump.append(system.dump(what))
        return to_dump
    
    def _read(self, first, last, step):
        
        # formats recognized by defaults
        if self.backend is None:
            if self.fmt == 'xyz':
                self._parser_xyz()
            elif self.fmt == 'rumd':
                self._parser_rumd()
            else:
                raise ValueError('"{}" format is not recognized natively. You may try again with a backend.'.format(self.fmt))
                
        # atooms backend
        elif self.backend == 'atooms':
            self._parser_atooms()
        
        # MDTraj backend
        elif self.backend == 'mdtraj':
            self._parser_mdtraj()
        
        # wrong backend
        else:
            raise ValueError('backend "{}" is not a recognized backend'.format(self.backend))
        
        # Sanity checks
        #  constant number of particles
        n_particles = set([sys.number_of_particles for sys in self._systems])
        assert(len(n_particles) == 1), 'the number of particles should be kept constant in the trajectory.'
        #  constant volume
        volumes = set([tuple(sys.cell.side) for sys in self._systems])
        assert(len(volumes) == 1), 'the volume of the cell should be kept constant in the trajectory.'
        
        # Slice the trajectory
        self._slice(first, last, step)
        
    def _parser_xyz(self):
        """
        Read the trajectory from a XYZ file and put 
        the different frames in a list of `System`.
        """

        def _system_info(info):
            """Information on the system (dimension, fields, etc.)"""
            import re
            default_fields = ['id', 'type', 'species', 'pos', 
                              'position', 'x', 'y', 'z']
            # loop over all properties
            for p in info:
                # search columns
                re_cols = re.search('^[C|c]olumns:(.+)$', p)
                # keep additional fields/columns for a later use
                if re_cols:
                    fields = re_cols.group(1).split(',')
                    fields = [field for field in fields if field not in default_fields]
                # search cell
                re_cell = re.search('^[C|c]ell:(.+)$', p)
                if re_cell:
                    cell = re_cell.group(1).split(',')
                    cell = [float(L) for L in cell]
            return cell, fields
        
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
                sides, other_fields = _system_info(frame_info)
                cell = Cell(sides)
                system = System(cell=cell)
                dimension = len(sides)

                # Look for additional fields
                if self.additional_fields:
                    starting_idx = dimension + 1
                    # community/cluster label
                    default_cluster_fields = ['cluster', 'community', 'label']
                    read_cluster_field = True in [cls_field in self.additional_fields for cls_field in default_cluster_fields]
                    if read_cluster_field:
                        cluster_field_mask = [cf in other_fields for cf in default_cluster_fields]
                        has_cluster_field = True in cluster_field_mask
                        if has_cluster_field:
                            cluster_field_name = default_cluster_fields[cluster_field_mask.index(True)]
                            cidx = other_fields.index(cluster_field_name)
                    # other additional fields
                    fields_to_read = []
                    fields_to_read_idx = []
                    for field in self.additional_fields:
                        if field not in default_cluster_fields:
                            fidx = other_fields.index(field)
                            fields_to_read.append(field)
                            fields_to_read_idx.append(fidx)
                
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

                    # create the Particle object
                    particle = Particle(position=p_pos, species=p_type)
                    # set the additional fields
                    if self.additional_fields:
                        if read_cluster_field and has_cluster_field:
                            particle.label = int(line[starting_idx+cidx])
                        for field_name, field_idx in zip(fields_to_read, fields_to_read_idx):
                            val = tipify(line[starting_idx+field_idx])
                            particle.__setattr__(field_name, val)                  
                    
                    # add the particle to the system
                    system.add_particle(particle)
                    
                # Add system to trajectory
                self.add_system(system)

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
                # TODO: the following line does not work on 3.6... but it does like this
                # re_cell = re.search('^(sim_box|boxLengths)=(\w+),(.+)$', p)
                # and then shifts groups below
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

    def _parser_atooms(self):
        
        try:
            from atooms.trajectory import Trajectory as AtoomsTrajectory
            
            _Trajectory = AtoomsTrajectory.formats[self.fmt]
            
            # Read additional fields if the trajectory format allows 
            if self.additional_fields:
                try:
                    atooms_traj = _Trajectory(self.filename, fields=['id', 'pos']+self.additional_fields)
                except TypeError:
                    print('This trajectory format does not support additional fields')
                    print('Warning: ignoring additional fields.')
                    self.additional_fields = []
                    atooms_traj = _Trajectory(self.filename)
            else:
                atooms_traj = _Trajectory(self.filename)

            # Fill the native Trajectory using atooms Trajectory
            for atooms_sys in atooms_traj:
                cell = Cell(atooms_sys.cell.side)
                system = System(cell=cell)
                for atooms_p in atooms_sys.particle:
                    pos = atooms_p.position.copy()
                    spe = atooms_p.species
                    particle = Particle(position=pos, species=spe)
                    # additional fields
                    for field in self.additional_fields:
                        value = atooms_p.__getattribute__(field)
                        particle.__setattr__(field, value)
                    system.add_particle(particle)
                self.add_system(system)
            atooms_traj.close()
                
        except ModuleNotFoundError:
            print('No `atooms` module found.')

    def _parser_mdtraj(self):
        
        try:
            import mdtraj as md
            md_traj = md.load(self.filename, top=self.fmt)
            for frame in range(md_traj.n_frames):
                cell = Cell(side=md_traj.unitcell_lengths[frame])
                system = System(cell=cell)
                for atom in range(md_traj.n_atoms):
                    pos = md_traj.xyz[frame, atom]
                    spe = md_traj[frame].topology.atom(atom).element.symbol
                    particle = Particle(position=pos, species=spe)
                    system.add_particle(particle)
                self.add_system(system)
        except ModuleNotFoundError:
            print('No `mdtraj` module found.')

    def _write(self, output_path, fmt='xyz', backend=None, additional_fields=[], precision=6):

        # formats recognized by defaults
        if backend is None:
            if fmt == 'xyz':
                self._write_xyz(output_path, additional_fields, precision)
            elif fmt == 'rumd':
                self._write_rumd(output_path, additional_fields, precision)
            else:
                raise ValueError('"{}" format is not recognized natively. You may try again with a backend.'.format(self.fmt))
                
        # atooms backend
        elif backend == 'atooms':
            self._write_atooms(output_path, fmt, additional_fields, precision)
        
        # MDTraj backend
        elif backend == 'mdtraj':
            self._write_mdtraj(output_path, fmt, additional_fields, precision)

    def _write_xyz(self, output_path, additional_fields, precision):
        if not output_path.endswith('.xyz'):
            output_path += '.xyz'
        with open(output_path, 'w') as file:
            for system in self._systems:
                file.write('{}\n'.format(system.number_of_particles))
                columns = 'columns:id,pos,'
                for field in additional_fields:
                    columns += '{},'.format(field)
                columns = columns[:-1]
                header = columns + ' cell:{}\n'
                file.write(header.format(','.join('{:.{}f}'.format(L, precision) for L in system.cell.side)))
                for particle in system.particle:
                    line = '{} '.format(particle.species)
                    line += '{} '.format(' '.join('{:.{}f}'.format(p_i, precision) for p_i in particle.position))
                    for field in additional_fields:
                        line += '{} '.format(particle.__getattribute__(field))
                    line += '\n'
                    file.write(line)

    #TODO: check if actually readable by RUMD
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
            
    def _write_atooms(self, output_path, fmt, fields, precision):
        try:
            from atooms.trajectory import Trajectory as AtoomsTrajectory
            from atooms.system import System as _System
            from atooms.system import Particle as _Particle
            from atooms.system import Cell as _Cell
            
            _Trajectory = AtoomsTrajectory.formats[fmt]

            # Write additional fields if the trajectory format allows 
            try:
                with _Trajectory(output_path, 'w', fields=['id', 'pos']+fields) as atooms_traj:
                    for n, system in enumerate(self._systems):
                        new_cell = _Cell(side=system.cell.side)
                        new_system = _System(cell=new_cell)
                        for particle in system.particle:
                            pos = particle.position
                            spe = particle.species
                            label = particle.label
                            new_particle = _Particle(species=spe, position=pos)
                            new_particle.label = label
                            new_system.particle.append(new_particle)
                        atooms_traj.write(new_system, step=n)
                        
            except TypeError:
                print('This trajectory format does not support additional fields (e.g. cluster labels)')
        
        except ModuleNotFoundError:
            print('No `atooms` module found.')        
            
    def _write_mdtraj(self, output_path, fmt, fields, precision):
        raise NotImplementedError('Writing output trajectories with the MDTraj backend is currently impossible')

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
