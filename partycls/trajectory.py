"""
This class is inspired by the `atooms` framework authored by Daniele Coslovich
See https://framagit.org/atooms/atooms 
"""

import numpy
from .system import System
from .particle import Particle
from .cell import Cell
from .core.utils import tipify

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
        Format of the trajectory. Needed when using "atooms" as a backend.
        
    backend : str, optional, default: None
        Name of a third-party package to use as backend when reading the input
        trajectory. Currently supports "atooms" and "mdtraj".
        
    top : str, mdtraj.Trajectory, or mdtraj.Topology, optional, defaut: None
        Topology information. Needed when using "mdtraj" as backend on a 
        trajectory file whose format requires topology information. See MDTraj
        documentation for more information.
         
    additional_fields : list of str, optional, default: []
        Additional fields (i.e. particle properties) to read from the 
        trajectory. Not all trajectory formats allow for additional fields.
        
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
        
    fmt : str
        Format of the original trajectory file.
        
    backend : str or None
        Name of the third-party package used to read the input trajectory file.
        
    additional_fields : list, default: None
        List of additional particle properties that were extracted from the
        original trajectory file.
    
    Examples
    --------
    
    >>> from partycls.trajectory import Trajectory
    >>> traj = Trajectory('trajectory.xyz', additional_fields=['mass'])
    """
    
    def __init__(self, filename, fmt=None, backend=None, top=None, additional_fields=None, first=0, last=None, step=1):
        self.filename = filename
        if backend is None and fmt is None:
            self.fmt = 'xyz'
        else:
            self.fmt = fmt
        self.backend = backend
        self.top = top
        if additional_fields is None:
            self.additional_fields = []
        else:
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
        
    def remove(self, frame):
        """
        Remove the system at position `frame` from the trajectory.

        Parameters
        ----------
        frame : int
            Index of the frame to remove from the trajectory.

        Returns
        -------
        None.

        """
        self._systems.pop(frame)
    
    def get_property(self, what, subset=None):
        """
        Return a list of numpy arrays with the system property specified by 
        `what`. The list size is the number of systems in the trajectory.

        Parameters
        ----------
        what : str
            Requested system property.
            
            `what` must be of the form 
            "particle.<attribute>" or "cell.<attribute>". 
            
            The following particle aliases are accepted:
            - 'position': 'particle.position'
            - 'pos': 'particle.position'
            - 'position[0]': 'particle.position[0]', 
            - 'pos[0]': 'particle.position[0]'
            - 'x': 'particle.position[0]'
            - 'position[1]': 'particle.position[1]',
            - 'pos[1]': 'particle.position[1]'
            - 'y': 'particle.position[1]'
            - 'position[2]': 'particle.position[2]'
            - 'pos[2]': 'particle.position[2]'
            - 'z': 'particle.position[2]'
            - 'species': 'particle.species'
            - 'spe': 'particle.species'
            - 'label': 'particle.label'
            - 'index': 'particle.index'
            - 'mass': 'particle.mass'
            - 'radius': 'particle.radius'

        subset : str, optional
            Subset of particles for which the property must be dumped. Must be 
            of the form "particle.<attribute>" unless "<attribute>" is an 
            alias. The default is None (all particles will be included).
            This is ignored if `what` is cell property.

        Returns
        -------
        to_dump : list
            List of the requested system property with length equal to the 
            number of frames in the trajectory. Each element of the list is a
            numpy.ndarray of the requested system property.
            
        Examples
        --------
        >>> traj = Trajectory('trajectory.xyz')
        >>> pos = traj.get_property('position')
        >>> spe = traj.get_property('species')
        >>> sides = traj.get_property('cell.side')
        
        """
        to_dump = []
        for system in self._systems:
            to_dump.append(system.get_property(what, subset))
        return to_dump

    def dump(self, what, subset=None):
        """
        Alias for the method get_property().
        """
        return self.get_property(what, subset=subset)
    
    def set_property(self, what, value, subset=None):
        """
        Set a property `what` to `value` for all the particles in the 
        trajectory or for a given subset of particles specified by `subset`.

        Parameters
        ----------
        what : str
            Name of the property to set. This is considered to be a particle
            property by default, unless it starts with "cell", e.g. 
            "cell.side".
            
        value : int, float, list, or numpy.ndarray
            Value(s) of the property to set. An instance of `int` or `float`
            will set the same value for all concerned particles. An instance
            of `list` or `numpy.ndarray` will assign a specific value to each
            particle. In this case, the shape of `value` should respect the
            number of frames in the trajectory and the number of concerned
            particles.
            
        subset : str, optional
            Particles to which the property must be set. The default is None.
            This is ignored if `what` is cell property.

        Returns
        -------
        None.

        Examples
        --------
        >>> traj.set_property('mass', 1.0)
        >>> traj.set_property('radius', 0.5, "species == 'A'")
        >>> labels = [[0, 1, 0], # 2 frames, 3 particles in the subset
                      [1, 1, 0]]
        >>> traj.set_property('label', labels, "species == 'B'")

        """
        if not isinstance(value, (list, numpy.ndarray)):
            for system in self._systems:
                system.set_property(what, value, subset=subset)
        else:
            assert len(value) == self.__len__(), '`value` should have the same length than the Trajectory.'
            for frame, system in enumerate(self._systems):
                system.set_property(what, value[frame], subset=subset)

    def show(self, frames=None, backend='matplotlib', color='species', **kwargs):
        """
        Show the frames on index `frames` of the trajectory and color particles
        according to an arbitrary property, such as species, cluster label, 
        etc. Current visualization backends are 'matplotlib' and '3dmol'.

        Parameters
        ----------
        frames : list of int, optional
            Indices of the frames to show. The default is None (shows all frames).
        backend : str, optional
            Name of the backend to use for visualization. 
            The default is 'matplotlib'.
        color : str, optional
            Name of the particle property to use as basis for coloring the 
            particles. This property must be defined for all the particles in 
            the system. The default is 'species'.
        **kwargs : additional keyworded arguments (backend-dependent).

        Raises
        ------
        ValueError
            In case of unknown `backend`.

        Returns
        -------
        list of Figure or View (backend-dependent)
        
        Examples
        --------
        >>> traj.show(frames=[0,1,2], color='label', backend='3dmol')
        >>> traj.show(frames=[0,1], color='energy', backend='matplotlib', cmap='viridis')
        >>> traj[0].show() # use the iterability of Trajectory objects

        """
        # show all frames (default)
        if frames is None:
            frames = range(len(self))
        # list of figures/views returned by each system
        snapshots = []
        for frame in frames:
            kwargs_f = kwargs.copy()
            if 'outfile' in kwargs:
                kwargs_f['outfile'] += '{:04}'.format(frame)
            snapshot = self._systems[frame].show(backend=backend, 
                                                 color=color,
                                                 **kwargs_f)
            snapshots.append(snapshot)
        return snapshots

    def fold(self):
        """
        Fold the particle positions into the central cell.

        Returns
        -------
        None.

        """        
        for system in self._systems:
            system.fold()
    
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
        
        # # Standardize the species by giving each a numeral ID
        self._make_species_numeral()
        
        # Sanity checks
        #  constant number of particles
        n_particles = set([len(sys.particle) for sys in self._systems])
        assert(len(n_particles) == 1), 'the number of particles should be kept constant in the trajectory.'
        
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
            default_fields = ['id', 'type', 'name', 'species', 'pos', 
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
                if not firstline:
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
                for _ in range(n_particles):
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
                    system.particle.append(particle)
                    
                # Add system to trajectory
                self._systems.append(system)

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
                re_cell = re.search('^(sim_box|boxLengths)=(.+)$', p)
                if re_cell:
                    cell = re_cell.group(2).split(',')
                    cell = [float(L) for L in cell if L[0].isdigit()]
            # look for community/cluster field
            cluster_field = 'community' in fields or 'cluster' in fields
            return cell, cluster_field
        
        with gzip.open(self.filename, mode='rt') as trajectory:
            while True:
                firstline = trajectory.readline()
                # Stop if EOF
                if not firstline:
                    break
                # Current frame
                n_particles = int(firstline)
                frame_info = trajectory.readline().split()
                sides, cluster_field = _system_info(frame_info)
                cell = Cell(sides)
                system = System(cell=cell)
                dimension = len(sides)            

                # Loop over particles
                for _ in range(n_particles):
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
                    system.particle.append(particle)
                    
                # Add system to trajectory
                self._systems.append(system)

    def _parser_atooms(self):
        
        try:
            from atooms.trajectory import Trajectory as AtoomsTrajectory
            
            supported = list(AtoomsTrajectory.formats.keys())
            assert self.fmt in supported, 'the current version of atooms only supports the following formats: {}'.format(supported)
            
            _Trajectory = AtoomsTrajectory.formats[self.fmt]
            
            # Read additional fields if the trajectory format allows 
            if self.additional_fields:
                try:
                    atooms_traj = _Trajectory(self.filename, mode='r', fields=['id', 'pos']+self.additional_fields)
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
                    system.particle.append(particle)
                self._systems.append(system)
            atooms_traj.close()
                
        except ModuleNotFoundError:
            raise ModuleNotFoundError('No `atooms` module found.')

    def _parser_mdtraj(self):
        
        try:
            import mdtraj as md
            try:
                md_traj = md.load(self.filename)
            except ValueError:
                md_traj = md.load(self.filename, top=self.top)
            for frame in range(md_traj.n_frames):
                input_cell = md_traj.unitcell_lengths
                assert input_cell is not None, 'cell dimensions are needed to read the trajectory'
                cell = Cell(side=input_cell[frame])
                system = System(cell=cell)
                for atom in range(md_traj.n_atoms):
                    pos = md_traj.xyz[frame, atom]
                    spe = md_traj[frame].topology.atom(atom).element.symbol
                    # virtual site
                    if spe == 'VS':
                        spe = md_traj[frame].topology.atom(atom).name
                    particle = Particle(position=pos, species=spe)
                    system.particle.append(particle)
                self._systems.append(system)
        except ModuleNotFoundError:
            raise ModuleNotFoundError('No `mdtraj` module found.')

    def write(self, output_path, fmt='xyz', backend=None, additional_fields=None, precision=6):
        """
        Write the current trajectory to a file.

        Parameters
        ----------
        output_path : str
            Name of the output trajectory file.
        fmt : str, optional
            Format of the output trajectory file. The default is 'xyz'.
        backend : str, optional
            Name of a third-party package to use when writing the output
            trajectory. The default is None.
        additional_fields : list of str, optional
            Additional fields (i.e. particle properties) to write in the output
            trajectory. Not all trajectory formats allow for additional fields. 
            The default is [].
        precision : int, optional
            Number of decimals when writing the output trajectory. 
            The default is 6.

        Raises
        ------
        ValueError
            - If `backend=None` and `fmt` is not recognized natively.
            - If `backend` is unknown.

        Returns
        -------
        None.

        """
        if additional_fields is None:
            additional_fields = []
        else:
            additional_fields = additional_fields

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
        
        # wrong backend
        else:
            raise ValueError('backend "{}" is not a recognized backend'.format(self.backend))

    def _write_xyz(self, output_path, additional_fields, precision):
        if not output_path.endswith('.xyz'):
            output_path += '.xyz'
        with open(output_path, 'w') as file:
            for system in self._systems:
                file.write('{}\n'.format(len(system.particle)))
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

    # /!\ RUMD does not seem to accept additional fields
    def _write_rumd(self, output_path, fields, precision):
        import gzip
        if not output_path.endswith('.xyz.gz'):
            output_path += '.xyz.gz'
        with gzip.open(output_path, 'wt') as file:
            for frame, system in enumerate(self._systems):
                file.write('{}\n'.format(len(system.particle)))
                header = 'ioformat=1 dt=0.001 timeStepIndex={} boxLengths={} '
                dimension = system.n_dimensions
                if dimension == 2:
                    columns = 'columns=type,x,y\n'
                if dimension == 3:
                    columns = 'columns=type,x,y,z\n'
                header = header + 'numTypes={} mass={} '
                header = header + columns
                file.write(header.format(frame,
                                         ','.join('{:.{}f}'.format(L, precision) for L in system.cell.side),
                                         len(system.distinct_species),
                                         ','.join('1.0' for s in system.distinct_species)))
                for particle in system.particle:
                    line = '{} '.format(particle.species_id-1)
                    line += '{} '.format(' '.join('{:.{}f}'.format(p_i, precision) for p_i in particle.position))
                    # no additional field
                    #line += '{} '.format(particle.label)
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
