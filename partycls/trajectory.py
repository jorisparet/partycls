"""
Physical trajectory.

This class is inspired by the framework `atooms <https://framagit.org/atooms/atooms>`_
authored by `Daniele Coslovich <https://www2.units.it/daniele.coslovich/>`_.
"""

import numpy
from .system import System
from .particle import Particle, aliases
from .cell import Cell
from .core.utils import tipify, _nearest_neighbors_methods_


class Trajectory:
    """
    A trajectory is composed by one or several frames, each frame being an 
    instance of ``System``. ``Trajectory`` instances are iterable. By default, only
    the positions and particle types are being read from the trajectory file.
    Additional particle properties in the file can be read using the 
    ``additional_fields`` parameter.
        
    Attributes
    ----------
    filename : str
        Name of the original trajectory file.
        
    fmt : str
        Format of the original trajectory file.
        
    backend : str
        Name of the third-party package used to read the input trajectory file.
        
    additional_fields : list
        List of additional particle properties that were extracted from the
        original trajectory file.
    """

    # accepted names for neighbors when reading/writing
    _default_neighbors_fields = ['neighbor', 'neighbors', 
                                 'neighbour', 'neighbours', 
                                 'nearest_neighbors',
                                 'nearest_neighbours']

    def __init__(self, filename, fmt=None, backend=None, top=None,
                 additional_fields=None, first=0, last=None, step=1):
        """
        Parameters
        ----------
        filename : str
            Path to the trajectory file to read.
            
        fmt : str, default: "xyz"
            Format of the trajectory. Needed when using ``"atooms"`` as a backend.
            
        backend : str, default: None
            Name of a third-party package to use as backend when reading the input
            trajectory. Currently supports ``"atooms"`` and ``"mdtraj"``.
            
        top : str, mdtraj.Trajectory, or mdtraj.Topology, defaut: None
            Topology information. Needed when using ``"mdtraj"`` as backend on a 
            trajectory file whose format requires topology information. See MDTraj
            documentation for more information.
            
        additional_fields : list, optional, default: None
            Additional fields (*i.e.* particle properties) to read from the 
            trajectory. Not all trajectory formats allow for additional fields.
            
        first : int, default: 0
            Index of the first frame to consider in the trajectory. Starts at zero.
            
        last : int, default: None
            Index of the last frame to consider in the trajectory. Default is the 
            last frame.
            
        step : int, default: 1
            Step between each frame to consider in the trajectory. For example,
            if ``step=2``, one out of every two frames is read.

        Examples
        --------
        >>> traj = Trajectory('trajectory.xyz', additional_fields=['mass'])
        >>> traj = Trajectory('trajectory.dat', fmt='lammps', backend='atooms')
        """
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
        # nearest neighbors
        self.nearest_neighbors_method = 'auto'
        self.nearest_neighbors_cutoffs = [None for pair in self._systems[0].pairs_of_species]

    @property 
    def nearest_neighbors_method(self):
        """
        Method used to identify the nearest neighbors of all the particles in
        the trajectory. Should be one of ``"auto"``, ``"fixed"``, ``"sann"`` or 
        ``"voronoi"``.
        """
        return self._nearest_neighbors_method

    @nearest_neighbors_method.setter
    def nearest_neighbors_method(self, value):
        value = value.lower()
        if value in _nearest_neighbors_methods_:
            self._nearest_neighbors_method = value
        else:
            raise ValueError('Invalid method for nearest neighbors. Should be one of {}.'.format(_nearest_neighbors_methods_))

    @property
    def nearest_neighbors_cutoffs(self):
        """
        List of cutoffs that delimit the first coordination shell. Cutoffs are 
        usually defined on the basis of the first minimum of the partial radial 
        distribution function of each pair of species, 
        :math:`g_{\\alpha\\beta}(r)`. The list must have the same length as the 
        number of pairs of species in the system (*e.g.* 2 species yield 4 
        possible pairs, 3 species yield 6 pairs, etc.).
        """
        return self._nearest_neighbors_cutoffs

    @nearest_neighbors_cutoffs.setter
    def nearest_neighbors_cutoffs(self, value):
        n_pairs = len(self._systems[0].pairs_of_species)
        if len(value) != n_pairs:
            raise ValueError("Incorrect number of cutoffs: {} were provided while {} are expected.".format(
                len(value), n_pairs))
        else:
            self._nearest_neighbors_cutoffs = value

    def remove(self, frame):
        """
        Remove the system at position ``frame`` from the trajectory.

        Parameters
        ----------
        frame : int
            Index of the frame to remove from the trajectory.

        Returns
        -------
        None
        """
        self._systems.pop(frame)

    def get_property(self, what, subset=None):
        """
        Return a list of numpy.ndarrays with the system property specified by 
        ``what``. The list size is the number of systems in the trajectory.

        Parameters
        ----------
        what : str
            Requested system property. ``what`` must be of the form 
            ``"particle.<attribute>"`` or ``"cell.<attribute>"``.
            The following particle aliases are accepted:

            - ``'position'`` : ``'particle.position'``
            - ``'pos'`` : ``'particle.position'``
            - ``'position[0]'`` : ``'particle.position[0]'``
            - ``'pos[0]'`` : ``'particle.position[0]'``
            - ``'x'`` : ``'particle.position[0]'``
            - ``'position[1]'`` : ``'particle.position[1]'``
            - ``'pos[1]'`` : ``'particle.position[1]'``
            - ``'y'`` : ``'particle.position[1]'``
            - ``'position[2]'`` : ``'particle.position[2]'``
            - ``'pos[2]'`` : ``'particle.position[2]'``
            - ``'z'`` : ``'particle.position[2]'``
            - ``'species'`` : ``'particle.species'``
            - ``'spe'`` : ``'particle.species'``
            - ``'label'`` : ``'particle.label'``
            - ``'mass'`` : ``'particle.mass'``
            - ``'radius'`` : ``'particle.radius'``
            - ``'nearest_neighbors'`` : ``'particle.nearest_neighbors'``
            - ``'neighbors'`` : ``particle.nearest_neighbors'``
            - ``'neighbours'`` : ``'particle.nearest_neighbors'``
            - ``'voronoi_signature'`` : ``'particle.voronoi_signature'``
            - ``'signature'`` : ``'particle.voronoi_signature'``

        subset : str, optional, default: ``None``
            Subset of particles for which the property must be dumped. Must be 
            of the form ``"particle.<attribute>"`` unless ``"<attribute>"`` is an 
            alias. The default is ``None`` (all particles will be included).
            This is ignored if ```what``` is cell property.

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
        Alias for the method ``get_property``.
        """
        return self.get_property(what, subset=subset)

    def set_property(self, what, value, subset=None):
        """
        Set a property ``what`` to ``value`` for all the particles in the 
        trajectory or for a given subset of particles specified by ``subset``.

        Parameters
        ----------
        what : str
            Name of the property to set. This is considered to be a particle
            property by default, unless it starts with ``"cell"``, *e.g.*
            ``"cell.side"``.
            
        value : int, float, list, numpy.ndarray
            Value(s) of the property to set. An instance of ``int`` or ``float``
            will set the same value for all concerned particles. An instance
            of ``list`` or ``numpy.ndarray`` will assign a specific value to each
            particle. In this case, the shape of ``value`` should respect the
            number of frames in the trajectory and the number of concerned
            particles.
            
        subset : str, default: None
            Particles to which the property must be set. The default is ``None``.
            This is ignored if ``what`` is a cell property.

        Returns
        -------
        None

        Examples
        --------
        >>> traj.set_property('mass', 1.0)
        >>> traj.set_property('radius', 0.5, subset="species == 'A'")
        >>> labels = [[0, 1, 0], # 2 frames, 3 particles in the subset
                      [1, 1, 0]]
        >>> traj.set_property('label', labels, subset="species == 'B'")
        """
        if not isinstance(value, (list, numpy.ndarray)):
            for system in self._systems:
                system.set_property(what, value, subset=subset)
        else:
            assert len(value) == self.__len__(), '`value` should have the same length as the Trajectory.'
            for frame, system in enumerate(self._systems):
                system.set_property(what, value[frame], subset=subset)

    def compute_nearest_neighbors(self, method=None, cutoffs=None, dr=0.1):
        """
        Compute the nearest neighbors for all the particles in the trajectory using
        the provided method. Neighbors are stored in the ``nearest_neighbors`` particle 
        property. Available methods are:

        - ``'auto'`` : read neighbors from the trajectory file, if explicitly requested with the ``additional_fields`` argument in the constructor.
        - ``'fixed'`` : use fixed cutoffs for each pair of species in the trajectory.
        - ``'sann'`` : solid-angle based nearest neighbor algorithm (see https://doi.org/10.1063/1.4729313).
        - ``'voronoi'`` : radical Voronoi tessellation method (uses particles' radii) (see https://doi.org/10.1016/0022-3093(82)90093-X)

        Parameters
        ----------
        method : str, default: None
            Method to identify the nearest neighbors. Must be one of 
            ``'auto'``, ``'fixed'``, ``'sann'``, or ``'voronoi'``. ``None`` defaults to ``'auto'``. If 
            method is ``'auto'``, neighbors are read directly from the trajectory file,
            if specified with the ``additional_fields`` argument in the constructor.
            If no neighbors are found, falls back to ``method='fixed'`` instead.

        cutoffs : list, default: None
            List containing the cutoffs distances for each pair of species
            in the trajectory (for method ``'fixed'`` and ``'sann'``). If ``None``, cutoffs
            will be computed automatically. For method ``'sann'``, cutoffs are
            required as a first guess to identify the nearest neighbors.

        dr : float, default: 0.1
            Radial grid spacing :math:`\Delta r` for computing the cutoffs on the basis
            of the first minimum of each partial radial distribution function
            in the trajectory, if cutoffs are not provided.

        Returns
        -------
        None

        Examples
        --------
        >>> traj.compute_nearest_neighbors(method='fixed', cutoffs=[1.5, 1.4, 1.4, 1.3])
        >>> traj.compute_nearest_neighbors(method='sann', cutoffs=[1.5, 1.4, 1.4, 1.3])
        >>> traj.compute_nearest_neighbors(method='voronoi')
        """

        # Convert method to Enum
        if method is not None:
            try:
                self.nearest_neighbors_method = method
            except ValueError:
                raise ValueError('Invalid method for nearest neighbors. Should be one of {}.'.format(_nearest_neighbors_methods_))

        # Read neighbors from the trajectory file if method is 'auto' (default)
        if self.nearest_neighbors_method == 'auto':
            neighbors = self.dump('nearest_neighbors') 
            # if no neighbors in the system, fallback to fixed-cutoffs method
            if None in neighbors[0]:
                self.compute_nearest_neighbors(method='fixed',
                                               cutoffs=cutoffs)
            return        

        # Compute cutoffs if not provided (for 'fixed' and 'sann')
        if self.nearest_neighbors_method in ['fixed', 'sann']:
            if cutoffs is None:
                if None in self.nearest_neighbors_cutoffs:
                    self.compute_nearest_neighbors_cutoffs()
            else:
                self.nearest_neighbors_cutoffs = cutoffs

        # Iterate over the systems
        for system in self._systems:
            system.compute_nearest_neighbors(self.nearest_neighbors_method,
                                             self.nearest_neighbors_cutoffs)

    def set_nearest_neighbors_cutoff(self, s_a, s_b, rcut, mirror=True):
        """
        Set the nearest-neighbor cutoff for the pair of species ``(s1, s2)``. The
        cutoff of the mirror pair ``(s2, s1)`` is set automatically if the ``mirror`` 
        parameter is ``True`` (default). Writes in the ``nearest_neighbors_cutoffs`` 
        list attribute.

        Parameters
        ----------
        s_a : str
            Symbol of the first species :math:`\\alpha`.
        s_b : str
            Symbol of the second species :math:`\\beta`.
        rcut : float
            Value of the cutoff for the pair :math:`(\\alpha,\\beta) r` = ``(s_a, s_b)``.
        mirror : bool, default: None
            Set the cutoff for the mirror pair ``(s_a, s_b)``. The default is ``True``.

        Returns
        -------
        None
        """
        pairs = self._systems[0].pairs_of_species
        idx_ab = pairs.index((s_a, s_b))
        self.nearest_neighbors_cutoffs[idx_ab] = rcut
        if mirror:
            idx_ba = pairs.index((s_b, s_a))
            self.nearest_neighbors_cutoffs[idx_ba] = rcut

    def compute_nearest_neighbors_cutoffs(self, dr=0.1):
        """
        Compute the nearest neighbors cutoffs on the basis of the first
        minimum of the partial radial distribution function 
        :math:`g_{\\alpha\\beta}(r)` between each pair of species 
        :math:`(\\alpha,\\beta)` in the trajectory. Sets the 
        ``nearest_neighbors_cutoffs`` list attribute.

        Parameters
        ----------
        dr : float, default: 0.1
            Bin width :math:`\Delta r` for the radial grid used to compute the partial
            radial distribution functions :math:`g_{\\alpha\\beta}(r)`.

        Returns
        -------
        None
        """
        from .descriptors import RadialDescriptor
        pairs = self._systems[0].pairs_of_species
        for pair in pairs:
            if self.nearest_neighbors_cutoffs[pairs.index(pair)] is None:
                s_a, s_b = pair
                # use the smallest side of the smallest box in case of
                #  non-constant volume trajectory
                L = numpy.min(self.dump('cell.side'))
                bounds = (0.0, L / 2)
                descriptor = RadialDescriptor(self, dr=dr, bounds=bounds)
                descriptor.add_filter("species == '{}'".format(s_a), group=0)
                descriptor.add_filter("species == '{}'".format(s_b), group=1)
                descriptor.compute()
                # grid and average descriptor
                r = descriptor.grid
                h_ab = descriptor.average
                # normalized g(r)
                g_ab = descriptor.normalize(h_ab, method="gr")
                # find the first minimum of g_ab(r)
                first_max = numpy.argmax(g_ab)
                first_min = numpy.argmin(g_ab[first_max:]) + first_max
                rcut = r[first_min]
                # set the cutoff
                self.set_nearest_neighbors_cutoff(s_a, s_b, rcut)

    def compute_voronoi_signatures(self):
        """
        Compute the Voronoi signatures of all the particles in the trajectory
        using the radical Voronoi tessellation method (see 
        https://doi.org/10.1016/0022-3093(82)90093-X).).
        
        Particle radii must be set using the ``set_property`` method if the 
        original trajectory file does not contain such information.

        Creates a ``voronoi_signature`` property for the particles.
        
        Returns
        -------
        None
        """
        for system in self._systems:
            system.compute_voronoi_signatures()

    def show(self, frames=None, backend='matplotlib', color='species', **kwargs):
        """
        Show the frames on index ``frames`` of the trajectory and color particles
        according to an arbitrary property, such as species, cluster label, 
        etc. Current visualization backends are ``"matplotlib"``, ``"ovito"``,
        and ``"3dmol"``.

        Parameters
        ----------
        frames : list, default: None
            Indices of the frames to show. The default is ``None`` (shows all frames).

        backend : str, default: "matplotlib"
            Name of the backend to use for visualization.

        color : str, default: "species"
            Name of the particle property to use as basis for coloring the 
            particles. This property must be defined for all the particles in 
            the system.

        **kwargs : additional keyworded arguments (backend-dependent).

        Raises
        ------
        ValueError
            In case of unknown ``backend``.

        Returns
        -------
        Backend-dependent
        
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

    def write(self, output_path, fmt='xyz', backend=None, additional_fields=None, precision=6):
        """
        Write the current trajectory to a file.

        Parameters
        ----------
        output_path : str
            Name of the output trajectory file.

        fmt : str, default: "xyz"
            Format of the output trajectory file.

        backend : str, default: None
            Name of a third-party package to use when writing the output
            trajectory.

        additional_fields : list, default: None
            Additional fields (*i.e.* particle properties) to write in the output
            trajectory. Not all trajectory formats allow for additional fields. 
            The default is to not write any additional particle property.

        precision : int, default: 6
            Number of decimals when writing the output trajectory.

        Raises
        ------
        ValueError
            - If ``backend=None`` and ``fmt`` is not recognized natively.
            - If ``backend`` is unknown.

        Returns
        -------
        None
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
                raise ValueError(
                    '"{}" format is not recognized natively. You may try again with a backend.'.format(self.fmt))

        # atooms backend
        elif backend == 'atooms':
            self._write_atooms(output_path, fmt, additional_fields, precision)

        # MDTraj backend
        elif backend == 'mdtraj':
            self._write_mdtraj(output_path, fmt, additional_fields, precision)

        # wrong backend
        else:
            raise ValueError('backend "{}" is not a recognized backend'.format(self.backend))

    def fold(self):
        """
        Fold the particles' positions into the central cell.

        Returns
        -------
        None
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
                raise ValueError(
                    '"{}" format is not recognized natively. You may try again with a backend.'.format(self.fmt))

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
                    read_cluster_field = True in [
                        cls_field in self.additional_fields for cls_field in default_cluster_fields]
                    if read_cluster_field:
                        cluster_field_mask = [cf in other_fields for cf in default_cluster_fields]
                        has_cluster_field = True in cluster_field_mask
                        if has_cluster_field:
                            cluster_field_name = default_cluster_fields[cluster_field_mask.index(True)]
                            cidx = other_fields.index(cluster_field_name)
                            
                    read_neigh_field = True in [
                        ngh_field in self.additional_fields for ngh_field in self._default_neighbors_fields]
                    if read_neigh_field:
                        neigh_field_mask = [nf in other_fields for nf in self._default_neighbors_fields]
                        has_neighbors_field = True in neigh_field_mask
                        if has_neighbors_field:
                            neighbors_field_name = self._default_neighbors_fields[neigh_field_mask.index(True)]
                            nidx = other_fields.index(neighbors_field_name)
                            
                    # other additional fields
                    fields_to_read = []
                    fields_to_read_idx = []
                    for field in self.additional_fields:
                        if field not in [*default_cluster_fields, *self._default_neighbors_fields]:
                            fidx = other_fields.index(field)
                            fields_to_read.append(field)
                            fields_to_read_idx.append(fidx)

                # Loop over particles
                for n in range(n_particles):
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
                    particle._index = n
                    # set the additional fields
                    if self.additional_fields:
                        if read_cluster_field and has_cluster_field:
                            particle.label = int(line[starting_idx + cidx])
                        if read_neigh_field and has_neighbors_field:
                            particle.nearest_neighbors = tipify(line[starting_idx + nidx])
                        for field_name, field_idx in zip(fields_to_read, fields_to_read_idx):
                            val = tipify(line[starting_idx + field_idx])
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
                for n in range(n_particles):
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
                    particle._index = n
                    system.particle.append(particle)

                # Add system to trajectory
                self._systems.append(system)

    def _parser_atooms(self):
        try:
            from atooms.trajectory import Trajectory as _Trajectory
        except ModuleNotFoundError:
            raise ModuleNotFoundError('No `atooms` module found.')

        # Fill the native Trajectory using atooms Trajectory
        atooms_traj = _Trajectory(self.filename, fmt=self.fmt)
        for atooms_sys in atooms_traj:
            cell = Cell(atooms_sys.cell.side)
            system = System(cell=cell)
            for n, atooms_p in enumerate(atooms_sys.particle):
                pos = atooms_p.position.copy()
                spe = atooms_p.species
                particle = Particle(position=pos, species=spe)
                particle._index = n
                # additional fields
                for field in self.additional_fields:
                    if field in self._default_neighbors_fields:
                        value = atooms_p.__getattribute__('neighbors')
                        particle.__setattr__('nearest_neighbors', value)
                    else:
                        value = atooms_p.__getattribute__(field)
                        particle.__setattr__(field, value)
                system.particle.append(particle)
            self._systems.append(system)
        atooms_traj.close()

    def _parser_mdtraj(self):
        try:
            import mdtraj as md
        except ModuleNotFoundError:
            raise ModuleNotFoundError('No `mdtraj` module found.')
        
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
                particle._index = atom
                system.particle.append(particle)
            self._systems.append(system)

    def _write_xyz(self, output_path, additional_fields, precision):
        from collections.abc import Iterable

        if not output_path.endswith('.xyz'):
            output_path += '.xyz'
        # deal with additional fields that are aliases
        processed_fields = []
        for field in additional_fields:
            if field in aliases.keys():
                processed_fields.append(aliases[field].split('.')[-1])
            else:
                processed_fields.append(field)
        with open(output_path, 'w') as file:
            for system in self._systems:
                file.write('{}\n'.format(len(system.particle)))
                columns = 'columns:id,pos,'
                for field in processed_fields:
                    columns += '{},'.format(field)
                columns = columns[:-1]
                header = columns + ' cell:{}\n'
                file.write(header.format(','.join('{:.{}f}'.format(L, precision) for L in system.cell.side)))
                for particle in system.particle:
                    line = '{} '.format(particle.species)
                    line += '{} '.format(' '.join('{:.{}f}'.format(p_i, precision) for p_i in particle.position))
                    for field in processed_fields:
                        attribute = particle.__getattribute__(field)
                        if isinstance(attribute, Iterable):
                            line += '{} '.format(','.join(map(str, attribute)))
                        else:
                            line += '{} '.format(attribute)
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
                    line = '{} '.format(particle.species_id - 1)
                    line += '{} '.format(' '.join('{:.{}f}'.format(p_i, precision) for p_i in particle.position))
                    # no additional field
                    # line += '{} '.format(particle.label)
                    line += '\n'
                    file.write(line)

    def _write_atooms(self, output_path, fmt, fields, precision):
        try:
            from atooms.trajectory import Trajectory as _Trajectory
            from atooms.system import System as _System
            from atooms.system import Particle as _Particle
            from atooms.system import Cell as _Cell
        except ModuleNotFoundError:
            print('No `atooms` module found.')

        # open an output trajectory
        atooms_traj = _Trajectory(output_path, fmt=fmt, mode='w')

        # additional fields
        variables = ['particle.species', 'particle.position']
        for field in fields:
            variables.append('particle.{}'.format(field))
        atooms_traj.variables = variables

        # write the output trajectory
        for n, system in enumerate(self._systems):
            new_cell = _Cell(side=system.cell.side)
            new_system = _System(cell=new_cell)
            for particle in system.particle:
                pos = particle.position
                spe = particle.species
                label = particle.label
                new_particle = _Particle(species=spe, position=pos)
                # additional fields
                for field in fields:
                    if field in self._default_neighbors_fields:
                        value = particle.__getattribute__('nearest_neighbors')
                    else:
                        value = particle.__getattribute__(field)
                    new_particle.__setattr__(field, value)
                new_particle.label = label
                new_system.particle.append(new_particle)
            atooms_traj.write(new_system, step=n)
        atooms_traj.close()

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
        frames = all_frames[first:last + 1:step]
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

    def __getitem__(self, item):
        return self._systems[item]

    def __len__(self):
        return len(self._systems)

    def __str__(self):
        rep = 'Trajectory(filename="{}", number_of_frames={})'
        return rep.format(self.filename, self.__len__())

    def __repr__(self):
        return self.__str__()