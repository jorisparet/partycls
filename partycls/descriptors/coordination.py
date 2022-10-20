import numpy
from .descriptor import StructuralDescriptor

class GlobalCoordinationDescriptor(StructuralDescriptor):
    """
    Global coordination descriptor.
    
    The global coordination number :math:`n_g(i)` of a particle :math:`i` is given by
    the number of its nearest neighbors that belong to ``group=1`` when computing
    the descriptor.
    
    Essentially, if no filter is applied on ``group=1``, :math:`n_g(i)` is the total
    number of neighbors of particle :math:`i`.

    The global coordination number can therefore be computed by considering an
    **arbitrary** subset of particles (*e.g.* particles whose radius is smaller
    than a certain value).

    The resulting feature vector for particle :math:`i` is given by

    .. math::
        X^\mathrm{GC}(i) = (\: n_g(i) \:) .
    
    .. note::
        Unlike most descriptors, this descriptor is **scalar**. Its feature vector
        :math:`X^\mathrm{GC}(i)` is thus composed of a single feature, and the 
        inherited ``grid`` attribute is therefore not relevant.

    See the tutorials for more details.

    Attributes
    ----------
    trajectory : Trajectory
        Trajectory on which the structural descriptor will be computed.
        
    active_filters : list
        All the active filters on both groups prior to the computation of the
        descriptor.
        
    dimension : int
        Spatial dimension of the descriptor (2 or 3).
        
    features : numpy.ndarray
        Array of all the structural features for the particles in group=0 in
        accordance with the defined filters (if any). This attribute is 
        initialized when the method ``compute`` is called (default value is ``None``).

    groups : tuple
        Composition of the groups: ``groups[0]`` and ``groups[1]`` contain lists of all
        the ``Particle`` instances in groups 0 and 1 respectively. Each element of 
        the tuple is a list of ``Particle`` in ``trajectory``, *e.g.* ``groups[0][0]``
        is the list of all the particles in the first frame of ``trajectory`` that 
        belong to group=0.

    verbose : bool
        Show progress information and warnings about the computation of the 
        descriptor when verbose is ``True``, and remain silent when verbose is 
        ``False``.

    neighbors_boost : float, default: 1.5
        Scaling factor to estimate the number of neighbors relative to a
        an ideal gas with the same density. This is used internally to set
        the dimensions of lists of neighbors. A too small number creates a
        risk of overfilling the lists of neighbors, and a too large number
        increases memory usage. This only works if the associated ``Trajectory``
        has valid cutoffs in the ``Trajectory.nearest_neighbors_cutoffs`` list
        attribute. This sets the value of the ``max_num_neighbors`` attribute
        during the computation of the descriptor.

    max_num_neighbors : int, default: 100
        Maximum number of neighbors. This is used internally to set the dimensions
        of lists of neighbors. This number is automatically adjusted to limit
        memory usage if the associated ``Trajectory`` has valid cutoffs in the 
        ``Trajectory.nearest_neighbors_cutoffs`` list attribute. The
        default value ``100`` is used if no cutoffs can be used to estimate a
        better value. The default value is sufficient in most cases, otherwise 
        this number can manually be increased **before** computing the descriptor.
    """

    name = 'global coordination'
    symbol = 'glco'

    def __init__(self, trajectory, accept_nans=True, verbose=False):
        """
        Parameters
        ----------
        trajectory : Trajectory
            Trajectory on which the structural descriptor will be computed.

        accept_nans: bool, default: True
            If ``False``, discard any row from the array of features that contains a 
            `NaN` element. If ``True``, keep `NaN` elements in the array of features.

        verbose : bool, default: False
            Show progress information and warnings about the computation of the 
            descriptor when verbose is ``True``, and remain silent when verbose 
            is ``False``.
        """
        StructuralDescriptor.__init__(self, trajectory,
                                      accept_nans=accept_nans,
                                      verbose=verbose)
        self.grid = numpy.zeros(1, dtype=numpy.float64)

    def compute(self):
        """
        Compute the global coordination number for the particles in group=0.
        Returns the data matrix and also overwrites the ``features`` attribute.

        Returns
        -------
        features : numpy.ndarray
            Data matrix with the global coordination number.
        """
        self._set_up(dtype=numpy.int64)
        self._manage_nearest_neighbors()
        self._filter_neighbors()
        n_frames = len(self.trajectory)
        # computation
        start = 0
        for n in self._trange(n_frames):
            npart = len(self.groups[0][n])
            feat_n = self._neighbors_number[n].copy()
            self.features[start: start+npart, 0] = feat_n
            start += npart
        self._handle_nans()
        return self.features


class ChemicalCoordinationDescriptor(StructuralDescriptor):
    """
    Chemical coordination descriptor.
    
    The chemical coordination number :math:`n_\\alpha(i)` of a particle :math:`i` 
    is given by the number of its nearest neighbors in ``group=1`` whose chemical
    species is :math:`\\alpha`.

    This is repeated for each chemical species :math:`\\alpha_{i=1 \dots n}` in the
    trajectory. The resulting feature vector for particle :math:`i` is given by

    .. math::
        X^\mathrm{CC}(i) = (\: n_{\\alpha_1}(i) \;\; \dots \;\; n_{\\alpha_n}(i) \:) .
    
    See the tutorials for more details.

    Attributes
    ----------
    trajectory : Trajectory
        Trajectory on which the structural descriptor will be computed.
        
    active_filters : list
        All the active filters on both groups prior to the computation of the
        descriptor.
        
    dimension : int
        Spatial dimension of the descriptor (2 or 3).

    grid : numpy.ndarray
        Grid of chemical species :math:`\{ \\alpha_i \}.`
        
    features : numpy.ndarray
        Array of all the structural features for the particles in group=0 in
        accordance with the defined filters (if any). This attribute is 
        initialized when the method ``compute`` is called (default value is ``None``).

    groups : tuple
        Composition of the groups: ``groups[0]`` and ``groups[1]`` contain lists of all
        the ``Particle`` instances in groups 0 and 1 respectively. Each element of 
        the tuple is a list of ``Particle`` in ``trajectory``, *e.g.* ``groups[0][0]``
        is the list of all the particles in the first frame of ``trajectory`` that 
        belong to group=0.

    verbose : bool
        Show progress information and warnings about the computation of the 
        descriptor when verbose is ``True``, and remain silent when verbose is 
        ``False``.

    neighbors_boost : float, default: 1.5
        Scaling factor to estimate the number of neighbors relative to a
        an ideal gas with the same density. This is used internally to set
        the dimensions of lists of neighbors. A too small number creates a
        risk of overfilling the lists of neighbors, and a too large number
        increases memory usage. This only works if the associated ``Trajectory``
        has valid cutoffs in the ``Trajectory.nearest_neighbors_cutoffs`` list
        attribute. This sets the value of the ``max_num_neighbors`` attribute
        during the computation of the descriptor.

    max_num_neighbors : int, default: 100
        Maximum number of neighbors. This is used internally to set the dimensions
        of lists of neighbors. This number is automatically adjusted to limit
        memory usage if the associated ``Trajectory`` has valid cutoffs in the 
        ``Trajectory.nearest_neighbors_cutoffs`` list attribute. The
        default value ``100`` is used if no cutoffs can be used to estimate a
        better value. The default value is sufficient in most cases, otherwise 
        this number can manually be increased **before** computing the descriptor.
    """

    name = 'chemical coordination'
    symbol = 'chco'

    def __init__(self, trajectory, accept_nans=True, verbose=False):
        """
        Parameters
        ----------
        trajectory : Trajectory
            Trajectory on which the structural descriptor will be computed.

        accept_nans: bool, default: True
            If ``False``, discard any row from the array of features that contains a 
            `NaN` element. If ``True``, keep `NaN` elements in the array of features.

        verbose : bool, default: False
            Show progress information and warnings about the computation of the 
            descriptor when verbose is ``True``, and remain silent when verbose 
            is ``False``.
        """
        StructuralDescriptor.__init__(self, trajectory,
                                      accept_nans=accept_nans,
                                      verbose=verbose)
        self.grid = self.trajectory[0].distinct_species

    def compute(self):
        """
        Compute the chemical coordination numbers for the particles in group=0.
        Returns the data matrix and also overwrites the ``features`` attribute.

        Returns
        -------
        features : numpy.ndarray
            Data matrix with the chemical coordination number.
        """
        self._set_up(dtype=numpy.int64)
        self._manage_nearest_neighbors()
        self._filter_neighbors()
        n_frames = len(self.trajectory)
        row = 0
        # all relevant arrays
        spe_all = self.trajectory.dump('species')
        nn_spe_i = numpy.empty_like(self.grid, dtype=numpy.int64)
        # computation
        for n in self._trange(n_frames):
            for i in range(len(self.groups[0][n])):
                nn_i = self._neighbors_number[n][i]
                neigh_i = self._neighbors[n][i,0:nn_i]
                for j, spe_j in enumerate(self.grid):
                    nn_spe_i[j] = numpy.count_nonzero(spe_all[n][neigh_i] == spe_j)
                self.features[row] = nn_spe_i
                row += 1
        self._handle_nans()
        return self.features