import numpy
from .descriptor import AngularStructuralDescriptor
from .realspace_wrap import compute

class CompactnessDescriptor(AngularStructuralDescriptor):
    """
    Compactness descriptor. Structural descriptor based on the compactness metric 
    (or packing efficiency) defined by Tong & Tanaka 
    (https://doi.org/10.1103/PhysRevX.8.011041).

    Measures the compactness of the local environment around a central particle
    and computes its deviation from an ideal packing configuration.

    This descriptor is scalar. Therefore, the `grid` attribute is not relevant.
    
    See the parent class for more details.

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
        Grid over which the structural features will be computed.
        
    features : numpy.ndarray
        Array of all the structural features for the particles in group=0 in
        accordance with the defined filters (if any). This attribute is 
        initialized when the method `compute` is called (default value is None).
    """

    name = 'compactness'
    symbol = 'compact'
    
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
        AngularStructuralDescriptor.__init__(self,
                                             trajectory,
                                             accept_nans=accept_nans,
                                             verbose=verbose)
        self.grid = numpy.zeros(1, dtype=numpy.float64)
        
    def compute(self):
        # set up
        self._set_up(dtype=numpy.float64)
        self._manage_nearest_neighbors()
        self._filter_neighbors()
        self._filter_subsidiary_neighbors()
        n_frames = len(self.trajectory)
        row = 0
        # all relevant arrays
        pos_all = self.trajectory.dump('position')
        radii = self.trajectory.dump('radius')
        box = self.trajectory.dump('cell.side')
        # computation
        for n in self._trange(n_frames):
            pos_all_n = pos_all[n].T
            for i in range(len(self.groups[0][n])):
                nn_i = self._neighbors_number[n][i]
                tetra_i = self.tetrahedra(i, 
                                          list(self._neighbors[n][i][0:nn_i]),
                                          self._subsidiary_neighbors[n][i])
                theta_i = compute.compactness(pos_all_n, tetra_i.T, radii[n], box[n])
                self.features[row] = theta_i
                row += 1
        self._handle_nans()
        return self.features
    
    def tetrahedra(self, i, neigh_i, neigh_neigh_i):
        """
        Return a 2D array that contains all the tetrahedra around
        particle `i`. Each row is a 4-array with the indices of the
        particles forming the tetrahedron. The first index is `i`
        by construction.

        Parameters
        ----------
        i : int
            Index of the central particle.

        neigh_i : list
            List of neighbors of the central particles.

        neigh_neigh_i: list
            Neighbors of the neighbors of the central particle.
            Each element is a list.

        Returns
        -------
        tetra_i : array
            The array that contains the indices of all the tetrahedra
            around particle `i`.
        """
        tetras = []
        for j_idx, j in enumerate(neigh_i):
            neigh_j = neigh_neigh_i[j_idx]
            common = set(neigh_i) & set(neigh_j)
            for k in common:
                for m in common:
                    if m <= k:
                        continue
                    k_idx = neigh_i.index(k)
                    if m in neigh_neigh_i[k_idx]:
                        tetras.append([i,j,k,m])
        return numpy.array(tetras, dtype=numpy.int64)
