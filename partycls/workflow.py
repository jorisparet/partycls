"""
Workflow for clustering analysis.

A workflow is a procedure that goes through various steps (some of which are 
optional) to perform a structural clustering on a trajectory.
"""

from .trajectory import Trajectory
from partycls.descriptors import *
from .clustering import KMeans, GaussianMixture, CommunityInference
from .dim_reduction import PCA, TSNE, LocallyLinearEmbedding, AutoEncoder
from .feature_scaling import ZScore, MinMax, MaxAbs, Robust
from partycls.core import __version__ as _version
from partycls.core import __code_extension__ as code_extension
import time
import datetime


class Workflow:
    """
    A Workflow is a clustering procedure that goes through the following steps:

    - compute a structural descriptor on a given trajectory ;
    - (*optional*) apply a feature scaling on the previously computed structural features ;
    - (*optional*) apply a dimensionality reduction on the (raw/scaled) features ;
    - run a clustering algorithm to partition particles into structurally different clusters ;
    
    Attributes
    ----------
    trajectory : Trajectory
        The trajectory file as read by the ``Trajectory`` class.
        
    descriptor : StructuralDescriptor
        Structural descriptor associated to the trajectory.
        
    scaling : ZScore, MinMax, MaxAbs or Robust
        Feature scaling method.
        
    dim_reduction : PCA, TSNE, LocallyLinearEmbedding or AutoEncoder
        Dimensionality reduction method.
        
    clustering : Clustering
        Clustering method.
        
    output_metadata : dict
        Dictionnary that controls the writing process and 
        the properties of all the output files.
        
    features: numpy.ndarray
        Raw features as computed by the associated structural descriptor.
        Initial value is ``None`` if features were not computed.
        
    scaled_features: numpy.ndarray
        Features after being rescaled by a feature scaling method.
        Equal to ``None`` if no scaling is applied to the features.
        
    reduced_features: numpy.ndarray
        Features in the reduced space after applying a dimensionality reduction
        technique. Equal to ``None`` if no reduction is applied to the features.
        
    naming_convention : str
        Base name for output files. 
        Default is ``"{filename}.{code}.{descriptor}.{clustering}"``, where each
        tag will be replaced by its value in the current instance of 
        ``Workflow`` (*e.g.* ``"traj.xyz.partycls.gr.kmeans"``).
        
        Base name can be changed using any combination of the available tags:

        - ``{filename}``
        - ``{code}``
        - ``{descriptor}``
        - ``{scaling}``
        - ``{dim_reduction}``
        - ``{clustering}``
    
        Example: ``"{filename}_descriptor-{descriptor}_scaling-{scaling}.{code}"``.
    """

    # Databases
    descriptor_db = {RadialDescriptor.symbol: RadialDescriptor,
                     'rad': RadialDescriptor,
                     'radial': RadialDescriptor,
                     'gr': RadialDescriptor,
                     BondAngleDescriptor.symbol: BondAngleDescriptor,
                     'ang': BondAngleDescriptor,
                     'angular': BondAngleDescriptor,
                     SmoothedBondOrientationalDescriptor.symbol: SmoothedBondOrientationalDescriptor,
                     BondOrientationalDescriptor.symbol: BondOrientationalDescriptor,
                     'boo': BondOrientationalDescriptor,
                     'bop': BondOrientationalDescriptor,
                     'steinhardt': SteinhardtDescriptor,
                     LocallyAveragedBondOrientationalDescriptor.symbol: LocallyAveragedBondOrientationalDescriptor,
                     'ld': LechnerDellagoDescriptor,
                     'lechner-dellago': LechnerDellagoDescriptor,
                     'lechner dellago': LechnerDellagoDescriptor,
                     SmoothedBondOrientationalDescriptor.symbol: SmoothedBondOrientationalDescriptor,
                     'sboo': SmoothedBondOrientationalDescriptor,
                     'sbop': SmoothedBondOrientationalDescriptor,
                     RadialBondOrientationalDescriptor.symbol: RadialBondOrientationalDescriptor,
                     'rboo': RadialBondOrientationalDescriptor,
                     'rbop': RadialBondOrientationalDescriptor,
                     'boattini': BoattiniDescriptor,
                     TetrahedralDescriptor.symbol: TetrahedralDescriptor,
                     CompactnessDescriptor.symbol: CompactnessDescriptor,
                     'tong-tanaka': TongTanakaDescriptor,
                     'tong tanaka': TongTanakaDescriptor,
                     CoordinationDescriptor.symbol: CoordinationDescriptor,
                     'coordination': CoordinationDescriptor}

    clustering_db = {'k-means': KMeans,
                     'kmeans': KMeans,
                     'gaussian mixture': GaussianMixture,
                     'gaussian-mixture': GaussianMixture,
                     'gmm': GaussianMixture,
                     'gm': GaussianMixture,
                     'community inference': CommunityInference,
                     'community-inference': CommunityInference,
                     'inference': CommunityInference,
                     'cinf': CommunityInference}

    scaling_db = {'standard': ZScore,
                  'zscore': ZScore,
                  'z-score': ZScore,
                  'minmax': MinMax,
                  'min-max': MinMax,
                  'maxabs': MaxAbs,
                  'max-abs': MaxAbs,
                  'robust': Robust}

    dim_reduction_db = {'pca': PCA,
                        'tsne': TSNE,
                        't-sne': TSNE,
                        'lle': LocallyLinearEmbedding,
                        'autoencoder': AutoEncoder,
                        'auto-encoder': AutoEncoder,
                        'ae': AutoEncoder}

    def __init__(self, trajectory, descriptor='gr', scaling=None, dim_reduction=None, clustering='kmeans'):
        """
        Parameters
        ----------
        trajectory : Trajectory
            An instance of ``Trajectory`` or a path to trajectory file to read, or 
            an instance of a class with compatible interface.
            
        descriptor : StructuralDescriptor, default: "gr"
            An instance of ``StructuralDescriptor``, the short name of a descriptor 
            (str), or an instance of a class with compatible interface. See the
            ``descriptor_db`` class attribute for compatible strings. Examples:
            
            - ``"gr"`` : radial distribution of particles around a central particle.
            - ``"ba"`` : angular distribution of pairs of nearest neighbors of a central particle.
            - ``"bo"`` : Steinhardt bond-orientational order parameter.
            - ``"ld"`` : Lechner-Dellago cond-orientational order parameter.
            
        scaling : method, default: None
            Feature scaling method. See the ``scaling_db`` class attribute for
            compatible strings. Examples:
            
            - ``"zscore"`` : standardize features by removing the mean and scaling to unit variance
            - ``"minmax"`` : scale and translate each feature individually such that it is in the given range on the training set, *e.g.* between zero and one
            - ``"maxabs"`` : scale and translate each feature individually such that the maximal absolute value of each feature in the training set will be 1.
            - ``"robust"`` : remove the median and scale the data according to the specified quantile range (default is between 25th quantile and 75th quantile)
                
        dim_reduction : method, default: None
            Dimensionality reduction method. See the ``dim_reduction_db`` class attribute
            for compatible strings. Examples:
            
            - ``"pca"`` : Principal Component Analysis
            - ``"tsne"`` : t-distributed Stochastic Neighbor Embedding
            - ``"lle"`` : Locally Linear Embedding
            - ``"ae"`` : neural network Auto-Encoder
            
        clustering : Clustering, default: 'kmeans'
            Clustering algorithm. See the ``clustering_db`` class attribute for 
            compatible strings. Examples:
            
            - ``"kmeans"`` : K-Means algorithm
            - ``"gmm"`` : Gaussian Mixture Model
            - ``"cinf"`` : Community Inference (see https://doi.org/10.1063/5.0004732).

        Example
        -------
        >>> wf = Workflow('trajectory.xyz', descriptor='ba', scaling='zscore')
        >>> wf.run()
        """
        # Trajectory
        if isinstance(trajectory, str):
            self.trajectory = Trajectory(trajectory)
        else:
            self.trajectory = trajectory

        # Descriptor
        if isinstance(descriptor, str):
            self.descriptor = self.descriptor_db[descriptor.lower()](self.trajectory)
        else:
            self.descriptor = descriptor
        self.features = self.descriptor.features

        # Feature scaling
        if isinstance(scaling, str):
            self.scaling = self.scaling_db[scaling.lower()]()
        else:
            self.scaling = scaling
        self.scaled_features = None

        # Dimensionality reduction
        if isinstance(dim_reduction, str):
            self.dim_reduction = self.dim_reduction_db[dim_reduction.lower()]()
        else:
            self.dim_reduction = dim_reduction
        self.reduced_features = None

        # Clustering
        if isinstance(clustering, str):
            self.clustering = self.clustering_db[clustering.lower()]()
        else:
            self.clustering = clustering

        # Default output metadata
        self.output_metadata = {'trajectory': {'enable': True,
                                               'writer': self.write_trajectory,
                                               'filename': None,
                                               'fmt': 'xyz',
                                               'backend': None,
                                               'additional_fields': [],
                                               'precision': 6},

                                'log': {'enable': True,
                                        'writer': self.write_log,
                                        'filename': None,
                                        'precision': 6},

                                'centroids': {'enable': True,
                                              'writer': self.write_centroids,
                                              'filename': None,
                                              'precision': 6},

                                'labels': {'enable': False,
                                           'writer': self.write_labels,
                                           'filename': None},

                                'dataset': {'enable': False,
                                            'writer': self.write_dataset,
                                            'filename': None,
                                            'precision': 6}}

        # Naming convention for output files
        self.naming_convention = '{filename}.{code}.{descriptor}.{clustering}'

        # Internal
        self._has_run = False
        self._start = None
        self._end = None
        self._time = None

    @property
    def labels(self):
        """
        Clustering labels.
        """
        return self.clustering.labels

    @property
    def fractions(self):
        """
        Fraction of particles in each cluster.
        """
        return self.clustering.fractions

    @property
    def populations(self):
        """
        Number of particles in each cluster.
        """
        return self.clustering.populations

    @property
    def centroids(self):
        """
        Centroid of each cluster.
        """
        return self.clustering.centroids

    # TODO: like scipy functions, we may return a dict() with the optimization results
    def run(self):
        """
        Compute the clustering and write the output files according to the
        defined Workflow :

        - compute the descriptor
        - (*optional*) apply feature scaling
        - (*optional*) apply dimensionality reduction
        - compute the clustering
        - (*optional*) write the output files     

        Raises
        ------
        ValueError
            If a community inference clustering is attempted with feature
            scaling or dimensionality reduction.

        Returns
        -------
        None
        """
        # Start the timer
        self._start = time.time()

        # Make sure the descriptor has been computed
        if self.descriptor.features is None:
            self.features = self.descriptor.compute()

        # Feature scaling
        if self.scaling is None:
            X_scaled = self.descriptor.features.copy()
        else:
            X_scaled = self.scaling.fit_transform(self.descriptor.features)
            self.scaled_features = X_scaled

        # Dimensionality reduction
        if self.dim_reduction is None:
            X_red = X_scaled.copy()
        else:
            X_red = self.dim_reduction.reduce(X_scaled)
            self.reduced_features = X_red

        # Clustering
        #  community inference needs to fit an instance of descriptor
        if self.clustering.symbol == 'cinf' and (self.scaling is not None or self.dim_reduction is not None):
            raise ValueError('community inference is not meant to run with feature-scaling or dimensionality reduction')
        elif self.clustering.symbol == 'cinf':
            self.clustering.fit(self.descriptor)
        else:
            self.clustering.fit(X_red)
        #  give its predicted label to each selected `Particle` in the trajectory.
        n = 0
        for frame in self.descriptor.groups[0]:
            for particle in frame:
                particle.label = self.labels[n]
                n += 1

        # Workflow has run at least once
        self._has_run = True
        # End the timer
        self._end = time.time()
        self._time = self._end - self._start

        # Outputs
        for filetype in self.output_metadata.keys():
            enable = self.output_metadata[filetype]['enable']
            if enable:
                writer = self.output_metadata[filetype]['writer']
                writer(**self.output_metadata[filetype])

    def set_output_metadata(self, what, **kwargs):
        """
        Change the output properties.

        Parameters
        ----------
        what : str
            Type of output file to change. Must be one of:
            
            - ``"trajectory"``
            - ``"log"``
            - ``"centroids"``
            - ``"labels"``
            - ``"dataset"``

        **kwargs : 
            Keywords arguments (specific to each type of file)

        Returns
        -------
        None
        
        Examples
        --------
        >>> wf = Workflow('trajectory.xyz')
        >>> wf.set_output_metadata('log', enable=False) # do not write the log file
        >>> wf.set_output_metadata('trajectory', filename='awesome_trajectory.xyz') # change the default output name
        >>> wf.run('dataset', enable=True, precision=8) # write the dataset and change the writing precision to 8 digits
        """
        for key, val in kwargs.items():
            self.output_metadata[what][key] = val

    def disable_output(self):
        """
        Disable all outputs.

        Returns
        -------
        None
        """
        for key in self.output_metadata.keys():
            self.output_metadata[key]['enable'] = False

    def write_trajectory(self, filename=None, fmt='xyz', backend=None, additional_fields=None, precision=6, **kwargs):
        """
        Write the trajectory file with cluster labels (default) and other
        additional fields (if any).        

        Parameters
        ----------
        filename : str, default: None
            Filename of the output trajectory. Uses a default naming convention
            if not specified. The default is None.

        fmt : str, default: "xyz"
            Output trajectory format.

        backend : str, default: None
            Name of the backend to use to write the trajectory. Must be either
            ``None``, ``"atooms"`` or ``"mdtraj"``.

        additional_fields : list, default: None
            Additional fields (*i.e.* particle properties) to write in the
            output trajectory. Note that all the ``Particle`` objects should
            have the specified properties as attributes.

        precision : int, default: 6
            Number of decimals when writing the output trajectory.

        Returns
        -------
        None

        Examples
        --------
        >>> wf = Workflow('trajectory.xyz')
        >>> wf.write_trajectory(fmt='rumd')
        >>> wf.write_trajectory(additional_field=['particle.mass']) # `Particle` must have the attribute `mass`.
        >>> wf.write_trajectory(filename='my_custom_name', precision=8)
        """
        if additional_fields is None:
            additional_fields = []
        if filename is None:
            filename = self._output_file(fmt)
        self.trajectory.write(filename, fmt=fmt, backend=backend,
                              additional_fields=['label'] + additional_fields,
                              precision=precision)

    # TODO: more info needed in the log?
    # TODO: log should be a string, so that it can be parsed/printed by python code
    def write_log(self, filename=None, precision=6, **kwargs):
        """
        Write a log file with all relevant information about the workflow.
        The log file can be written only if the workflow has been run at
        least once with the method `Workflow.run`.

        Parameters
        ----------
        filename : str, default: None
            Filename of the log file. Uses a default naming convention
            if not specified.

        precision : int, default: 6
            Number of decimals when writing the log file.

        Returns
        -------
        None
        """
        if filename is None:
            filename = self._output_file('log')
        if self._has_run:
            filters = self.descriptor.active_filters
            fractions = self.clustering.fractions
            n_init = self.clustering.n_init
            with open(filename, 'w') as file:
                file.write('# title: workflow log \n')
                file.write(self._get_header())
                file.write('\nExecution time: {:.{}f}s \n'.format(self._time, precision))
                file.write('\nNumber of repetitions: {} \n'.format(n_init))
                # fractions
                file.write('\nFractions (k, f_k): \n')
                for k, f_k in enumerate(fractions):
                    file.write('- {} {:.{}f} \n'.format(k, f_k, precision))
                # filters
                if len(filters) == 0:
                    file.write('\nNo filter applied \n')
                else:
                    file.write('\nApplied filters: \n')
                    for fltr in filters:
                        file.write('- group: {}, filter: {} \n'.format(fltr[1], fltr[0]))

    def write_centroids(self, filename=None, precision=6, **kwargs):
        """
        Write the coordinates of the clusters' centroids using the raw features
        from the descriptor (i.e. nor scaled or reduced).

        Parameters
        ----------
        filename : str, default: None
            Filename of the centroids file. Uses a default naming convention
            if not specified.

        precision : int, default: 6
            Number of decimals when writing the centroids file.

        Returns
        -------
        None
        """
        if filename is None:
            filename = self._output_file('centroids')
        with open(filename, 'w') as file:
            kind = self.descriptor.name
            file.write('# title: centroids ({} features)\n'.format(kind))
            file.write(self._get_header())
            file.write('# columns: feature, centroids \n')
            grid = self.descriptor.grid
            C_k = self.clustering.centroids(self.descriptor.features)
            n_clusters = self.clustering.n_clusters
            for n, g_i, in enumerate(grid):
                g = '{:.{}f} '.format(g_i, precision)
                line = g + ''.join(['{:.{}f} '.format(C_k[k, n], precision) for k in range(n_clusters)]) + '\n'
                file.write(line)

    def write_labels(self, filename=None, **kwargs):
        """
        Write the clusters' labels only.

        Parameters
        ----------
        filename : str, default: None
            Filename of the labels file. Uses a default naming convention
            if not specified.

        Returns
        -------
        None
        """
        if filename is None:
            filename = self._output_file('labels')
        with open(filename, 'w') as file:
            file.write("# title: clusters' labels\n")
            file.write(self._get_header())
            for ki in self.labels:
                file.write('{} \n'.format(ki))

    def write_dataset(self, filename=None, precision=6, **kwargs):
        """
        Write the full raw dataset from the descriptor as an array (*i.e.* all 
        the individual raw features of each particle).

        Parameters
        ----------
        filename : str, default: None
            Filename of the dataset file. Uses a default naming convention
            if not specified.

        precision : int, default: 6
            Number of decimals when writing the dataset file.

        Returns
        -------
        None
        """
        if filename is None:
            filename = self._output_file('dataset')
        with open(filename, 'w') as file:
            file.write('# title: data set matrix \n')
            file.write(self._get_header())
            for vector in self.descriptor.features:
                line = ''.join(['{:.{}f} '.format(vi, precision) for vi in vector]) + '\n'
                file.write(line)

    def _get_header(self):
        # Time
        now = datetime.datetime.now()
        date = '# date: {Y}-{M:02}-{D:02} {h}:{m:02}:{s:02} \n'.format(Y=now.year,
                                                                       M=now.month,
                                                                       D=now.day,
                                                                       h=now.hour,
                                                                       m=now.minute,
                                                                       s=now.second)

        # Version
        version = '# version: {} \n'.format(_version)

        # Parent
        parent = '# parent: {} \n'.format(self.trajectory.filename)

        # Feature scaling
        scaling = '# feature scaling: {} \n'
        if self.scaling is not None:
            scaling = scaling.format(self.scaling.full_name)
        else:
            scaling = scaling.format('none')

        # Dimensionality reduction
        dim_reduction = '# dimensionality reduction: {} \n'
        if self.dim_reduction is not None:
            dim_reduction = dim_reduction.format(self.dim_reduction.full_name)
        else:
            dim_reduction = dim_reduction.format('none')

        # Clustering method
        clustering = '# clustering method: {} \n'.format(self.clustering.full_name)

        # Number of communities/clusters
        n_clusters = '# requested number of clusters: {} \n'.format(self.clustering.n_clusters)

        # Assemble header
        header = ''.join([date, version, parent, scaling, dim_reduction, clustering, n_clusters])
        return header

    def _output_file(self, kind):
        """Returns path of output file."""
        scaling_symbol = self.scaling.symbol if self.scaling is not None else ''
        dim_reduction_symbol = self.dim_reduction.symbol if self.dim_reduction is not None else ''
        base_name = self.naming_convention.format(filename=self.trajectory.filename,
                                                  code=code_extension,
                                                  descriptor=self.descriptor.symbol,
                                                  scaling=scaling_symbol,
                                                  dim_reduction=dim_reduction_symbol,
                                                  clustering=self.clustering.symbol)
        return '{base}.{kind}'.format(base=base_name, kind=kind)

    def __str__(self):
        rep = 'Workflow(filename="{}", descriptor="{}", scaling="{}", dim_reduction="{}", clustering="{}", has_run={})'
        rep = rep.format(self.trajectory.filename,
                         self.descriptor.symbol,
                         self.scaling.symbol if self.scaling is not None else None,
                         self.dim_reduction.symbol if self.dim_reduction is not None else None,
                         self.clustering.symbol,
                         self._has_run)
        return rep

    def __repr__(self):
        return self.__str__()
