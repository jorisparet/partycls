from pysc.trajectory import Trajectory
from pysc.descriptor import AngularDescriptor, RadialDescriptor, BondOrientationalDescriptor, LechnerDellagoDescriptor
from .clustering import KMeans, GaussianMixture, CommunityInference
from .dim_redux import PCA, TSNE, AutoEncoder
from .feature_scaling import ZScore, MinMax
import time, datetime

# Code name
_output_path = '{filename}.pysc.{mode}.{method}.{kind}'

# Code version
ci_version = '0.1'

# Databases
descriptor_db = {'gr': RadialDescriptor,
                 'ba': AngularDescriptor,
                 'bo': BondOrientationalDescriptor,
                 'ld': LechnerDellagoDescriptor}
scaling_db = {'zscore': ZScore,
              'minmax': MinMax}
dim_redux_db = {'pca': PCA,
                'tsne': TSNE,
                'ae': AutoEncoder}
clustering_db = {'kmeans': KMeans,
                 'gmm': GaussianMixture,
                 'cinf': CommunityInference}

class Optimization:
    """
    An optimization is a workflow that goes through the following steps:
    - compute a structural descriptor on a given trajectory ;
    - (optional) apply a feature scaling on the previously computed structural features ;
    - (optional) apply a dimensionality reduction on the (raw/scaled) features ;
    - run a clustering algorithm to partition particles into structurally different clusters ;
    
    Parameters
    ----------
    
    trajectory : Trajectory, or str
        An instance of `Trajectory` or path to trajectory file to read.
        
    descriptor : {'gr', 'ba', 'bo', 'ld', or an instance of StructuralDescriptor}
        Structural descriptor to be computed on the trajectory.
        
        'gr' : radial distribution of particles around a central particle.
        
        'ba' : angular distribution of pairs of nearest neighbors of a central particle.
        
        'bo' : Steinhardt bond-orientational order parameter (see https://doi.org/10.1103/PhysRevB.28.784)
        
        'ld' : Lechner-Dellago cond-orientational order parameter (see https://doi.org/10.1063/1.2977970)
        
    scaling : {'zscore', 'minmax', None, or an object with the proper interface}, optional, default: None
        Feature scaling method.
        
        'zscore' : standardize features by removing the mean and scaling to unit variance
            
        'minmax' : scale and translate each feature individually such that it 
        is in the given range on the training set, e.g. between zero and one
            
    dim_redux : {'pca', 'tsne', 'lle', 'ae', None, or an object with the proper interface}, optional, default: None
        Dimensionality reduction method.
        
        'pca' : Principal Component Analysis.
        
        'tsne' : t-distributed Stochastic Neighbor Embedding.
        
        'lle' : Locally Linear Embedding.
        
        'ae' : neural network Auto-Encoder.
        
    clustering : {'kmeans', 'gmm', 'cinf'}, default: 'kmeans'
        Clustering algorithm.
        
        'kmeans' : K-Means algorithm.
        
        'gmm' : Gaussian Mixture Model.
        
        'cinf' : Community Inference (see https://doi.org/10.1063/5.0004732).
    
    Attributes
    ----------
    
    trajectory : Trajectory
        The trajectory file as read by the Trajectory class.
        
    descriptor : StructuralDescriptor
        Structural descriptor associated to the trajectory.
        
    scaling : ZScore, MinMax
        Feature scaling.
        
    dim_redux : PCA, TSNE, LocallyLinearEmbedding, AutoEncoder
        Dimensionality reduction.
        
    clustering : Clustering
        Clustering method.
        
    output_metadata : dict
        Dictionnary that controls the writing process and 
        the properties of all the output files.
        
    Examples
    --------
    
    >>> from pysc import Optimization
    >>> opt = Optimization('trajectory.xyz', descriptor='ba', scaling='zscore')
    >>> opt.run()
    """
    
    def __init__(self, trajectory, descriptor='gr', scaling=None, dim_redux=None, clustering='kmeans'):
        
        # Trajectory
        if isinstance(trajectory, Trajectory):
            self.trajectory = trajectory
        elif isinstance(trajectory, str):
            self.trajectory = Trajectory(trajectory)
        else:
            raise TypeError('`trajectory` should be an instance of `str` or `Trajectory`.')

        # Descriptor
        if isinstance(descriptor, str):
            self.descriptor = descriptor_db[descriptor](self.trajectory)
        else:
            self.descriptor = descriptor

        # Feature scaling
        if isinstance(scaling, str):
            self.scaling = scaling_db[scaling]()
        else:
            self.scaling = scaling

        # Dimensionality reduction
        if isinstance(dim_redux, str):
            self.dim_redux = dim_redux_db[dim_redux]()
        else:
            self.dim_redux = dim_redux
            
        # Clustering
        if isinstance(clustering, str):
            self.clustering = clustering_db[clustering]()
        else:
            self.clustering = clustering
            
        # Default output metadata
        self.output_metadata = {'trajectory': {'enable':True,
                                               'writer':self.write_trajectory,
                                               'filename':None,
                                               'fmt':'xyz',
                                               'additional_fields':[],
                                               'precision':6},
    
                                'log': {'enable':True,
                                        'writer':self.write_log,
                                        'filename':None,
                                        'precision':6},
                                        
                                'centroids': {'enable':True,
                                              'writer':self.write_centroids,
                                              'filename':None,
                                              'precision':6},
                                
                                'labels': {'enable':False,
                                           'writer':self.write_labels,
                                           'filename':None},
                                           
                                'dataset': {'enable':False,
                                            'writer':self.write_dataset,
                                            'filename':None,
                                            'precision':6}}
                                
        # Internal
        self._has_run = False
        self._start = None
        self._end = None
        self._time = None
        
    def run(self):
        """
        Compute the clustering and write the output files according to the
        defined workflow :
        - compute the descriptor ;
        - (optional) apply feature scaling ;
        - (optional) apply dimensionality reduction ;
        - compute the clustering ;
        - (optional) write the output files ;
        """
        # Start the timer
        self._start = time.time()

        # Make sure the descriptor has been computed
        if self.descriptor.features is None:
            self.descriptor.compute()

        # Feature scaling
        if self.scaling is None:
            X_scaled = self.descriptor.features.copy()
        else:
            X_scaled = self.scaling.fit_transform(self.descriptor.features)

        # Dimensionality reduction
        if self.dim_redux is None:
            X_red = X_scaled.copy()
        else:
            X_red = self.dim_redux.reduce(X_scaled)

        # Clustering
        if self.clustering.symbol == 'cinf':
            # TODO: fix that. Not pretty and not practical!
            print('Warning: community inference is using the whole descriptor.')
            self.clustering.fit(self.descriptor)
        else:
            self.clustering.fit(X_red)
        #  give its predicted label to each selected `Particle` in the trajectory.
        n = 0
        for frame in self.descriptor._groups[0]:
            for particle in frame:
                particle.label = self.labels[n]
                n += 1
                
        # Optimization has run at least once
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

    @property
    def labels(self):
        return self.clustering.labels
    
    @property
    def fractions(self):
        return self.clustering.fractions
    
    @property
    def populations(self):
        return self.clustering.populations
    
    @property
    def centroids(self):
        return self.clustering.centroids

    def set_output_metadata(self, what, **kwargs):
        """
        Change the output properties.
        
        Parameters
        ----------
        what : {'trajectory', 'log', 'centroids', 'labels', or 'dataset'}
            Type of output file to change.
            
        kwargs : keywords arguments (specific to each type of file)
        
        Examples
        --------
        >>> opt = Optimisation('trajectory.xyz')
        >>> opt.set_output_metadata('log', enable=False) # do not write the log file
        >>> opt.set_output_metadata('trajectory', filename='awesome_trajectory.xyz') # change the default output name
        >>> opt.run('dataset', enable=True, precision=8) # write the dataset and change the writing precision to 8 digits
        """
        for key, val in kwargs.items():
            self.output_metadata[what][key] = val
        
    def disable_output(self):
        """
        Disable all outputs.
        """
        for key in self.output_metadata.keys():
            self.output_metadata[key]['enable'] = False

    def write_trajectory(self, filename=None, fmt='xyz', additional_fields=[], precision=6, **kwargs):
        """
        Write the trajectory file with cluster labels (default) and other
        additional fields (if any).
        
        Parameters
        ----------
        filename : str, optional, default: None
            Filename of the output trajectory. Uses a default naming convention
            if not specified.
            
        fmt : {'xyz', 'rumd'}, default: 'xyz'
            Output trajectory format.
            
        additional_fields : list of str, optional, default: []
            Additional fields (i.e. particle properties) to write in the
            output trajectory. Note that all the `Particle` objects should
            have the specified properties as attributes.
            
        precision : int, optional, default: 6
            Number of decimals when writing the output trajectory.
            
        Examples
        --------
        >>> opt = Optimisation('trajectory.xyz')
        >>> opt.write_trajectory(fmt='rumd')
        >>> opt.write_trajectory(additional_field=['particle.mass']) # `Particle` must have the attribute `mass`.
        >>> opt.write_trajectory(filename='my_custom_name', precision=8)
        """
        if filename is None:
            filename = self._output_file(fmt)
        self.trajectory._write(filename, fmt=fmt, 
                               additional_fields=['label']+additional_fields, 
                               precision=precision)

    # TODO: more info needed in the log?
    def write_log(self, filename=None, precision=6, **kwargs):
        """
        Write a log file with all relevant information about the optimization.
        The log file can be written only if the optimization has been run at
        least once with the method `Optimization.run`.
        
        Parameters
        ----------
        filename : str, optional, default: None
            Filename of the log file. Uses a default naming convention
            if not specified.   
            
        precision : int, optional, default: 6
            Number of decimals when writing the log file. 
        """
        if filename is None:
            filename = self._output_file('log')
        if self._has_run:
            filters = self.descriptor.active_filters
            fractions = self.clustering.fractions
            n_init = self.clustering.n_init
            with open(filename, 'w') as file:
                file.write('# title: optimization log \n')
                file.write(self._get_header())
                file.write('\nOptimization time: {:.{}f}s \n'.format(self._time, precision))
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
        filename : str, optional, default: None
            Filename of the centroids file. Uses a default naming convention
            if not specified.   
            
        precision : int, optional, default: 6
            Number of decimals when writing the centroids file. 
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
        filename : str, optional, default: None
            Filename of the labels file. Uses a default naming convention
            if not specified.   
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
        Write the full raw dataset from the descriptor as an array (i.e. all 
        the individual raw features of each particle).
        
        Parameters
        ----------
        filename : str, optional, default: None
            Filename of the dataset file. Uses a default naming convention
            if not specified.   
            
        precision : int, optional, default: 6
            Number of decimals when writing the dataset file. 
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
        version = '# version: {} \n'.format(ci_version)
        
        # Parent
        parent = '# parent: {} \n'.format(self.trajectory.filename)
        
        # Feature scaling
        scaling = '# feature scaling: {} \n'
        if self.scaling is not None:
            scaling = scaling.format(self.scaling.full_name)
        else:
            scaling = scaling.format('none')
            
        # Dimensionality reduction
        dim_redux = '# dimensionality reduction: {} \n'
        if self.dim_redux is not None:
            dim_redux = dim_redux.format(self.dim_redux.full_name)
        else:
            dim_redux = dim_redux.format('none')
            
        # Clustering method
        clustering = '# clustering method: {} \n'.format(self.clustering.full_name)

        # Number of communities/clusters
        n_clusters = '# requested number of clusters: {} \n'.format(self.clustering.n_clusters)
        
        # Assemble header
        header = ''.join([date, version, parent, scaling, dim_redux, clustering, n_clusters])
        return header
    
    def _output_file(self, kind):
        """Returns path of output file."""
        return _output_path.format(filename=self.trajectory.filename,
                                   mode=self.descriptor.symbol,
                                   method=self.clustering.symbol,
                                   kind=kind)
    
    def __str__(self):
        rep = 'Optimization(filename="{}", descriptor="{}", scaling="{}", dim_redux="{}", clustering="{}", has_run={})'
        rep = rep.format(self.trajectory.filename,
                         self.descriptor.symbol,
                         self.scaling.symbol if self.scaling is not None else None,
                         self.dim_redux.symbol if self.dim_redux is not None else None,
                         self.clustering.symbol,
                         self._has_run)
        return rep
    
    def __repr__(self):
        return self.__str__()
