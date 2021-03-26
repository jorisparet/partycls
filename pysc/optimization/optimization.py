from pysc.trajectory import Trajectory
from pysc.descriptor import AngularDescriptor, RadialDescriptor, BondOrientationalDescriptor, LechnerDellagoDescriptor
from pysc.processing import KMeans, GaussianMixture, CommunityInference
from pysc.processing import PCA, TSNE, AutoEncoder
from pysc.processing import ZScore, MinMax
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
                                               'fields':[],
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
        Run the optimization.
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
        Later.
        """
        for key, val in kwargs.items():
            self.output_metadata[what][key] = val
        
    def disable_output(self):
        """
        Disable all outputs.
        """
        for key in self.output_metadata.keys():
            self.output_metadata[key]['enable'] = False

    def write_trajectory(self, filename=None, fmt='xyz', fields=[], precision=6, **kwargs):
        """
        Write trajectory to file.
        """
        if filename is None:
            filename = self._output_file(fmt)
        self.trajectory._write(filename, fmt=fmt, fields=fields, precision=precision)

    # TODO: more info needed in the log?
    def write_log(self, filename=None, precision=6, **kwargs):
        """
        Write a log file with all relevant information about the optimization.
        The log file can be written only if the optimization has been run.
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
        Later.
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
        if filename is None:
            filename = self._output_file('labels')
        with open(filename, 'w') as file:
            file.write("# title: clusters' labels\n")
            file.write(self._get_header())
            for ki in self.labels:
                file.write('{} \n'.format(ki))
                
    def write_dataset(self, filename=None, precision=6, **kwargs):
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
