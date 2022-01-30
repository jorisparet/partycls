URLS=[
"partycls/index.html",
"partycls/helpers.html",
"partycls/particle.html",
"partycls/core/index.html",
"partycls/core/utils.html",
"partycls/system.html",
"partycls/clustering.html",
"partycls/descriptor/index.html",
"partycls/descriptor/bo.html",
"partycls/descriptor/descriptor.html",
"partycls/descriptor/gr.html",
"partycls/descriptor/ba.html",
"partycls/descriptor/realspace_wrap.html",
"partycls/descriptor/dscribe.html",
"partycls/trajectory.html",
"partycls/feature_scaling.html",
"partycls/workflow.html",
"partycls/cell.html",
"partycls/dim_reduction.html"
];
INDEX=[
{
"ref":"partycls",
"url":0,
"doc":"partycls is a Python package for cluster analysis of systems of interacting particles. By grouping particles that share similar structural or dynamical properties, partycls enables rapid and unsupervised exploration of the system's relevant features. It provides descriptors suitable for applications in condensed matter physics, such as structural analysis of disordered or partially ordered materials, and integrates the necessary tools of unsupervised learning into a streamlined workflow."
},
{
"ref":"partycls.helpers",
"url":1,
"doc":"Various helper functions for visualization, cluster analysis, etc."
},
{
"ref":"partycls.helpers.AMI",
"url":1,
"doc":"Adjusted Mutual Information between two clusterings. Adjusted Mutual Information (AMI) is an adjustment of the Mutual Information (MI) score to account for chance. It accounts for the fact that the MI is generally higher for two clusterings with a larger number of clusters, regardless of whether there is actually more information shared. For two clusterings :math: U and :math: V , the AMI is given as AMI(U, V) = [MI(U, V) - E(MI(U, V ] / [avg(H(U), H(V - E(MI(U, V ] This metric is independent of the absolute values of the labels: a permutation of the class or cluster label values won't change the score value in any way. This metric is furthermore symmetric: switching  label_true with  label_pred will return the same score value. This can be useful to measure the agreement of two independent label assignments strategies on the same dataset when the real ground truth is not known. Be mindful that this function is an order of magnitude slower than other metrics, such as the Adjusted Rand Index. Read more in the :ref: User Guide   . Parameters      labels_true : int array, shape = [n_samples] A clustering of the data into disjoint subsets. labels_pred : int array-like of shape (n_samples,) A clustering of the data into disjoint subsets. average_method : str, default='arithmetic' How to compute the normalizer in the denominator. Possible options are 'min', 'geometric', 'arithmetic', and 'max'.  versionadded 0.20  versionchanged 0.22 The default value of  average_method changed from 'max' to 'arithmetic'. Returns    - ami: float (upperlimited by 1.0) The AMI returns a value of 1 when the two partitions are identical (ie perfectly matched). Random partitions (independent labellings) have an expected AMI around 0 on average hence can be negative. See Also     adjusted_rand_score : Adjusted Rand Index. mutual_info_score : Mutual Information (not adjusted for chance). Examples     Perfect labelings are both homogeneous and complete, hence have score 1.0 >>> from sklearn.metrics.cluster import adjusted_mutual_info_score >>> adjusted_mutual_info_score([0, 0, 1, 1], [0, 0, 1, 1])  .  doctest: +SKIP 1.0 >>> adjusted_mutual_info_score([0, 0, 1, 1], [1, 1, 0, 0])  .  doctest: +SKIP 1.0 If classes members are completely split across different clusters, the assignment is totally in-complete, hence the AMI is null >>> adjusted_mutual_info_score([0, 0, 0, 0], [0, 1, 2, 3])  .  doctest: +SKIP 0.0 References       [1]  Vinh, Epps, and Bailey, (2010). Information Theoretic Measures for Clusterings Comparison: Variants, Properties, Normalization and Correction for Chance, JMLR   _  [2]  Wikipedia entry for the Adjusted Mutual Information   _",
"func":1
},
{
"ref":"partycls.helpers.ARI",
"url":1,
"doc":"Rand index adjusted for chance. The Rand Index computes a similarity measure between two clusterings by considering all pairs of samples and counting pairs that are assigned in the same or different clusters in the predicted and true clusterings. The raw RI score is then \"adjusted for chance\" into the ARI score using the following scheme ARI = (RI - Expected_RI) / (max(RI) - Expected_RI) The adjusted Rand index is thus ensured to have a value close to 0.0 for random labeling independently of the number of clusters and samples and exactly 1.0 when the clusterings are identical (up to a permutation). ARI is a symmetric measure adjusted_rand_score(a, b)  adjusted_rand_score(b, a) Read more in the :ref: User Guide   . Parameters      labels_true : int array, shape = [n_samples] Ground truth class labels to be used as a reference labels_pred : array-like of shape (n_samples,) Cluster labels to evaluate Returns    - ARI : float Similarity score between -1.0 and 1.0. Random labelings have an ARI close to 0.0. 1.0 stands for perfect match. Examples     Perfectly matching labelings have a score of 1 even >>> from sklearn.metrics.cluster import adjusted_rand_score >>> adjusted_rand_score([0, 0, 1, 1], [0, 0, 1, 1]) 1.0 >>> adjusted_rand_score([0, 0, 1, 1], [1, 1, 0, 0]) 1.0 Labelings that assign all classes members to the same clusters are complete but may not always be pure, hence penalized >>> adjusted_rand_score([0, 0, 1, 2], [0, 0, 1, 1]) 0.57 . ARI is symmetric, so labelings that have pure clusters with members coming from the same classes but unnecessary splits are penalized >>> adjusted_rand_score([0, 0, 1, 1], [0, 0, 1, 2]) 0.57 . If classes members are completely split across different clusters, the assignment is totally incomplete, hence the ARI is very low >>> adjusted_rand_score([0, 0, 0, 0], [0, 1, 2, 3]) 0.0 References       [Hubert1985] L. Hubert and P. Arabie, Comparing Partitions, Journal of Classification 1985 https: link.springer.com/article/10.1007%2FBF01908075  [Steinley2004] D. Steinley, Properties of the Hubert-Arabie adjusted Rand index, Psychological Methods 2004  [wk] https: en.wikipedia.org/wiki/Rand_index Adjusted_Rand_index See Also     adjusted_mutual_info_score : Adjusted Mutual Information.",
"func":1
},
{
"ref":"partycls.helpers.silhouette_samples",
"url":1,
"doc":"Compute the Silhouette Coefficient for each sample. The Silhouette Coefficient is a measure of how well samples are clustered with samples that are similar to themselves. Clustering models with a high Silhouette Coefficient are said to be dense, where samples in the same cluster are similar to each other, and well separated, where samples in different clusters are not very similar to each other. The Silhouette Coefficient is calculated using the mean intra-cluster distance ( a ) and the mean nearest-cluster distance ( b ) for each sample. The Silhouette Coefficient for a sample is  (b - a) / max(a, b) . Note that Silhouette Coefficient is only defined if number of labels is 2    . Parameters      X : array-like of shape (n_samples_a, n_samples_a) if metric  \"precomputed\" or (n_samples_a, n_features) otherwise An array of pairwise distances between samples, or a feature array. labels : array-like of shape (n_samples,) Label values for each sample. metric : str or callable, default='euclidean' The metric to use when calculating distance between instances in a feature array. If metric is a string, it must be one of the options allowed by :func: sklearn.metrics.pairwise.pairwise_distances . If  X is the distance array itself, use \"precomputed\" as the metric. Precomputed distance matrices must have 0 along the diagonal.  kwds : optional keyword parameters Any further parameters are passed directly to the distance function. If using a  scipy.spatial.distance metric, the parameters are still metric dependent. See the scipy docs for usage examples. Returns    - silhouette : array-like of shape (n_samples,) Silhouette Coefficients for each sample. References       [1]  Peter J. Rousseeuw (1987). \"Silhouettes: a Graphical Aid to the Interpretation and Validation of Cluster Analysis\". Computational and Applied Mathematics 20: 53-65.   _  [2]  Wikipedia entry on the Silhouette Coefficient   _",
"func":1
},
{
"ref":"partycls.helpers.silhouette_score",
"url":1,
"doc":"Compute the mean Silhouette Coefficient of all samples. The Silhouette Coefficient is calculated using the mean intra-cluster distance ( a ) and the mean nearest-cluster distance ( b ) for each sample. The Silhouette Coefficient for a sample is  (b - a) / max(a, b) . To clarify,  b is the distance between a sample and the nearest cluster that the sample is not a part of. Note that Silhouette Coefficient is only defined if number of labels is  2   . Parameters      X : array-like of shape (n_samples_a, n_samples_a) if metric  \"precomputed\" or (n_samples_a, n_features) otherwise An array of pairwise distances between samples, or a feature array. labels : array-like of shape (n_samples,) Predicted labels for each sample. metric : str or callable, default='euclidean' The metric to use when calculating distance between instances in a feature array. If metric is a string, it must be one of the options allowed by :func: metrics.pairwise.pairwise_distances   . If  X is the distance array itself, use  metric=\"precomputed\" . sample_size : int, default=None The size of the sample to use when computing the Silhouette Coefficient on a random subset of the data. If  sample_size is None , no sampling is used. random_state : int, RandomState instance or None, default=None Determines random number generation for selecting a subset of samples. Used when  sample_size is not None . Pass an int for reproducible results across multiple function calls. See :term: Glossary   .  kwds : optional keyword parameters Any further parameters are passed directly to the distance function. If using a scipy.spatial.distance metric, the parameters are still metric dependent. See the scipy docs for usage examples. Returns    - silhouette : float Mean Silhouette Coefficient for all samples. References       [1]  Peter J. Rousseeuw (1987). \"Silhouettes: a Graphical Aid to the Interpretation and Validation of Cluster Analysis\". Computational and Applied Mathematics 20: 53-65.   _  [2]  Wikipedia entry on the Silhouette Coefficient   _",
"func":1
},
{
"ref":"partycls.helpers.show_matplotlib",
"url":1,
"doc":"Make a snapshot of the  system using matplotlib. The figure is returned for further customization or visualization in jupyter notebooks. Parameters      system : System An instance of  System . color : str Particle property to use for color coding, e.g. 'species', 'label'. view : str, optional View type, i.e. face of the box to show. Only works for a 3D system. The default is 'top'. palette : list, optional List of colors when coloring particles according to a discrete property, such as 'species' or 'label'. A default palette will be used if not specified. The default is None. cmap : str, optional Name of a matplotlib colormap to use when coloring particles according to a continuous property such as 'velocity' or 'energy'. List of available colormap can be found in  matplotlib.cm.cmaps_listed . The default is 'viridis'. outfile : str, optional Output filename to save the snapshot. The default is None (not saved). linewidth : int or float, optional The default is 0.5. alpha : int or float, optional Transparency parameter. The default is 1.0. show : bool, optional Show the snapshot when calling the function. The default is False. Returns    - fig : matplotlib.figure.Figure Figure of the snapshot.",
"func":1
},
{
"ref":"partycls.helpers.show_ovito",
"url":1,
"doc":"Make a snapshot of the  system using Ovito. The image is returned for further customization or visualization in jupyter notebooks. Parameters      system : System An instance of  System . color : str Particle property to use for color coding, e.g. 'species', 'label'. view : str, optional View type, i.e. face of the box to show. Only works for a 3D system. The default is 'top'. palette : list, optional List of colors when coloring particles according to a discrete property, such as 'species' or 'label'. Colors must be expressed in RGB format through tuples (e.g. palette=[(0,0,1), (1,0,0)]). A default palette will be used if not specified. The default is None. cmap : str, optional Name of a matplotlib colormap to use when coloring particles according to a continuous property such as 'velocity' or 'energy'. List of available colormap can be found in  matplotlib.cm.cmaps_listed . The default is 'viridis'. outfile : str, optional Output filename to save the snapshot. The default is None (not saved). size : tuple, optional Size of the image to render. The default is (640, 480). zoom : bool, optional Zoom on the simulation box. The default is True. Returns    - Image Rendered image.",
"func":1
},
{
"ref":"partycls.helpers.show_3dmol",
"url":1,
"doc":"Visualize the  system using 3dmol http: 3dmol.csb.pitt.edu/ The py3Dmol view is returned for further customization or visualization in jupyter notebooks. Parameters      system : System An instance of  System . color : str Particle property to use for color coding, e.g. 'species', 'label'. This property must be a string or an integer. palette : list, optional List of colors when coloring particles according to a discrete property, such as 'species' or 'label'. A default palette will be used if not specified. The default is None. Raises    ValueError If the  color parameter refers to a float particle property. Returns    - view : py3Dmol.view py3Dmol view.",
"func":1
},
{
"ref":"partycls.helpers.shannon_entropy",
"url":1,
"doc":"Shannon entropy of distribution p(x). Parameters      px : list or numpy.array Distribution p(x). dx : float, optional Differential of x. The default is 1.0. Returns    - S : float Shannon entropy.",
"func":1
},
{
"ref":"partycls.helpers.merge_clusters",
"url":1,
"doc":"Merge clusters into  n_clusters_min new clusters based on the probabilities that particles initially belong to each of the original clusters with a certain probability and using an entropy criterion. See https: doi.org/10.1198/jcgs.2010.08111 (Baudry et al.) Parameters      weights : list or numpy.ndarray Probabilities that each particle belongs to each cluster. If there are N particles, then the length of the list (or first dimension of the array) must be N. If there are K original clusters, each element of  weights (or the first dimension of the array) must be K.  weights[i][j] (list) or  weights[i,k] (array) is the probability that particle  i belongs to cluster  k before merging. For each particle, sum(weights[i]) = 1. n_clusters_min : int, optional Final number of clusters after merging. The default is 2. epsilon_ : float Small number (close to zero). This is needed as a replacement for zero when computing a logarithm to avoid errors. The default is 1e-15. Returns    - new_weights : numpy.ndarray New weights after merging. Same shape and interpretation as the  weights input parameter. new_labels : list New discrete labels based on the weights after merging.",
"func":1
},
{
"ref":"partycls.helpers.sort_clusters",
"url":1,
"doc":"Make a consistent labeling of the clusters based on their centroids by computing an associated numerical value as sorting criterion. By default, the labeling is based on the Shannon entropy of each cluster. Parameters      labels : list Original labels. centroids : numpy.ndarray Cluster centroids. func : function, optional Function used to associate a numerical value to each cluster, to be used as sorting criterion. This function must accept a list or a one dimensional array as parameter (this parameter being the coordinates of a given centroid). The default is shannon_entropy. Returns    - new_labels : list New labels based on centroid entropies. new_centroids : numpy.ndarray Centroids arranged in order of descending entropies.",
"func":1
},
{
"ref":"partycls.particle",
"url":2,
"doc":"Point particles in a cartesian reference frame. This class is inspired by the  atooms framework authored by Daniele Coslovich See https: framagit.org/atooms/atooms"
},
{
"ref":"partycls.particle.Particle",
"url":2,
"doc":"A particle is defined by its position, its type, and an optional cluster label (default is -1). Parameters      position : list of float or float array, optional, default: None The position of the particle. If not given, it will be set to [0.0, 0.0, 0.0]. species : str, optional, default: \"A\" Particle type / species. label : int, optional, default: -1 Cluster label of the particle. Default is -1 (i.e. not belonging to any cluster). radius : float, optional, defaut: 0.5 Particle radius. Attributes      position : float array The position of the particle. species : str Particle type / species. label : int Cluster label of the particle. radius : float Particle radius. index : int A unique index to identify the particle. Examples     >>> p = Particle([0.0, 0.0, 0.0], species='A') >>> p = Particle([0.0, 0.0], species='B')"
},
{
"ref":"partycls.particle.Particle.fold",
"url":2,
"doc":"Fold the particle position into the central cell. Parameters      cell : Cell Simulation cell. Returns    - None",
"func":1
},
{
"ref":"partycls.core",
"url":3,
"doc":""
},
{
"ref":"partycls.core.utils",
"url":4,
"doc":""
},
{
"ref":"partycls.core.utils.tipify",
"url":4,
"doc":"Convert a string  s into the best matching type. Parameters      s : str String to convert Returns    - int, float, or str Best-matching type for the input string  s .",
"func":1
},
{
"ref":"partycls.core.utils.standardize_condition",
"url":4,
"doc":"Check that the condition is correctly formated (i.e  _operator_  ). Parameters      condition : str condition. Raises    ValueError If condition is not valid or if the  is not recognized). Returns    - condition : str A standardized condition.",
"func":1
},
{
"ref":"partycls.system",
"url":5,
"doc":"The physical system at hand. The system of interest in a classical atomistic simulations is composed of interacting point particles, usually enclosed in a simulation cell. This class is inspired by the  atooms framework authored by Daniele Coslovich See https: framagit.org/atooms/atooms"
},
{
"ref":"partycls.system.System",
"url":5,
"doc":"A system is composed of a collection of particles that lie within an orthorhombic cell. Parameters      particle : list of  Particle , optional, default: None A list of instances of  Particle . cell : Cell, optional, default: None The cell (simulation box). Attributes      particle : list of  Particle All the particles in the system. cell : Cell The cell where all the particles lie. Examples     >>> p = [Particle(position=[0.0, 0.0, 0.0], species='A'), Particle(position=[1.0, 1.0, 1.0], species='B')] >>> c = Cell([5.0, 5.0, 5.0]) >>> sys = System(particle=p, cell=c)"
},
{
"ref":"partycls.system.System.n_dimensions",
"url":5,
"doc":"Number of spatial dimensions, guessed from the length of  particle[0].position ."
},
{
"ref":"partycls.system.System.density",
"url":5,
"doc":"Number density of the system. It will raise a ValueException if  cell is None."
},
{
"ref":"partycls.system.System.distinct_species",
"url":5,
"doc":"Sorted numpy array of all the distinct species in the system."
},
{
"ref":"partycls.system.System.pairs_of_species",
"url":5,
"doc":"List of all the possible pairs of species."
},
{
"ref":"partycls.system.System.pairs_of_species_id",
"url":5,
"doc":"List of all the possible pairs of species ID."
},
{
"ref":"partycls.system.System.chemical_fractions",
"url":5,
"doc":"Numpy array with the chemical fractions of each species in the system."
},
{
"ref":"partycls.system.System.get_property",
"url":5,
"doc":"Return a numpy array with the system property specified by  what . If  what is a particle property, return the property for all particles in the system, or for a given subset of particles specified by  subset . Parameters      what : str Requested system property.  what must be of the form \"particle. \" or \"cell. \". The following particle aliases are accepted: - 'position': 'particle.position' - 'pos': 'particle.position' - 'position[0]': 'particle.position[0]', - 'pos[0]': 'particle.position[0]' - 'x': 'particle.position[0]' - 'position[1]': 'particle.position[1]', - 'pos[1]': 'particle.position[1]' - 'y': 'particle.position[1]' - 'position[2]': 'particle.position[2]' - 'pos[2]': 'particle.position[2]' - 'z': 'particle.position[2]' - 'species': 'particle.species' - 'spe': 'particle.species' - 'label': 'particle.label' - 'index': 'particle.index' - 'mass': 'particle.mass' - 'radius': 'particle.radius' subset : str, optional Subset of particles for which the property must be dumped. Must be of the form \"particle. \" unless \" \" is an alias. The default is None (all particles will be included). This is ignored if  what is cell property. Returns    - to_dump : numpy.ndarray Array of the requested system property. Examples     >>> traj = Trajectory('trajectory.xyz') >>> sys = traj[0] >>> pos_0 = sys.get_property('position') >>> spe_0 = sys.get_property('species') >>> sides = sys.get_property('cell.side')",
"func":1
},
{
"ref":"partycls.system.System.dump",
"url":5,
"doc":"Alias for the method get_property.",
"func":1
},
{
"ref":"partycls.system.System.set_property",
"url":5,
"doc":"Set a system property  what to  value . If  what is a particle property, set the property for all the particles in the system or for a given subset of particles specified by  subset . Parameters      what : str Name of the property to set. This is considered to be a particle property by default, unless it starts with \"cell\", e.g. \"cell.side\". value : int, float, list, or numpy.ndarray Value(s) of the property to set. An instance of  int or  float will set the same value for all concerned particles. An instance of  list or  numpy.ndarray will assign a specific value to each particle. In this case, the size of  value should respect the number of concerned particles. subset : str, optional Particles to which the property must be set. The default is None. This is ignored if  what is cell property. Returns    - None. Examples     >>> sys.set_property('mass', 1.0) >>> sys.set_property('radius', 0.5, \"species  'A'\") >>> labels = [0, 1, 0]  3 particles in the subset >>> sys.set_property('label', labels, \"species  'B'\") >>> sys.set_property('cell.side[0]', 2.0)",
"func":1
},
{
"ref":"partycls.system.System.show",
"url":5,
"doc":"Show a snapshot of the system and color particles according to an arbitrary property, such as species, cluster label, etc. Current visualization backends are 'matplotlib', 'ovito' and '3dmol'. Parameters      backend : str, optional Name of the backend to use for visualization. The default is 'matplotlib'. color : str, optional Name of the particle property to use as basis for coloring the particles. This property must be defined for all the particles in the system. The default is 'species'.  kwargs : additional keyworded arguments (backend-dependent). Raises    ValueError In case of unknown  backend . Returns    - Figure or View (backend-dependent) Examples     >>> sys.show(frame=0, color='label', backend='3dmol') >>> sys.show(frame=1, color='energy', backend='matplotlib', cmap='viridis')",
"func":1
},
{
"ref":"partycls.system.System.fold",
"url":5,
"doc":"Fold the particle positions into the central cell. Returns    - None.",
"func":1
},
{
"ref":"partycls.clustering",
"url":6,
"doc":"Clustering algorithms."
},
{
"ref":"partycls.clustering.Clustering",
"url":6,
"doc":"Base class for clustering methods. If a scikit-learn compatible backend is available (see  backend parameter below), it will be used within Strategy. Parameters      n_clusters : int, optional Requested number of clusters. The default is 2. n_init : int, optional Number of times the clustering will be run with different seeds. The default is 1. backend : scikit-learn compatible backend, optional Backend used for the clustering method. If provided, it must be an object implementing an sklearn compatible interface, with a  fit() method and a  labels_ attribute. Duck typing is assumed. The default is None. Attributes      n_clusters : int, optional Number of clusters. n_init : int, optional Number of times the clustering is run. labels : list of int Cluster labels. The default is None. Initialized after the  fit method is called."
},
{
"ref":"partycls.clustering.Clustering.fit",
"url":6,
"doc":"Run a scikit-learn compatible clustering backend (if available) on  X . Subclasses implementing a specific clustering algorithm must override this method.",
"func":1
},
{
"ref":"partycls.clustering.Clustering.fractions",
"url":6,
"doc":""
},
{
"ref":"partycls.clustering.Clustering.populations",
"url":6,
"doc":""
},
{
"ref":"partycls.clustering.Clustering.centroids",
"url":6,
"doc":"Central feature vector of each cluster. Each object in the dataset over which the clustering was performed is assigned a discrete label. This label represents the index of the nearest cluster center to which this object belongs. The centroid (i.e. the cluster center), is thus the average feature vector of all the objects in the cluster. Cluster memberships of the objects are stored in the  labels attribute. Coordinates of the centroids can then be calculated for an arbitrary dataset  X , provided it has the same shape as the original dataset used for the clustering. Parameters      X : numpy.ndarray Array of features (dataset) for which to compute the centroids. Returns    - C_k : numpy.ndarray Cluster centroids. C_k[n] is the coordinates of the n-th cluster center.",
"func":1
},
{
"ref":"partycls.clustering.KMeans",
"url":6,
"doc":"KMeans clustering. This class relies on the class  KMeans from the machine learning package \"scikit-learn\". An instance of sklearn.cluster.KMeans is created when calling the  fit method, and is then accessible through the  backend attribute for later use. See scikit's documentation for more information on the original class."
},
{
"ref":"partycls.clustering.KMeans.fit",
"url":6,
"doc":"Run the K-Means algorithm on  X . The predicted labels are updated in the attribute  labels of the current instance of  KMeans .",
"func":1
},
{
"ref":"partycls.clustering.KMeans.centroids",
"url":6,
"doc":"Central feature vector of each cluster. Each object in the dataset over which the clustering was performed is assigned a discrete label. This label represents the index of the nearest cluster center to which this object belongs. The centroid (i.e. the cluster center), is thus the average feature vector of all the objects in the cluster. Cluster memberships of the objects are stored in the  labels attribute. Coordinates of the centroids can then be calculated for an arbitrary dataset  X , provided it has the same shape as the original dataset used for the clustering. Parameters      X : numpy.ndarray Array of features (dataset) for which to compute the centroids. Returns    - C_k : numpy.ndarray Cluster centroids. C_k[n] is the coordinates of the n-th cluster center.",
"func":1
},
{
"ref":"partycls.clustering.GaussianMixture",
"url":6,
"doc":"Gaussian Mixture. This class relies on the class  GaussianMixture from the machine learning package \"scikit-learn\". An instance of sklearn.mixture.GaussianMixture is created when calling the  fit method, and is then accessible through the  backend attribute for later use. See scikit's documentation for more information on the original class."
},
{
"ref":"partycls.clustering.GaussianMixture.fit",
"url":6,
"doc":"Run the EM algorithm on  X using a mixture of Gaussians. The predicted labels are updated in the attribute  labels of the current instance of  GaussianMixture .",
"func":1
},
{
"ref":"partycls.clustering.GaussianMixture.centroids",
"url":6,
"doc":"Central feature vector of each cluster. Each object in the dataset over which the clustering was performed is assigned a discrete label. This label represents the index of the nearest cluster center to which this object belongs. The centroid (i.e. the cluster center), is thus the average feature vector of all the objects in the cluster. Cluster memberships of the objects are stored in the  labels attribute. Coordinates of the centroids can then be calculated for an arbitrary dataset  X , provided it has the same shape as the original dataset used for the clustering. Parameters      X : numpy.ndarray Array of features (dataset) for which to compute the centroids. Returns    - C_k : numpy.ndarray Cluster centroids. C_k[n] is the coordinates of the n-th cluster center.",
"func":1
},
{
"ref":"partycls.clustering.CommunityInference",
"url":6,
"doc":"Community Inference is a hard clustering method based on information theory. See \"Paret et al. https: doi.org/10.1063/5.0004732\" for more details."
},
{
"ref":"partycls.clustering.CommunityInference.fit",
"url":6,
"doc":"Community inference algorithm.",
"func":1
},
{
"ref":"partycls.clustering.CommunityInference.centroids",
"url":6,
"doc":"Central feature vector of each cluster. Each object in the dataset over which the clustering was performed is assigned a discrete label. This label represents the index of the nearest cluster center to which this object belongs. The centroid (i.e. the cluster center), is thus the average feature vector of all the objects in the cluster. Cluster memberships of the objects are stored in the  labels attribute. Coordinates of the centroids can then be calculated for an arbitrary dataset  X , provided it has the same shape as the original dataset used for the clustering. Parameters      X : numpy.ndarray Array of features (dataset) for which to compute the centroids. Returns    - C_k : numpy.ndarray Cluster centroids. C_k[n] is the coordinates of the n-th cluster center.",
"func":1
},
{
"ref":"partycls.descriptor",
"url":7,
"doc":"Structural descriptors."
},
{
"ref":"partycls.descriptor.bo",
"url":8,
"doc":""
},
{
"ref":"partycls.descriptor.bo.BondOrientationalDescriptor",
"url":8,
"doc":"Structural descriptor based on bond order parameters as defined by Steinhardt et al. (https: doi.org/10.1103%2FPhysRevB.28.784). See the parent class for more details. Parameters      trajectory : str or an instance of  Trajectory . Trajectory on which the structural descriptor will be computed. lmin : int, default: 0 Minimum degree. This set the lower bound of the grid. lmax : int, default: 8 Minimum degree. This set the upper bound of the grid. orders: list, default: None Specific values of orders to compute, e.g. orders=[4,6]. This has the priority over  lmin and  lmax . Attributes      trajectory : Trajectory Trajectory on which the structural descriptor will be computed. active_filters : list of str All the active filters on both groups prior to the computation of the descriptor. dimension : int Spatial dimension of the descriptor (2 or 3). grid : array Grid over which the structural features will be computed. features : ndarray Array of all the structural features for the particles in group=0 in accordance with the defined filters (if any). This attribute is initialized when the method  compute is called (default value is None). cutoffs : list of float List of cutoff distances to identify the nearest neighbors using the fixed-cutoff ('FC') method. nearest_neighbors_method : str, default: 'FC' Nearest neighbor method, 'FC' or 'SANN'. Examples:     - >>> D = BondOrientationalDescriptor('trajectory.xyz', orders=[4,6]) >>> D.nearest_neighbors_method = 'FC' >>> D.compute()"
},
{
"ref":"partycls.descriptor.bo.BondOrientationalDescriptor.name",
"url":8,
"doc":""
},
{
"ref":"partycls.descriptor.bo.BondOrientationalDescriptor.symbol",
"url":8,
"doc":""
},
{
"ref":"partycls.descriptor.bo.BondOrientationalDescriptor.orders",
"url":8,
"doc":""
},
{
"ref":"partycls.descriptor.bo.BondOrientationalDescriptor.compute",
"url":8,
"doc":"",
"func":1
},
{
"ref":"partycls.descriptor.bo.BondOrientationalDescriptor.set_cutoff",
"url":9,
"doc":"Set the nearest-neighbor cutoff for the pair of species (s1, s2). The cutoff of the mirror pair (s2, s1) is set automatically if the  mirror parameter is True (default). Parameters      s1 : str Symbol of the first species. s2 : str Symbol of the second species. rcut : float Value of the cutoff for the pair (s1,s2). mirror : bool, optional Set the cutoff for the mirror pair (s2,s1). The default is True. Returns    - None.",
"func":1
},
{
"ref":"partycls.descriptor.bo.BondOrientationalDescriptor.nearest_neighbors",
"url":9,
"doc":"Compute the nearest neighbors of particles in group=0 using one of the following methods: - \"Fixed cutoff\" (method='FC'): uses the partial radial distribution functions to compute the cutoffs between each possible pair of species (s1, s2) ; - \"Solid-Angle based Nearest Neighbors\" (method='SANN'): see van Meel et al. (https: doi.org/10.1063/1.4729313) ; Parameters      method : str, optional Method to identify nearest neighbors. Must be 'FC' or 'SANN'. The default is 'FC'. Returns    - None.",
"func":1
},
{
"ref":"partycls.descriptor.bo.BondOrientationalDescriptor.add_filter",
"url":9,
"doc":"Add a filter on the group (0 or 1) to select the subset of particles that respects the provided condition. Parameters      condition : str The condition should have the following format:  _operator_  where: -  is a particle property (accepts aliases) ; - _operator_ is a logical operator ( =, >) ; -  is the corresponding value of  with the proper type ; group : int, optional Index of the group to which the filter must be applied. The default is 0. Returns    - None. Examples:     - >>> S = StructuralDescriptor('trajectory.xyz') >>> S.add_filter(\"particle.radius  >> S.add_filter(\"species  'A'\", group=1) >>> S.add_filter(\"x < 0\", group=0)  particles on the left side of the box",
"func":1
},
{
"ref":"partycls.descriptor.bo.BondOrientationalDescriptor.clear_filters",
"url":9,
"doc":"Clear all active filters on  group . All particles are included again in  group . Parameters      group : int, optional Index of the group on which to clear the filters. The default is 0. Returns    - None.",
"func":1
},
{
"ref":"partycls.descriptor.bo.BondOrientationalDescriptor.clear_all_filters",
"url":9,
"doc":"Clear all active filters in both groups. All particles are included in both groups again. Returns    - None.",
"func":1
},
{
"ref":"partycls.descriptor.bo.BondOrientationalDescriptor.group_size",
"url":9,
"doc":"Return the number of particles in  group .",
"func":1
},
{
"ref":"partycls.descriptor.bo.BondOrientationalDescriptor.get_group_property",
"url":9,
"doc":"Return a list of numpy arrays with the properties of the particles in group  group . The list size is the number of systems in the trajectory. Parameters      what : str Requested particle property.  what must be a particle property or an alias. The following particle aliases are accepted: - 'position': 'particle.position' - 'pos': 'particle.position' - 'position[0]': 'particle.position[0]', - 'pos[0]': 'particle.position[0]' - 'x': 'particle.position[0]' - 'position[1]': 'particle.position[1]', - 'pos[1]': 'particle.position[1]' - 'y': 'particle.position[1]' - 'position[2]': 'particle.position[2]' - 'pos[2]': 'particle.position[2]' - 'z': 'particle.position[2]' - 'species': 'particle.species' - 'spe': 'particle.species' - 'label': 'particle.label' - 'index': 'particle.index' - 'mass': 'particle.mass' - 'radius': 'particle.radius' Returns    - to_dump : list List of the requested particle property with length equal to the number of frames in the trajectory. Each element of the list is a numpy.ndarray of the requested particle property. Examples     >>> traj = Trajectory('trajectory.xyz') >>> D = StructuralDescriptor(traj) >>> D.get_group_property('position', 0) >>> D.get_group_property('x', 1) >>> D.get_group_property('energy', 0)",
"func":1
},
{
"ref":"partycls.descriptor.bo.BondOrientationalDescriptor.dump",
"url":9,
"doc":"Alias for the method get_group_property.",
"func":1
},
{
"ref":"partycls.descriptor.bo.BondOrientationalDescriptor.group_fraction",
"url":9,
"doc":"Return the fraction of particles inside  group over the whole trajectory.",
"func":1
},
{
"ref":"partycls.descriptor.bo.BondOrientationalDescriptor.size",
"url":9,
"doc":"Total number of particles in the descriptor (i.e. in group=0)."
},
{
"ref":"partycls.descriptor.bo.BondOrientationalDescriptor.n_features",
"url":9,
"doc":"Number of features of the descriptor."
},
{
"ref":"partycls.descriptor.bo.BondOrientationalDescriptor.average",
"url":9,
"doc":"Average feature vector of the descriptor."
},
{
"ref":"partycls.descriptor.bo.BondOrientationalDescriptor.normalize",
"url":9,
"doc":"Generic normalization function for child classes. Returns the input distribution unchanged.",
"func":1
},
{
"ref":"partycls.descriptor.bo.LechnerDellagoDescriptor",
"url":8,
"doc":"Structural descriptor based on locally averaged bond order parameters as defined by Lechner & Dellago (https: doi.org/10.1063/1.2977970). See the parent class for more details. Parameters      trajectory : str or an instance of  Trajectory . Trajectory on which the structural descriptor will be computed. lmin : int, default: 0 Minimum degree. This set the lower bound of the grid. lmax : int, default: 8 Minimum degree. This set the upper bound of the grid. orders: list, default: None Specific values of orders to compute, e.g. orders=[4,6]. This has the priority over  lmin and  lmax . Attributes      trajectory : Trajectory Trajectory on which the structural descriptor will be computed. active_filters : list of str All the active filters on both groups prior to the computation of the descriptor. dimension : int Spatial dimension of the descriptor (2 or 3). grid : array Grid over which the structural features will be computed. features : ndarray Array of all the structural features for the particles in group=0 in accordance with the defined filters (if any). This attribute is initialized when the method  compute is called (default value is None). cutoffs : list of float List of cutoff distances to identify the nearest neighbors using the fixed-cutoff ('FC') method. nearest_neighbors_method : str, default: 'FC' Nearest neighbor method, 'FC' or 'SANN'. Examples:     - >>> D = LechnerDellagoDescriptor('trajectory.xyz', orders=[4,6]) >>> D.nearest_neighbors_method = 'FC' >>> D.add_filter(\"species  'A'\") >>> D.compute()"
},
{
"ref":"partycls.descriptor.bo.LechnerDellagoDescriptor.name",
"url":8,
"doc":""
},
{
"ref":"partycls.descriptor.bo.LechnerDellagoDescriptor.symbol",
"url":8,
"doc":""
},
{
"ref":"partycls.descriptor.bo.LechnerDellagoDescriptor.compute",
"url":8,
"doc":"",
"func":1
},
{
"ref":"partycls.descriptor.bo.LechnerDellagoDescriptor.set_cutoff",
"url":9,
"doc":"Set the nearest-neighbor cutoff for the pair of species (s1, s2). The cutoff of the mirror pair (s2, s1) is set automatically if the  mirror parameter is True (default). Parameters      s1 : str Symbol of the first species. s2 : str Symbol of the second species. rcut : float Value of the cutoff for the pair (s1,s2). mirror : bool, optional Set the cutoff for the mirror pair (s2,s1). The default is True. Returns    - None.",
"func":1
},
{
"ref":"partycls.descriptor.bo.LechnerDellagoDescriptor.nearest_neighbors",
"url":9,
"doc":"Compute the nearest neighbors of particles in group=0 using one of the following methods: - \"Fixed cutoff\" (method='FC'): uses the partial radial distribution functions to compute the cutoffs between each possible pair of species (s1, s2) ; - \"Solid-Angle based Nearest Neighbors\" (method='SANN'): see van Meel et al. (https: doi.org/10.1063/1.4729313) ; Parameters      method : str, optional Method to identify nearest neighbors. Must be 'FC' or 'SANN'. The default is 'FC'. Returns    - None.",
"func":1
},
{
"ref":"partycls.descriptor.bo.LechnerDellagoDescriptor.add_filter",
"url":9,
"doc":"Add a filter on the group (0 or 1) to select the subset of particles that respects the provided condition. Parameters      condition : str The condition should have the following format:  _operator_  where: -  is a particle property (accepts aliases) ; - _operator_ is a logical operator ( =, >) ; -  is the corresponding value of  with the proper type ; group : int, optional Index of the group to which the filter must be applied. The default is 0. Returns    - None. Examples:     - >>> S = StructuralDescriptor('trajectory.xyz') >>> S.add_filter(\"particle.radius  >> S.add_filter(\"species  'A'\", group=1) >>> S.add_filter(\"x < 0\", group=0)  particles on the left side of the box",
"func":1
},
{
"ref":"partycls.descriptor.bo.LechnerDellagoDescriptor.clear_filters",
"url":9,
"doc":"Clear all active filters on  group . All particles are included again in  group . Parameters      group : int, optional Index of the group on which to clear the filters. The default is 0. Returns    - None.",
"func":1
},
{
"ref":"partycls.descriptor.bo.LechnerDellagoDescriptor.clear_all_filters",
"url":9,
"doc":"Clear all active filters in both groups. All particles are included in both groups again. Returns    - None.",
"func":1
},
{
"ref":"partycls.descriptor.bo.LechnerDellagoDescriptor.group_size",
"url":9,
"doc":"Return the number of particles in  group .",
"func":1
},
{
"ref":"partycls.descriptor.bo.LechnerDellagoDescriptor.get_group_property",
"url":9,
"doc":"Return a list of numpy arrays with the properties of the particles in group  group . The list size is the number of systems in the trajectory. Parameters      what : str Requested particle property.  what must be a particle property or an alias. The following particle aliases are accepted: - 'position': 'particle.position' - 'pos': 'particle.position' - 'position[0]': 'particle.position[0]', - 'pos[0]': 'particle.position[0]' - 'x': 'particle.position[0]' - 'position[1]': 'particle.position[1]', - 'pos[1]': 'particle.position[1]' - 'y': 'particle.position[1]' - 'position[2]': 'particle.position[2]' - 'pos[2]': 'particle.position[2]' - 'z': 'particle.position[2]' - 'species': 'particle.species' - 'spe': 'particle.species' - 'label': 'particle.label' - 'index': 'particle.index' - 'mass': 'particle.mass' - 'radius': 'particle.radius' Returns    - to_dump : list List of the requested particle property with length equal to the number of frames in the trajectory. Each element of the list is a numpy.ndarray of the requested particle property. Examples     >>> traj = Trajectory('trajectory.xyz') >>> D = StructuralDescriptor(traj) >>> D.get_group_property('position', 0) >>> D.get_group_property('x', 1) >>> D.get_group_property('energy', 0)",
"func":1
},
{
"ref":"partycls.descriptor.bo.LechnerDellagoDescriptor.dump",
"url":9,
"doc":"Alias for the method get_group_property.",
"func":1
},
{
"ref":"partycls.descriptor.bo.LechnerDellagoDescriptor.group_fraction",
"url":9,
"doc":"Return the fraction of particles inside  group over the whole trajectory.",
"func":1
},
{
"ref":"partycls.descriptor.bo.LechnerDellagoDescriptor.size",
"url":9,
"doc":"Total number of particles in the descriptor (i.e. in group=0)."
},
{
"ref":"partycls.descriptor.bo.LechnerDellagoDescriptor.n_features",
"url":9,
"doc":"Number of features of the descriptor."
},
{
"ref":"partycls.descriptor.bo.LechnerDellagoDescriptor.average",
"url":9,
"doc":"Average feature vector of the descriptor."
},
{
"ref":"partycls.descriptor.bo.LechnerDellagoDescriptor.normalize",
"url":9,
"doc":"Generic normalization function for child classes. Returns the input distribution unchanged.",
"func":1
},
{
"ref":"partycls.descriptor.gr",
"url":10,
"doc":""
},
{
"ref":"partycls.descriptor.gr.RadialDescriptor",
"url":10,
"doc":"Structural descriptor based on radial correlations between particles. See the parent class for more details. Parameters      trajectory : str or an instance of  Trajectory . Trajectory on which the structural descriptor will be computed. dr : float Bin width. n_shells : int, default: 3 Number of coordination shells (based on the RDF of group=0). This sets the upper bound for the distance up to which correlations are computed. bounds : tuple, default: None Lower and upper bounds to describe the radial correlations. If set, this has the priority over  n_shells . Attributes      trajectory : Trajectory Trajectory on which the structural descriptor will be computed. active_filters : list of str All the active filters on both groups prior to the computation of the descriptor. dimension : int Spatial dimension of the descriptor (2 or 3). grid : array Grid over which the structural features will be computed. features : ndarray Array of all the structural features for the particles in group=0 in accordance with the defined filters (if any). This attribute is initialized when the method  compute is called (default value is None). Examples:     - >>> D = RadialDescriptor('trajectory.xyz', bounds=(0.0,3.0 >>> D.add_filter(\"species  'A'\") >>> D.compute()"
},
{
"ref":"partycls.descriptor.gr.RadialDescriptor.name",
"url":10,
"doc":""
},
{
"ref":"partycls.descriptor.gr.RadialDescriptor.symbol",
"url":10,
"doc":""
},
{
"ref":"partycls.descriptor.gr.RadialDescriptor.n_shells",
"url":10,
"doc":"Upper bound for correlation expressed in number of coordinations shells."
},
{
"ref":"partycls.descriptor.gr.RadialDescriptor.bounds",
"url":10,
"doc":"Lower and upper bounds to describe the radial correlations."
},
{
"ref":"partycls.descriptor.gr.RadialDescriptor.dr",
"url":10,
"doc":"Grid spacing."
},
{
"ref":"partycls.descriptor.gr.RadialDescriptor.compute",
"url":10,
"doc":"Compute the radial correlations for the particles in group=0 in the range of distances given by  bounds . Returns    - features : numpy.ndarray Radial correlations.",
"func":1
},
{
"ref":"partycls.descriptor.gr.RadialDescriptor.normalize",
"url":10,
"doc":"Normalize a radial distribution. Parameters      distribution : array Distribution to normalize. method : str, optional Normalization method: - method='r2': returns r^2  g(r) (default); - method='gr' : returns the standard g(r) ; Raises    ValueError If  method is invalid. Returns    - array Normalized distribution.",
"func":1
},
{
"ref":"partycls.descriptor.gr.RadialDescriptor.add_filter",
"url":9,
"doc":"Add a filter on the group (0 or 1) to select the subset of particles that respects the provided condition. Parameters      condition : str The condition should have the following format:  _operator_  where: -  is a particle property (accepts aliases) ; - _operator_ is a logical operator ( =, >) ; -  is the corresponding value of  with the proper type ; group : int, optional Index of the group to which the filter must be applied. The default is 0. Returns    - None. Examples:     - >>> S = StructuralDescriptor('trajectory.xyz') >>> S.add_filter(\"particle.radius  >> S.add_filter(\"species  'A'\", group=1) >>> S.add_filter(\"x < 0\", group=0)  particles on the left side of the box",
"func":1
},
{
"ref":"partycls.descriptor.gr.RadialDescriptor.clear_filters",
"url":9,
"doc":"Clear all active filters on  group . All particles are included again in  group . Parameters      group : int, optional Index of the group on which to clear the filters. The default is 0. Returns    - None.",
"func":1
},
{
"ref":"partycls.descriptor.gr.RadialDescriptor.clear_all_filters",
"url":9,
"doc":"Clear all active filters in both groups. All particles are included in both groups again. Returns    - None.",
"func":1
},
{
"ref":"partycls.descriptor.gr.RadialDescriptor.group_size",
"url":9,
"doc":"Return the number of particles in  group .",
"func":1
},
{
"ref":"partycls.descriptor.gr.RadialDescriptor.get_group_property",
"url":9,
"doc":"Return a list of numpy arrays with the properties of the particles in group  group . The list size is the number of systems in the trajectory. Parameters      what : str Requested particle property.  what must be a particle property or an alias. The following particle aliases are accepted: - 'position': 'particle.position' - 'pos': 'particle.position' - 'position[0]': 'particle.position[0]', - 'pos[0]': 'particle.position[0]' - 'x': 'particle.position[0]' - 'position[1]': 'particle.position[1]', - 'pos[1]': 'particle.position[1]' - 'y': 'particle.position[1]' - 'position[2]': 'particle.position[2]' - 'pos[2]': 'particle.position[2]' - 'z': 'particle.position[2]' - 'species': 'particle.species' - 'spe': 'particle.species' - 'label': 'particle.label' - 'index': 'particle.index' - 'mass': 'particle.mass' - 'radius': 'particle.radius' Returns    - to_dump : list List of the requested particle property with length equal to the number of frames in the trajectory. Each element of the list is a numpy.ndarray of the requested particle property. Examples     >>> traj = Trajectory('trajectory.xyz') >>> D = StructuralDescriptor(traj) >>> D.get_group_property('position', 0) >>> D.get_group_property('x', 1) >>> D.get_group_property('energy', 0)",
"func":1
},
{
"ref":"partycls.descriptor.gr.RadialDescriptor.dump",
"url":9,
"doc":"Alias for the method get_group_property.",
"func":1
},
{
"ref":"partycls.descriptor.gr.RadialDescriptor.group_fraction",
"url":9,
"doc":"Return the fraction of particles inside  group over the whole trajectory.",
"func":1
},
{
"ref":"partycls.descriptor.gr.RadialDescriptor.size",
"url":9,
"doc":"Total number of particles in the descriptor (i.e. in group=0)."
},
{
"ref":"partycls.descriptor.gr.RadialDescriptor.n_features",
"url":9,
"doc":"Number of features of the descriptor."
},
{
"ref":"partycls.descriptor.gr.RadialDescriptor.average",
"url":9,
"doc":"Average feature vector of the descriptor."
},
{
"ref":"partycls.descriptor.descriptor",
"url":9,
"doc":""
},
{
"ref":"partycls.descriptor.descriptor.StructuralDescriptor",
"url":9,
"doc":"Base class for structural descriptors. The descriptor is calculated for the provided trajectory  trajectory . This can be: - an object implementing the  Trajectory interface ; - the path to a trajectory file in a format recognized by partyclsl ; A structural descriptor S(x) is a collection of N individual empirical correlation functions {s_i(x)} at the particle level, defined over a grid {x_j} of M features. These are stored in the  features array as a matrix usually refered to as the \"data set\": s_0(x_0) s_0(x_1)  . s_0(x_M) s_1(x_0) s_1(x_1)  . s_1(x_M)  .  .  . s_N(x_0) s_N(x_1)  . s_N(x_M) The  features array is None by default and is computed only when the  compute() method is called. The correlations can be calculated between two arbitrary subsets of particles called \"groups\": - group 0 is the main group, i.e. particles for which the correlations are being calculated ; - group 1 is the secondary group, i.e. particles that are being considered when calculating the correlations ; These groups are formed by adding filters on particles' properties (species, radius, position, etc.). Parameters      trajectory : str or an instance of  Trajectory . Trajectory on which the structural descriptor will be computed. Attributes      trajectory : Trajectory Trajectory on which the structural descriptor will be computed. active_filters : list of str All the active filters on both groups prior to the computation of the descriptor. dimension : int Spatial dimension of the descriptor (2 or 3). grid : array Grid over which the structural features will be computed. features : ndarray Array of all the structural features for the particles in group=0 in accordance with the defined filters (if any). This attribute is initialized when the method  compute is called (default value is None). groups : tuple Composition of the groups: groups[0] and groups[1] contain lists of all the  Particle objects in groups 0 and 1 respectively. Each element of the tuple is a list of  Particle in  trajectory , e.g. groups[0][0] is the list of all the particles in the first frame of  trajectory that belong to group=0. Examples:     - >>> D = StructuralDescriptor('trajectory.xyz') >>> D.add_filter(\"species  'A'\", group=0) >>> D.add_filter(\"species  'B'\", group=1) >>> D.active_filters [(\"particle.species  'A'\", 0), (\"particle.species  'B'\", 1)] >>> D.clear_filters(0) >>> D.active_filters [(\"particle.species  'B'\", 1)]"
},
{
"ref":"partycls.descriptor.descriptor.StructuralDescriptor.add_filter",
"url":9,
"doc":"Add a filter on the group (0 or 1) to select the subset of particles that respects the provided condition. Parameters      condition : str The condition should have the following format:  _operator_  where: -  is a particle property (accepts aliases) ; - _operator_ is a logical operator ( =, >) ; -  is the corresponding value of  with the proper type ; group : int, optional Index of the group to which the filter must be applied. The default is 0. Returns    - None. Examples:     - >>> S = StructuralDescriptor('trajectory.xyz') >>> S.add_filter(\"particle.radius  >> S.add_filter(\"species  'A'\", group=1) >>> S.add_filter(\"x < 0\", group=0)  particles on the left side of the box",
"func":1
},
{
"ref":"partycls.descriptor.descriptor.StructuralDescriptor.clear_filters",
"url":9,
"doc":"Clear all active filters on  group . All particles are included again in  group . Parameters      group : int, optional Index of the group on which to clear the filters. The default is 0. Returns    - None.",
"func":1
},
{
"ref":"partycls.descriptor.descriptor.StructuralDescriptor.clear_all_filters",
"url":9,
"doc":"Clear all active filters in both groups. All particles are included in both groups again. Returns    - None.",
"func":1
},
{
"ref":"partycls.descriptor.descriptor.StructuralDescriptor.group_size",
"url":9,
"doc":"Return the number of particles in  group .",
"func":1
},
{
"ref":"partycls.descriptor.descriptor.StructuralDescriptor.get_group_property",
"url":9,
"doc":"Return a list of numpy arrays with the properties of the particles in group  group . The list size is the number of systems in the trajectory. Parameters      what : str Requested particle property.  what must be a particle property or an alias. The following particle aliases are accepted: - 'position': 'particle.position' - 'pos': 'particle.position' - 'position[0]': 'particle.position[0]', - 'pos[0]': 'particle.position[0]' - 'x': 'particle.position[0]' - 'position[1]': 'particle.position[1]', - 'pos[1]': 'particle.position[1]' - 'y': 'particle.position[1]' - 'position[2]': 'particle.position[2]' - 'pos[2]': 'particle.position[2]' - 'z': 'particle.position[2]' - 'species': 'particle.species' - 'spe': 'particle.species' - 'label': 'particle.label' - 'index': 'particle.index' - 'mass': 'particle.mass' - 'radius': 'particle.radius' Returns    - to_dump : list List of the requested particle property with length equal to the number of frames in the trajectory. Each element of the list is a numpy.ndarray of the requested particle property. Examples     >>> traj = Trajectory('trajectory.xyz') >>> D = StructuralDescriptor(traj) >>> D.get_group_property('position', 0) >>> D.get_group_property('x', 1) >>> D.get_group_property('energy', 0)",
"func":1
},
{
"ref":"partycls.descriptor.descriptor.StructuralDescriptor.dump",
"url":9,
"doc":"Alias for the method get_group_property.",
"func":1
},
{
"ref":"partycls.descriptor.descriptor.StructuralDescriptor.group_fraction",
"url":9,
"doc":"Return the fraction of particles inside  group over the whole trajectory.",
"func":1
},
{
"ref":"partycls.descriptor.descriptor.StructuralDescriptor.size",
"url":9,
"doc":"Total number of particles in the descriptor (i.e. in group=0)."
},
{
"ref":"partycls.descriptor.descriptor.StructuralDescriptor.n_features",
"url":9,
"doc":"Number of features of the descriptor."
},
{
"ref":"partycls.descriptor.descriptor.StructuralDescriptor.average",
"url":9,
"doc":"Average feature vector of the descriptor."
},
{
"ref":"partycls.descriptor.descriptor.StructuralDescriptor.compute",
"url":9,
"doc":"",
"func":1
},
{
"ref":"partycls.descriptor.descriptor.StructuralDescriptor.normalize",
"url":9,
"doc":"Generic normalization function for child classes. Returns the input distribution unchanged.",
"func":1
},
{
"ref":"partycls.descriptor.descriptor.AngularStructuralDescriptor",
"url":9,
"doc":"Base class for angular structural descriptors. See the parent class for more details. Descriptors that exploit angular correlations and require nearest-neighbors information will inherit of this class. Two methods to identify nearest-neighbors are available: - \"Fixed cutoff\" (symbol: 'FC'): uses the partial radial distribution functions to compute the cutoffs between each possible pair of species (s1, s2) ; - \"Solid-Angle based Nearest Neighbors\" (symbol: 'SANN'): see van Meel et al. (https: doi.org/10.1063/1.4729313) The nearest-neighbors method can be changed by modifying the attribute  nearest_neighbors_method to 'FC' (default) or 'SANN'. When using the 'FC' method, it is also possible to specify the cutoffs manually for a pair of species (s1, s2) by using the method  set_cutoff . The cutoffs that were not set manually will be computed automatically. Parameters      trajectory : str or an instance of  Trajectory . Trajectory on which the structural descriptor will be computed. Attributes      cutoffs : list of float List of cutoff distances to identify the nearest neighbors using the fixed-cutoff ('FC') method. nearest_neighbors_method : str, default: 'FC' Nearest neighbor method, 'FC' or 'SANN'. neighbors : list Lists of nearest neighbors for all the particles in group=0. Empty by default and filled when calling the method  nearest_neighbors ."
},
{
"ref":"partycls.descriptor.descriptor.AngularStructuralDescriptor.set_cutoff",
"url":9,
"doc":"Set the nearest-neighbor cutoff for the pair of species (s1, s2). The cutoff of the mirror pair (s2, s1) is set automatically if the  mirror parameter is True (default). Parameters      s1 : str Symbol of the first species. s2 : str Symbol of the second species. rcut : float Value of the cutoff for the pair (s1,s2). mirror : bool, optional Set the cutoff for the mirror pair (s2,s1). The default is True. Returns    - None.",
"func":1
},
{
"ref":"partycls.descriptor.descriptor.AngularStructuralDescriptor.nearest_neighbors",
"url":9,
"doc":"Compute the nearest neighbors of particles in group=0 using one of the following methods: - \"Fixed cutoff\" (method='FC'): uses the partial radial distribution functions to compute the cutoffs between each possible pair of species (s1, s2) ; - \"Solid-Angle based Nearest Neighbors\" (method='SANN'): see van Meel et al. (https: doi.org/10.1063/1.4729313) ; Parameters      method : str, optional Method to identify nearest neighbors. Must be 'FC' or 'SANN'. The default is 'FC'. Returns    - None.",
"func":1
},
{
"ref":"partycls.descriptor.descriptor.AngularStructuralDescriptor.add_filter",
"url":9,
"doc":"Add a filter on the group (0 or 1) to select the subset of particles that respects the provided condition. Parameters      condition : str The condition should have the following format:  _operator_  where: -  is a particle property (accepts aliases) ; - _operator_ is a logical operator ( =, >) ; -  is the corresponding value of  with the proper type ; group : int, optional Index of the group to which the filter must be applied. The default is 0. Returns    - None. Examples:     - >>> S = StructuralDescriptor('trajectory.xyz') >>> S.add_filter(\"particle.radius  >> S.add_filter(\"species  'A'\", group=1) >>> S.add_filter(\"x < 0\", group=0)  particles on the left side of the box",
"func":1
},
{
"ref":"partycls.descriptor.descriptor.AngularStructuralDescriptor.clear_filters",
"url":9,
"doc":"Clear all active filters on  group . All particles are included again in  group . Parameters      group : int, optional Index of the group on which to clear the filters. The default is 0. Returns    - None.",
"func":1
},
{
"ref":"partycls.descriptor.descriptor.AngularStructuralDescriptor.clear_all_filters",
"url":9,
"doc":"Clear all active filters in both groups. All particles are included in both groups again. Returns    - None.",
"func":1
},
{
"ref":"partycls.descriptor.descriptor.AngularStructuralDescriptor.group_size",
"url":9,
"doc":"Return the number of particles in  group .",
"func":1
},
{
"ref":"partycls.descriptor.descriptor.AngularStructuralDescriptor.get_group_property",
"url":9,
"doc":"Return a list of numpy arrays with the properties of the particles in group  group . The list size is the number of systems in the trajectory. Parameters      what : str Requested particle property.  what must be a particle property or an alias. The following particle aliases are accepted: - 'position': 'particle.position' - 'pos': 'particle.position' - 'position[0]': 'particle.position[0]', - 'pos[0]': 'particle.position[0]' - 'x': 'particle.position[0]' - 'position[1]': 'particle.position[1]', - 'pos[1]': 'particle.position[1]' - 'y': 'particle.position[1]' - 'position[2]': 'particle.position[2]' - 'pos[2]': 'particle.position[2]' - 'z': 'particle.position[2]' - 'species': 'particle.species' - 'spe': 'particle.species' - 'label': 'particle.label' - 'index': 'particle.index' - 'mass': 'particle.mass' - 'radius': 'particle.radius' Returns    - to_dump : list List of the requested particle property with length equal to the number of frames in the trajectory. Each element of the list is a numpy.ndarray of the requested particle property. Examples     >>> traj = Trajectory('trajectory.xyz') >>> D = StructuralDescriptor(traj) >>> D.get_group_property('position', 0) >>> D.get_group_property('x', 1) >>> D.get_group_property('energy', 0)",
"func":1
},
{
"ref":"partycls.descriptor.descriptor.AngularStructuralDescriptor.dump",
"url":9,
"doc":"Alias for the method get_group_property.",
"func":1
},
{
"ref":"partycls.descriptor.descriptor.AngularStructuralDescriptor.group_fraction",
"url":9,
"doc":"Return the fraction of particles inside  group over the whole trajectory.",
"func":1
},
{
"ref":"partycls.descriptor.descriptor.AngularStructuralDescriptor.size",
"url":9,
"doc":"Total number of particles in the descriptor (i.e. in group=0)."
},
{
"ref":"partycls.descriptor.descriptor.AngularStructuralDescriptor.n_features",
"url":9,
"doc":"Number of features of the descriptor."
},
{
"ref":"partycls.descriptor.descriptor.AngularStructuralDescriptor.average",
"url":9,
"doc":"Average feature vector of the descriptor."
},
{
"ref":"partycls.descriptor.descriptor.AngularStructuralDescriptor.normalize",
"url":9,
"doc":"Generic normalization function for child classes. Returns the input distribution unchanged.",
"func":1
},
{
"ref":"partycls.descriptor.descriptor.DummyDescriptor",
"url":9,
"doc":"Base class for structural descriptors. The descriptor is calculated for the provided trajectory  trajectory . This can be: - an object implementing the  Trajectory interface ; - the path to a trajectory file in a format recognized by partyclsl ; A structural descriptor S(x) is a collection of N individual empirical correlation functions {s_i(x)} at the particle level, defined over a grid {x_j} of M features. These are stored in the  features array as a matrix usually refered to as the \"data set\": s_0(x_0) s_0(x_1)  . s_0(x_M) s_1(x_0) s_1(x_1)  . s_1(x_M)  .  .  . s_N(x_0) s_N(x_1)  . s_N(x_M) The  features array is None by default and is computed only when the  compute() method is called. The correlations can be calculated between two arbitrary subsets of particles called \"groups\": - group 0 is the main group, i.e. particles for which the correlations are being calculated ; - group 1 is the secondary group, i.e. particles that are being considered when calculating the correlations ; These groups are formed by adding filters on particles' properties (species, radius, position, etc.). Parameters      trajectory : str or an instance of  Trajectory . Trajectory on which the structural descriptor will be computed. Attributes      trajectory : Trajectory Trajectory on which the structural descriptor will be computed. active_filters : list of str All the active filters on both groups prior to the computation of the descriptor. dimension : int Spatial dimension of the descriptor (2 or 3). grid : array Grid over which the structural features will be computed. features : ndarray Array of all the structural features for the particles in group=0 in accordance with the defined filters (if any). This attribute is initialized when the method  compute is called (default value is None). groups : tuple Composition of the groups: groups[0] and groups[1] contain lists of all the  Particle objects in groups 0 and 1 respectively. Each element of the tuple is a list of  Particle in  trajectory , e.g. groups[0][0] is the list of all the particles in the first frame of  trajectory that belong to group=0. Examples:     - >>> D = StructuralDescriptor('trajectory.xyz') >>> D.add_filter(\"species  'A'\", group=0) >>> D.add_filter(\"species  'B'\", group=1) >>> D.active_filters [(\"particle.species  'A'\", 0), (\"particle.species  'B'\", 1)] >>> D.clear_filters(0) >>> D.active_filters [(\"particle.species  'B'\", 1)]"
},
{
"ref":"partycls.descriptor.descriptor.DummyDescriptor.name",
"url":9,
"doc":""
},
{
"ref":"partycls.descriptor.descriptor.DummyDescriptor.symbol",
"url":9,
"doc":""
},
{
"ref":"partycls.descriptor.descriptor.DummyDescriptor.normalize",
"url":9,
"doc":"Generic normalization function for child classes. Returns the input distribution unchanged.",
"func":1
},
{
"ref":"partycls.descriptor.descriptor.DummyDescriptor.add_filter",
"url":9,
"doc":"Add a filter on the group (0 or 1) to select the subset of particles that respects the provided condition. Parameters      condition : str The condition should have the following format:  _operator_  where: -  is a particle property (accepts aliases) ; - _operator_ is a logical operator ( =, >) ; -  is the corresponding value of  with the proper type ; group : int, optional Index of the group to which the filter must be applied. The default is 0. Returns    - None. Examples:     - >>> S = StructuralDescriptor('trajectory.xyz') >>> S.add_filter(\"particle.radius  >> S.add_filter(\"species  'A'\", group=1) >>> S.add_filter(\"x < 0\", group=0)  particles on the left side of the box",
"func":1
},
{
"ref":"partycls.descriptor.descriptor.DummyDescriptor.clear_filters",
"url":9,
"doc":"Clear all active filters on  group . All particles are included again in  group . Parameters      group : int, optional Index of the group on which to clear the filters. The default is 0. Returns    - None.",
"func":1
},
{
"ref":"partycls.descriptor.descriptor.DummyDescriptor.clear_all_filters",
"url":9,
"doc":"Clear all active filters in both groups. All particles are included in both groups again. Returns    - None.",
"func":1
},
{
"ref":"partycls.descriptor.descriptor.DummyDescriptor.group_size",
"url":9,
"doc":"Return the number of particles in  group .",
"func":1
},
{
"ref":"partycls.descriptor.descriptor.DummyDescriptor.get_group_property",
"url":9,
"doc":"Return a list of numpy arrays with the properties of the particles in group  group . The list size is the number of systems in the trajectory. Parameters      what : str Requested particle property.  what must be a particle property or an alias. The following particle aliases are accepted: - 'position': 'particle.position' - 'pos': 'particle.position' - 'position[0]': 'particle.position[0]', - 'pos[0]': 'particle.position[0]' - 'x': 'particle.position[0]' - 'position[1]': 'particle.position[1]', - 'pos[1]': 'particle.position[1]' - 'y': 'particle.position[1]' - 'position[2]': 'particle.position[2]' - 'pos[2]': 'particle.position[2]' - 'z': 'particle.position[2]' - 'species': 'particle.species' - 'spe': 'particle.species' - 'label': 'particle.label' - 'index': 'particle.index' - 'mass': 'particle.mass' - 'radius': 'particle.radius' Returns    - to_dump : list List of the requested particle property with length equal to the number of frames in the trajectory. Each element of the list is a numpy.ndarray of the requested particle property. Examples     >>> traj = Trajectory('trajectory.xyz') >>> D = StructuralDescriptor(traj) >>> D.get_group_property('position', 0) >>> D.get_group_property('x', 1) >>> D.get_group_property('energy', 0)",
"func":1
},
{
"ref":"partycls.descriptor.descriptor.DummyDescriptor.dump",
"url":9,
"doc":"Alias for the method get_group_property.",
"func":1
},
{
"ref":"partycls.descriptor.descriptor.DummyDescriptor.group_fraction",
"url":9,
"doc":"Return the fraction of particles inside  group over the whole trajectory.",
"func":1
},
{
"ref":"partycls.descriptor.descriptor.DummyDescriptor.size",
"url":9,
"doc":"Total number of particles in the descriptor (i.e. in group=0)."
},
{
"ref":"partycls.descriptor.descriptor.DummyDescriptor.n_features",
"url":9,
"doc":"Number of features of the descriptor."
},
{
"ref":"partycls.descriptor.descriptor.DummyDescriptor.average",
"url":9,
"doc":"Average feature vector of the descriptor."
},
{
"ref":"partycls.descriptor.ba",
"url":11,
"doc":""
},
{
"ref":"partycls.descriptor.ba.BondAngleDescriptor",
"url":11,
"doc":"Structural descriptor based on bond angles between particles. See the parent class for more details. Parameters      trajectory : str or an instance of  Trajectory . Trajectory on which the structural descriptor will be computed. dtheta : float Bin width in degrees. Attributes      trajectory : Trajectory Trajectory on which the structural descriptor will be computed. active_filters : list of str All the active filters on both groups prior to the computation of the descriptor. dimension : int Spatial dimension of the descriptor (2 or 3). grid : array Grid over which the structural features will be computed. features : ndarray Array of all the structural features for the particles in group=0 in accordance with the defined filters (if any). This attribute is initialized when the method  compute is called (default value is None). cutoffs : list of float List of cutoff distances to identify the nearest neighbors using the fixed-cutoff ('FC') method. nearest_neighbors_method : str, default: 'FC' Nearest neighbor method, 'FC' or 'SANN'. Examples:     - >>> D = BondAngleDescriptor('trajectory.xyz', dtheta=2.0) >>> D.nearest_neighbors_method = 'SANN' >>> D.add_filter(\"species  'A'\") >>> D.compute()"
},
{
"ref":"partycls.descriptor.ba.BondAngleDescriptor.name",
"url":11,
"doc":""
},
{
"ref":"partycls.descriptor.ba.BondAngleDescriptor.symbol",
"url":11,
"doc":""
},
{
"ref":"partycls.descriptor.ba.BondAngleDescriptor.dtheta",
"url":11,
"doc":""
},
{
"ref":"partycls.descriptor.ba.BondAngleDescriptor.compute",
"url":11,
"doc":"",
"func":1
},
{
"ref":"partycls.descriptor.ba.BondAngleDescriptor.normalize",
"url":11,
"doc":"Normalize a bond angle distribution. Parameters      distribution : array Distribution to normalize. method : str, optional Normalization method: - method='sin': by construction, the probability density of has a sinusoidal enveloppe in 3D for uniformly distributed points on a sphere (default) ; - method='pdf' : gives a flat probability density for uniformly distributed points on a sphere ; Raises    ValueError If  method is invalid. Returns    - array Normalized distribution.",
"func":1
},
{
"ref":"partycls.descriptor.ba.BondAngleDescriptor.set_cutoff",
"url":9,
"doc":"Set the nearest-neighbor cutoff for the pair of species (s1, s2). The cutoff of the mirror pair (s2, s1) is set automatically if the  mirror parameter is True (default). Parameters      s1 : str Symbol of the first species. s2 : str Symbol of the second species. rcut : float Value of the cutoff for the pair (s1,s2). mirror : bool, optional Set the cutoff for the mirror pair (s2,s1). The default is True. Returns    - None.",
"func":1
},
{
"ref":"partycls.descriptor.ba.BondAngleDescriptor.nearest_neighbors",
"url":9,
"doc":"Compute the nearest neighbors of particles in group=0 using one of the following methods: - \"Fixed cutoff\" (method='FC'): uses the partial radial distribution functions to compute the cutoffs between each possible pair of species (s1, s2) ; - \"Solid-Angle based Nearest Neighbors\" (method='SANN'): see van Meel et al. (https: doi.org/10.1063/1.4729313) ; Parameters      method : str, optional Method to identify nearest neighbors. Must be 'FC' or 'SANN'. The default is 'FC'. Returns    - None.",
"func":1
},
{
"ref":"partycls.descriptor.ba.BondAngleDescriptor.add_filter",
"url":9,
"doc":"Add a filter on the group (0 or 1) to select the subset of particles that respects the provided condition. Parameters      condition : str The condition should have the following format:  _operator_  where: -  is a particle property (accepts aliases) ; - _operator_ is a logical operator ( =, >) ; -  is the corresponding value of  with the proper type ; group : int, optional Index of the group to which the filter must be applied. The default is 0. Returns    - None. Examples:     - >>> S = StructuralDescriptor('trajectory.xyz') >>> S.add_filter(\"particle.radius  >> S.add_filter(\"species  'A'\", group=1) >>> S.add_filter(\"x < 0\", group=0)  particles on the left side of the box",
"func":1
},
{
"ref":"partycls.descriptor.ba.BondAngleDescriptor.clear_filters",
"url":9,
"doc":"Clear all active filters on  group . All particles are included again in  group . Parameters      group : int, optional Index of the group on which to clear the filters. The default is 0. Returns    - None.",
"func":1
},
{
"ref":"partycls.descriptor.ba.BondAngleDescriptor.clear_all_filters",
"url":9,
"doc":"Clear all active filters in both groups. All particles are included in both groups again. Returns    - None.",
"func":1
},
{
"ref":"partycls.descriptor.ba.BondAngleDescriptor.group_size",
"url":9,
"doc":"Return the number of particles in  group .",
"func":1
},
{
"ref":"partycls.descriptor.ba.BondAngleDescriptor.get_group_property",
"url":9,
"doc":"Return a list of numpy arrays with the properties of the particles in group  group . The list size is the number of systems in the trajectory. Parameters      what : str Requested particle property.  what must be a particle property or an alias. The following particle aliases are accepted: - 'position': 'particle.position' - 'pos': 'particle.position' - 'position[0]': 'particle.position[0]', - 'pos[0]': 'particle.position[0]' - 'x': 'particle.position[0]' - 'position[1]': 'particle.position[1]', - 'pos[1]': 'particle.position[1]' - 'y': 'particle.position[1]' - 'position[2]': 'particle.position[2]' - 'pos[2]': 'particle.position[2]' - 'z': 'particle.position[2]' - 'species': 'particle.species' - 'spe': 'particle.species' - 'label': 'particle.label' - 'index': 'particle.index' - 'mass': 'particle.mass' - 'radius': 'particle.radius' Returns    - to_dump : list List of the requested particle property with length equal to the number of frames in the trajectory. Each element of the list is a numpy.ndarray of the requested particle property. Examples     >>> traj = Trajectory('trajectory.xyz') >>> D = StructuralDescriptor(traj) >>> D.get_group_property('position', 0) >>> D.get_group_property('x', 1) >>> D.get_group_property('energy', 0)",
"func":1
},
{
"ref":"partycls.descriptor.ba.BondAngleDescriptor.dump",
"url":9,
"doc":"Alias for the method get_group_property.",
"func":1
},
{
"ref":"partycls.descriptor.ba.BondAngleDescriptor.group_fraction",
"url":9,
"doc":"Return the fraction of particles inside  group over the whole trajectory.",
"func":1
},
{
"ref":"partycls.descriptor.ba.BondAngleDescriptor.size",
"url":9,
"doc":"Total number of particles in the descriptor (i.e. in group=0)."
},
{
"ref":"partycls.descriptor.ba.BondAngleDescriptor.n_features",
"url":9,
"doc":"Number of features of the descriptor."
},
{
"ref":"partycls.descriptor.ba.BondAngleDescriptor.average",
"url":9,
"doc":"Average feature vector of the descriptor."
},
{
"ref":"partycls.descriptor.realspace_wrap",
"url":12,
"doc":"This module 'realspace_wrap' is auto-generated with f2py (version:2). Functions: Fortran 90/95 modules: compute  - pi,pbc(),radial_histogram(),angular_histogram(),pbc_(),cartesian_to_spherical(),factorial(),plm(),ylm(),qlm(),rotational_invariant(),ql(),qbarlm(),qbarl(),find_cutoff(),smoothed_qlm(),smoothed_ql(),nearest_neighbors(),sann(),sort()."
},
{
"ref":"partycls.descriptor.dscribe",
"url":13,
"doc":""
},
{
"ref":"partycls.descriptor.dscribe.DscribeDescriptor",
"url":13,
"doc":"Adapter for generic DScribe descriptors, without chemical species information. Essentially, all the particles are considered as hydrogen atoms."
},
{
"ref":"partycls.descriptor.dscribe.DscribeDescriptor.compute",
"url":13,
"doc":"",
"func":1
},
{
"ref":"partycls.descriptor.dscribe.DscribeDescriptor.normalize",
"url":13,
"doc":"Generic normalization function for child classes. Returns the input distribution unchanged.",
"func":1
},
{
"ref":"partycls.descriptor.dscribe.DscribeDescriptor.add_filter",
"url":9,
"doc":"Add a filter on the group (0 or 1) to select the subset of particles that respects the provided condition. Parameters      condition : str The condition should have the following format:  _operator_  where: -  is a particle property (accepts aliases) ; - _operator_ is a logical operator ( =, >) ; -  is the corresponding value of  with the proper type ; group : int, optional Index of the group to which the filter must be applied. The default is 0. Returns    - None. Examples:     - >>> S = StructuralDescriptor('trajectory.xyz') >>> S.add_filter(\"particle.radius  >> S.add_filter(\"species  'A'\", group=1) >>> S.add_filter(\"x < 0\", group=0)  particles on the left side of the box",
"func":1
},
{
"ref":"partycls.descriptor.dscribe.DscribeDescriptor.clear_filters",
"url":9,
"doc":"Clear all active filters on  group . All particles are included again in  group . Parameters      group : int, optional Index of the group on which to clear the filters. The default is 0. Returns    - None.",
"func":1
},
{
"ref":"partycls.descriptor.dscribe.DscribeDescriptor.clear_all_filters",
"url":9,
"doc":"Clear all active filters in both groups. All particles are included in both groups again. Returns    - None.",
"func":1
},
{
"ref":"partycls.descriptor.dscribe.DscribeDescriptor.group_size",
"url":9,
"doc":"Return the number of particles in  group .",
"func":1
},
{
"ref":"partycls.descriptor.dscribe.DscribeDescriptor.get_group_property",
"url":9,
"doc":"Return a list of numpy arrays with the properties of the particles in group  group . The list size is the number of systems in the trajectory. Parameters      what : str Requested particle property.  what must be a particle property or an alias. The following particle aliases are accepted: - 'position': 'particle.position' - 'pos': 'particle.position' - 'position[0]': 'particle.position[0]', - 'pos[0]': 'particle.position[0]' - 'x': 'particle.position[0]' - 'position[1]': 'particle.position[1]', - 'pos[1]': 'particle.position[1]' - 'y': 'particle.position[1]' - 'position[2]': 'particle.position[2]' - 'pos[2]': 'particle.position[2]' - 'z': 'particle.position[2]' - 'species': 'particle.species' - 'spe': 'particle.species' - 'label': 'particle.label' - 'index': 'particle.index' - 'mass': 'particle.mass' - 'radius': 'particle.radius' Returns    - to_dump : list List of the requested particle property with length equal to the number of frames in the trajectory. Each element of the list is a numpy.ndarray of the requested particle property. Examples     >>> traj = Trajectory('trajectory.xyz') >>> D = StructuralDescriptor(traj) >>> D.get_group_property('position', 0) >>> D.get_group_property('x', 1) >>> D.get_group_property('energy', 0)",
"func":1
},
{
"ref":"partycls.descriptor.dscribe.DscribeDescriptor.dump",
"url":9,
"doc":"Alias for the method get_group_property.",
"func":1
},
{
"ref":"partycls.descriptor.dscribe.DscribeDescriptor.group_fraction",
"url":9,
"doc":"Return the fraction of particles inside  group over the whole trajectory.",
"func":1
},
{
"ref":"partycls.descriptor.dscribe.DscribeDescriptor.size",
"url":9,
"doc":"Total number of particles in the descriptor (i.e. in group=0)."
},
{
"ref":"partycls.descriptor.dscribe.DscribeDescriptor.n_features",
"url":9,
"doc":"Number of features of the descriptor."
},
{
"ref":"partycls.descriptor.dscribe.DscribeDescriptor.average",
"url":9,
"doc":"Average feature vector of the descriptor."
},
{
"ref":"partycls.descriptor.dscribe.DscribeChemicalDescriptor",
"url":13,
"doc":"Adapter for generic DScribe descriptors, with chemical species information."
},
{
"ref":"partycls.descriptor.dscribe.DscribeChemicalDescriptor.normalize",
"url":13,
"doc":"Generic normalization function for child classes. Returns the input distribution unchanged.",
"func":1
},
{
"ref":"partycls.descriptor.dscribe.DscribeChemicalDescriptor.add_filter",
"url":9,
"doc":"Add a filter on the group (0 or 1) to select the subset of particles that respects the provided condition. Parameters      condition : str The condition should have the following format:  _operator_  where: -  is a particle property (accepts aliases) ; - _operator_ is a logical operator ( =, >) ; -  is the corresponding value of  with the proper type ; group : int, optional Index of the group to which the filter must be applied. The default is 0. Returns    - None. Examples:     - >>> S = StructuralDescriptor('trajectory.xyz') >>> S.add_filter(\"particle.radius  >> S.add_filter(\"species  'A'\", group=1) >>> S.add_filter(\"x < 0\", group=0)  particles on the left side of the box",
"func":1
},
{
"ref":"partycls.descriptor.dscribe.DscribeChemicalDescriptor.clear_filters",
"url":9,
"doc":"Clear all active filters on  group . All particles are included again in  group . Parameters      group : int, optional Index of the group on which to clear the filters. The default is 0. Returns    - None.",
"func":1
},
{
"ref":"partycls.descriptor.dscribe.DscribeChemicalDescriptor.clear_all_filters",
"url":9,
"doc":"Clear all active filters in both groups. All particles are included in both groups again. Returns    - None.",
"func":1
},
{
"ref":"partycls.descriptor.dscribe.DscribeChemicalDescriptor.group_size",
"url":9,
"doc":"Return the number of particles in  group .",
"func":1
},
{
"ref":"partycls.descriptor.dscribe.DscribeChemicalDescriptor.get_group_property",
"url":9,
"doc":"Return a list of numpy arrays with the properties of the particles in group  group . The list size is the number of systems in the trajectory. Parameters      what : str Requested particle property.  what must be a particle property or an alias. The following particle aliases are accepted: - 'position': 'particle.position' - 'pos': 'particle.position' - 'position[0]': 'particle.position[0]', - 'pos[0]': 'particle.position[0]' - 'x': 'particle.position[0]' - 'position[1]': 'particle.position[1]', - 'pos[1]': 'particle.position[1]' - 'y': 'particle.position[1]' - 'position[2]': 'particle.position[2]' - 'pos[2]': 'particle.position[2]' - 'z': 'particle.position[2]' - 'species': 'particle.species' - 'spe': 'particle.species' - 'label': 'particle.label' - 'index': 'particle.index' - 'mass': 'particle.mass' - 'radius': 'particle.radius' Returns    - to_dump : list List of the requested particle property with length equal to the number of frames in the trajectory. Each element of the list is a numpy.ndarray of the requested particle property. Examples     >>> traj = Trajectory('trajectory.xyz') >>> D = StructuralDescriptor(traj) >>> D.get_group_property('position', 0) >>> D.get_group_property('x', 1) >>> D.get_group_property('energy', 0)",
"func":1
},
{
"ref":"partycls.descriptor.dscribe.DscribeChemicalDescriptor.dump",
"url":9,
"doc":"Alias for the method get_group_property.",
"func":1
},
{
"ref":"partycls.descriptor.dscribe.DscribeChemicalDescriptor.group_fraction",
"url":9,
"doc":"Return the fraction of particles inside  group over the whole trajectory.",
"func":1
},
{
"ref":"partycls.descriptor.dscribe.DscribeChemicalDescriptor.size",
"url":9,
"doc":"Total number of particles in the descriptor (i.e. in group=0)."
},
{
"ref":"partycls.descriptor.dscribe.DscribeChemicalDescriptor.n_features",
"url":9,
"doc":"Number of features of the descriptor."
},
{
"ref":"partycls.descriptor.dscribe.DscribeChemicalDescriptor.average",
"url":9,
"doc":"Average feature vector of the descriptor."
},
{
"ref":"partycls.trajectory",
"url":14,
"doc":"Physical trajectory. This class is inspired by the  atooms framework authored by Daniele Coslovich See https: framagit.org/atooms/atooms"
},
{
"ref":"partycls.trajectory.Trajectory",
"url":14,
"doc":"A trajectory is composed by one or several frames, each frame being an instance of  System . Trajectory instances are iterable. By default, only the positions and particle types are being read from the trajectory file. Additional particle properties in the file can be read using the  additional_fields parameter. Parameters      filename : str Path to the trajectory file to read. fmt : str, optional, default: 'xyz' Format of the trajectory. Needed when using \"atooms\" as a backend. backend : str, optional, default: None Name of a third-party package to use as backend when reading the input trajectory. Currently supports \"atooms\" and \"mdtraj\". top : str, mdtraj.Trajectory, or mdtraj.Topology, optional, defaut: None Topology information. Needed when using \"mdtraj\" as backend on a trajectory file whose format requires topology information. See MDTraj documentation for more information. additional_fields : list of str, optional, default: [] Additional fields (i.e. particle properties) to read from the trajectory. Not all trajectory formats allow for additional fields. first : int, optional, default: 0 Index of the first frame to consider in the trajectory. Starts at zero. last : int, optional, default: None Index of the last frame to consider in the trajectory. Default is the last frame. step : int, optional, default: 1 Step between each frame to consider in the trajectory. For example, if  step=2 , one every two frames is read. Attributes      filename : str Name of the original trajectory file. fmt : str Format of the original trajectory file. backend : str or None Name of the third-party package used to read the input trajectory file. additional_fields : list, default: None List of additional particle properties that were extracted from the original trajectory file. Examples     >>> from partycls.trajectory import Trajectory >>> traj = Trajectory('trajectory.xyz', additional_fields=['mass'])"
},
{
"ref":"partycls.trajectory.Trajectory.remove",
"url":14,
"doc":"Remove the system at position  frame from the trajectory. Parameters      frame : int Index of the frame to remove from the trajectory. Returns    - None.",
"func":1
},
{
"ref":"partycls.trajectory.Trajectory.get_property",
"url":14,
"doc":"Return a list of numpy arrays with the system property specified by  what . The list size is the number of systems in the trajectory. Parameters      what : str Requested system property.  what must be of the form \"particle. \" or \"cell. \". The following particle aliases are accepted: - 'position': 'particle.position' - 'pos': 'particle.position' - 'position[0]': 'particle.position[0]', - 'pos[0]': 'particle.position[0]' - 'x': 'particle.position[0]' - 'position[1]': 'particle.position[1]', - 'pos[1]': 'particle.position[1]' - 'y': 'particle.position[1]' - 'position[2]': 'particle.position[2]' - 'pos[2]': 'particle.position[2]' - 'z': 'particle.position[2]' - 'species': 'particle.species' - 'spe': 'particle.species' - 'label': 'particle.label' - 'index': 'particle.index' - 'mass': 'particle.mass' - 'radius': 'particle.radius' subset : str, optional Subset of particles for which the property must be dumped. Must be of the form \"particle. \" unless \" \" is an alias. The default is None (all particles will be included). This is ignored if  what is cell property. Returns    - to_dump : list List of the requested system property with length equal to the number of frames in the trajectory. Each element of the list is a numpy.ndarray of the requested system property. Examples     >>> traj = Trajectory('trajectory.xyz') >>> pos = traj.get_property('position') >>> spe = traj.get_property('species') >>> sides = traj.get_property('cell.side')",
"func":1
},
{
"ref":"partycls.trajectory.Trajectory.dump",
"url":14,
"doc":"Alias for the method get_property().",
"func":1
},
{
"ref":"partycls.trajectory.Trajectory.set_property",
"url":14,
"doc":"Set a property  what to  value for all the particles in the trajectory or for a given subset of particles specified by  subset . Parameters      what : str Name of the property to set. This is considered to be a particle property by default, unless it starts with \"cell\", e.g. \"cell.side\". value : int, float, list, or numpy.ndarray Value(s) of the property to set. An instance of  int or  float will set the same value for all concerned particles. An instance of  list or  numpy.ndarray will assign a specific value to each particle. In this case, the shape of  value should respect the number of frames in the trajectory and the number of concerned particles. subset : str, optional Particles to which the property must be set. The default is None. This is ignored if  what is cell property. Returns    - None. Examples     >>> traj.set_property('mass', 1.0) >>> traj.set_property('radius', 0.5, \"species  'A'\") >>> labels =  0, 1, 0],  2 frames, 3 particles in the subset [1, 1, 0 >>> traj.set_property('label', labels, \"species  'B'\")",
"func":1
},
{
"ref":"partycls.trajectory.Trajectory.show",
"url":14,
"doc":"Show the frames on index  frames of the trajectory and color particles according to an arbitrary property, such as species, cluster label, etc. Current visualization backends are 'matplotlib', 'ovito', and '3dmol'. Parameters      frames : list of int, optional Indices of the frames to show. The default is None (shows all frames). backend : str, optional Name of the backend to use for visualization. The default is 'matplotlib'. color : str, optional Name of the particle property to use as basis for coloring the particles. This property must be defined for all the particles in the system. The default is 'species'.  kwargs : additional keyworded arguments (backend-dependent). Raises    ValueError In case of unknown  backend . Returns    - list of Figure or View (backend-dependent) Examples     >>> traj.show(frames=[0,1,2], color='label', backend='3dmol') >>> traj.show(frames=[0,1], color='energy', backend='matplotlib', cmap='viridis') >>> traj[0].show()  use the iterability of Trajectory objects",
"func":1
},
{
"ref":"partycls.trajectory.Trajectory.fold",
"url":14,
"doc":"Fold the particle positions into the central cell. Returns    - None.",
"func":1
},
{
"ref":"partycls.trajectory.Trajectory.write",
"url":14,
"doc":"Write the current trajectory to a file. Parameters      output_path : str Name of the output trajectory file. fmt : str, optional Format of the output trajectory file. The default is 'xyz'. backend : str, optional Name of a third-party package to use when writing the output trajectory. The default is None. additional_fields : list of str, optional Additional fields (i.e. particle properties) to write in the output trajectory. Not all trajectory formats allow for additional fields. The default is []. precision : int, optional Number of decimals when writing the output trajectory. The default is 6. Raises    ValueError - If  backend=None and  fmt is not recognized natively. - If  backend is unknown. Returns    - None.",
"func":1
},
{
"ref":"partycls.feature_scaling",
"url":15,
"doc":"Feature scaling techniques, to be performed on a dataset stored in a numpy array."
},
{
"ref":"partycls.feature_scaling.ZScore",
"url":15,
"doc":"Standardize features by removing the mean and scaling to unit variance The standard score of a sample  x is calculated as: z = (x - u) / s where  u is the mean of the training samples or zero if  with_mean=False , and  s is the standard deviation of the training samples or one if  with_std=False . Centering and scaling happen independently on each feature by computing the relevant statistics on the samples in the training set. Mean and standard deviation are then stored to be used on later data using :meth: transform . Standardization of a dataset is a common requirement for many machine learning estimators: they might behave badly if the individual features do not more or less look like standard normally distributed data (e.g. Gaussian with 0 mean and unit variance). For instance many elements used in the objective function of a learning algorithm (such as the RBF kernel of Support Vector Machines or the L1 and L2 regularizers of linear models) assume that all features are centered around 0 and have variance in the same order. If a feature has a variance that is orders of magnitude larger that others, it might dominate the objective function and make the estimator unable to learn from other features correctly as expected. This scaler can also be applied to sparse CSR or CSC matrices by passing  with_mean=False to avoid breaking the sparsity structure of the data. Read more in the :ref: User Guide   . Parameters      copy : bool, default=True If False, try to avoid a copy and do inplace scaling instead. This is not guaranteed to always work inplace; e.g. if the data is not a NumPy array or scipy.sparse CSR matrix, a copy may still be returned. with_mean : bool, default=True If True, center the data before scaling. This does not work (and will raise an exception) when attempted on sparse matrices, because centering them entails building a dense matrix which in common use cases is likely to be too large to fit in memory. with_std : bool, default=True If True, scale the data to unit variance (or equivalently, unit standard deviation). Attributes      scale_ : ndarray of shape (n_features,) or None Per feature relative scaling of the data to achieve zero mean and unit variance. Generally this is calculated using  np.sqrt(var_) . If a variance is zero, we can't achieve unit variance, and the data is left as-is, giving a scaling factor of 1.  scale_ is equal to  None when  with_std=False .  versionadded 0.17  scale_ mean_ : ndarray of shape (n_features,) or None The mean value for each feature in the training set. Equal to  None when  with_mean=False . var_ : ndarray of shape (n_features,) or None The variance for each feature in the training set. Used to compute  scale_ . Equal to  None when  with_std=False . n_samples_seen_ : int or ndarray of shape (n_features,) The number of samples processed by the estimator for each feature. If there are no missing samples, the  n_samples_seen will be an integer, otherwise it will be an array of dtype int. If  sample_weights are used it will be a float (if no missing data) or an array of dtype float that sums the weights seen so far. Will be reset on new calls to fit, but increments across  partial_fit calls. Examples     >>> from sklearn.preprocessing import StandardScaler >>> data =  0, 0], [0, 0], [1, 1], [1, 1 >>> scaler = StandardScaler() >>> print(scaler.fit(data StandardScaler() >>> print(scaler.mean_) [0.5 0.5] >>> print(scaler.transform(data  -1. -1.] [-1. -1.] [ 1. 1.] [ 1. 1. >>> print(scaler.transform( 2, 2   3. 3. See Also     scale : Equivalent function without the estimator API. :class: ~sklearn.decomposition.PCA : Further removes the linear correlation across features with 'whiten=True'. Notes   - NaNs are treated as missing values: disregarded in fit, and maintained in transform. We use a biased estimator for the standard deviation, equivalent to  numpy.std(x, ddof=0) . Note that the choice of  ddof is unlikely to affect model performance. For a comparison of the different scalers, transformers, and normalizers, see :ref: examples/preprocessing/plot_all_scaling.py   ."
},
{
"ref":"partycls.feature_scaling.ZScore.symbol",
"url":15,
"doc":""
},
{
"ref":"partycls.feature_scaling.ZScore.full_name",
"url":15,
"doc":""
},
{
"ref":"partycls.feature_scaling.ZScore.scale",
"url":15,
"doc":"Standardize features by removing the mean and scaling to unit variance. Parameters      X : numpy.ndarray Original features. Returns    - numpy.ndarray Scaled features.",
"func":1
},
{
"ref":"partycls.feature_scaling.MinMax",
"url":15,
"doc":"Transform features by scaling each feature to a given range. This estimator scales and translates each feature individually such that it is in the given range on the training set, e.g. between zero and one. The transformation is given by X_std = (X - X.min(axis=0 / (X.max(axis=0) - X.min(axis=0 X_scaled = X_std  (max - min) + min where min, max = feature_range. This transformation is often used as an alternative to zero mean, unit variance scaling. Read more in the :ref: User Guide   . Parameters      feature_range : tuple (min, max), default=(0, 1) Desired range of transformed data. copy : bool, default=True Set to False to perform inplace row normalization and avoid a copy (if the input is already a numpy array). clip: bool, default=False Set to True to clip transformed values of held-out data to provided  feature range .  versionadded 0.24 Attributes      min_ : ndarray of shape (n_features,) Per feature adjustment for minimum. Equivalent to  min - X.min(axis=0)  self.scale_ scale_ : ndarray of shape (n_features,) Per feature relative scaling of the data. Equivalent to  (max - min) / (X.max(axis=0) - X.min(axis=0   versionadded 0.17  scale_ attribute. data_min_ : ndarray of shape (n_features,) Per feature minimum seen in the data  versionadded 0.17  data_min_ data_max_ : ndarray of shape (n_features,) Per feature maximum seen in the data  versionadded 0.17  data_max_ data_range_ : ndarray of shape (n_features,) Per feature range  (data_max_ - data_min_) seen in the data  versionadded 0.17  data_range_ n_samples_seen_ : int The number of samples processed by the estimator. It will be reset on new calls to fit, but increments across  partial_fit calls. Examples     >>> from sklearn.preprocessing import MinMaxScaler >>> data =  -1, 2], [-0.5, 6], [0, 10], [1, 18 >>> scaler = MinMaxScaler() >>> print(scaler.fit(data MinMaxScaler() >>> print(scaler.data_max_) [ 1. 18.] >>> print(scaler.transform(data  0. 0. ] [0.25 0.25] [0.5 0.5 ] [1. 1.  >>> print(scaler.transform( 2, 2   1.5 0.  See Also     minmax_scale : Equivalent function without the estimator API. Notes   - NaNs are treated as missing values: disregarded in fit, and maintained in transform. For a comparison of the different scalers, transformers, and normalizers, see :ref: examples/preprocessing/plot_all_scaling.py   ."
},
{
"ref":"partycls.feature_scaling.MinMax.symbol",
"url":15,
"doc":""
},
{
"ref":"partycls.feature_scaling.MinMax.full_name",
"url":15,
"doc":""
},
{
"ref":"partycls.feature_scaling.MinMax.scale",
"url":15,
"doc":"Transform features by scaling each feature to a given range (default is [0,1]). Parameters      X : numpy.ndarray Original features. Returns    - numpy.ndarray Scaled features.",
"func":1
},
{
"ref":"partycls.feature_scaling.MaxAbs",
"url":15,
"doc":"Scale each feature by its maximum absolute value. This estimator scales and translates each feature individually such that the maximal absolute value of each feature in the training set will be 1.0. It does not shift/center the data, and thus does not destroy any sparsity. This scaler can also be applied to sparse CSR or CSC matrices.  versionadded 0.17 Parameters      copy : bool, default=True Set to False to perform inplace scaling and avoid a copy (if the input is already a numpy array). Attributes      scale_ : ndarray of shape (n_features,) Per feature relative scaling of the data.  versionadded 0.17  scale_ attribute. max_abs_ : ndarray of shape (n_features,) Per feature maximum absolute value. n_samples_seen_ : int The number of samples processed by the estimator. Will be reset on new calls to fit, but increments across  partial_fit calls. Examples     >>> from sklearn.preprocessing import MaxAbsScaler >>> X =  1., -1., 2.],  . [ 2., 0., 0.],  . [ 0., 1., -1. >>> transformer = MaxAbsScaler().fit(X) >>> transformer MaxAbsScaler() >>> transformer.transform(X) array( 0.5, -1. , 1. ], [ 1. , 0. , 0. ], [ 0. , 1. , -0.5 ) See Also     maxabs_scale : Equivalent function without the estimator API. Notes   - NaNs are treated as missing values: disregarded in fit, and maintained in transform. For a comparison of the different scalers, transformers, and normalizers, see :ref: examples/preprocessing/plot_all_scaling.py   ."
},
{
"ref":"partycls.feature_scaling.MaxAbs.symbol",
"url":15,
"doc":""
},
{
"ref":"partycls.feature_scaling.MaxAbs.full_name",
"url":15,
"doc":""
},
{
"ref":"partycls.feature_scaling.MaxAbs.scale",
"url":15,
"doc":"Scale each feature by its maximum absolute value. Parameters      X : numpy.ndarray Original features. Returns    - numpy.ndarray Scaled features.",
"func":1
},
{
"ref":"partycls.feature_scaling.Robust",
"url":15,
"doc":"Scale features using statistics that are robust to outliers. This Scaler removes the median and scales the data according to the quantile range (defaults to IQR: Interquartile Range). The IQR is the range between the 1st quartile (25th quantile) and the 3rd quartile (75th quantile). Centering and scaling happen independently on each feature by computing the relevant statistics on the samples in the training set. Median and interquartile range are then stored to be used on later data using the  transform method. Standardization of a dataset is a common requirement for many machine learning estimators. Typically this is done by removing the mean and scaling to unit variance. However, outliers can often influence the sample mean / variance in a negative way. In such cases, the median and the interquartile range often give better results.  versionadded 0.17 Read more in the :ref: User Guide   . Parameters      with_centering : bool, default=True If True, center the data before scaling. This will cause  transform to raise an exception when attempted on sparse matrices, because centering them entails building a dense matrix which in common use cases is likely to be too large to fit in memory. with_scaling : bool, default=True If True, scale the data to interquartile range. quantile_range : tuple (q_min, q_max), 0.0  >> from sklearn.preprocessing import RobustScaler >>> X =  1., -2., 2.],  . [ -2., 1., 3.],  . [ 4., 1., -2. >>> transformer = RobustScaler().fit(X) >>> transformer RobustScaler() >>> transformer.transform(X) array( 0. , -2. , 0. ], [-1. , 0. , 0.4], [ 1. , 0. , -1.6 ) See Also     robust_scale : Equivalent function without the estimator API. :class: ~sklearn.decomposition.PCA Further removes the linear correlation across features with 'whiten=True'. Notes   - For a comparison of the different scalers, transformers, and normalizers, see :ref: examples/preprocessing/plot_all_scaling.py   . https: en.wikipedia.org/wiki/Median https: en.wikipedia.org/wiki/Interquartile_range"
},
{
"ref":"partycls.feature_scaling.Robust.symbol",
"url":15,
"doc":""
},
{
"ref":"partycls.feature_scaling.Robust.full_name",
"url":15,
"doc":""
},
{
"ref":"partycls.feature_scaling.Robust.scale",
"url":15,
"doc":"Scale features using statistics that are robust to outliers. Parameters      X : numpy.ndarray Original features. Returns    - numpy.ndarray Scaled features.",
"func":1
},
{
"ref":"partycls.workflow",
"url":16,
"doc":"Workflow for clustering analysis. A workflow is a procedure that goes through various steps (some of which are optional) to perform a structural clustering on a trajectory."
},
{
"ref":"partycls.workflow.Workflow",
"url":16,
"doc":"A workflow is a clustering procedure that goes through the following steps: - compute a structural descriptor on a given trajectory ; - (optional) apply a feature scaling on the previously computed structural features ; - (optional) apply a dimensionality reduction on the (raw/scaled) features ; - run a clustering algorithm to partition particles into structurally different clusters ; Parameters      trajectory : Trajectory, or str An instance of  Trajectory a path to trajectory file to read, or an instance of a class with compatible interface. descriptor : str, or an instance of StructuralDescriptor Structural descriptor to be computed on the trajectory. See the  descriptor_db class attribute for compatible strings. Examples : 'gr' : radial distribution of particles around a central particle. 'ba' : angular distribution of pairs of nearest neighbors of a central particle. 'bo' : Steinhardt bond-orientational order parameter (see https: doi.org/10.1103/PhysRevB.28.784) 'ld' : Lechner-Dellago cond-orientational order parameter (see https: doi.org/10.1063/1.2977970) scaling : str, None, or an object with the proper interface, optional, default: None Feature scaling method. See the  scaling_db class attribute for compatible strings. Examples: 'zscore' : standardize features by removing the mean and scaling to unit variance 'minmax' : scale and translate each feature individually such that it is in the given range on the training set, e.g. between zero and one 'maxabs' : scale and translate each feature individually such that the maximal absolute value of each feature in the training set will be 1. 'robust' : remove the median and scale the data according to the specified quantile range (default is between 25th quantile and 75th quantile). dim_reduction : str, None, or an object with the proper interface, optional, default: None Dimensionality reduction method. See the  dim_reduction_db class attribute for compatible strings. Examples: 'pca' : Principal Component Analysis. 'tsne' : t-distributed Stochastic Neighbor Embedding. 'lle' : Locally Linear Embedding. 'ae' : neural network Auto-Encoder. clustering : str, an instance of  Clustering , or an object with the proper interface, optional, default: 'kmeans' Clustering algorithm. See the  clustering_db class attribute for compatible strings. Examples: 'kmeans' : K-Means algorithm. 'gmm' : Gaussian Mixture Model. 'cinf' : Community Inference (see https: doi.org/10.1063/5.0004732). Attributes      trajectory : Trajectory The trajectory file as read by the Trajectory class. descriptor : StructuralDescriptor Structural descriptor associated to the trajectory. scaling : ZScore, MinMax, MaxAbs, Robust Feature scaling. dim_reduction : PCA, TSNE, LocallyLinearEmbedding, AutoEncoder Dimensionality reduction. clustering : Clustering Clustering method. output_metadata : dict Dictionnary that controls the writing process and the properties of all the output files. features: numpy.ndarray Raw features as computed by the associated structural descriptor. Initial value is None if features were not computed. scaled_features: numpy.ndarray Features after being rescaled by a feature scaling method. Equal to None if no scaling is applied to the features. reduced_features: numpy.ndarray Features in the reduced space after applying a dimensionality reduction technique. Equal to None if no reduction is applied to the features. naming_convention : str Base name for output files. Default is '{filename}.{code}.{descriptor}.{clustering}', where each tag will be replaced by its value in the current instance of  Workflow (e.g. \"traj.xyz.partycls.gr.kmeans\"). Base name can be changed using any combination of the available tags: {filename}, {code}, {descriptor}, {scaling}, {dim_reduction}, {clustering}. Example: \"{filename}_descriptor-{descriptor}_scaling-{scaling}.{code}\". Examples     >>> from partycls import Workflow >>> wf = Workflow('trajectory.xyz', descriptor='ba', scaling='zscore') >>> wf.run()"
},
{
"ref":"partycls.workflow.Workflow.descriptor_db",
"url":16,
"doc":""
},
{
"ref":"partycls.workflow.Workflow.clustering_db",
"url":16,
"doc":""
},
{
"ref":"partycls.workflow.Workflow.scaling_db",
"url":16,
"doc":""
},
{
"ref":"partycls.workflow.Workflow.dim_reduction_db",
"url":16,
"doc":""
},
{
"ref":"partycls.workflow.Workflow.run",
"url":16,
"doc":"Compute the clustering and write the output files according to the defined workflow : - compute the descriptor ; - (optional) apply feature scaling ; - (optional) apply dimensionality reduction ; - compute the clustering ; - (optional) write the output files ; Raises    ValueError If a community inference clustering is attempted with feature scaling or dimensionality reduction. Returns    - None.",
"func":1
},
{
"ref":"partycls.workflow.Workflow.labels",
"url":16,
"doc":""
},
{
"ref":"partycls.workflow.Workflow.fractions",
"url":16,
"doc":""
},
{
"ref":"partycls.workflow.Workflow.populations",
"url":16,
"doc":""
},
{
"ref":"partycls.workflow.Workflow.centroids",
"url":16,
"doc":""
},
{
"ref":"partycls.workflow.Workflow.set_output_metadata",
"url":16,
"doc":"Change the output properties. Parameters      what : {'trajectory', 'log', 'centroids', 'labels', or 'dataset'} Type of output file to change  kwargs : keywords arguments (specific to each type of file) DESCRIPTION. Returns    - None. Examples     >>> wf = Workflow('trajectory.xyz') >>> wf.set_output_metadata('log', enable=False)  do not write the log file >>> wf.set_output_metadata('trajectory', filename='awesome_trajectory.xyz')  change the default output name >>> wf.run('dataset', enable=True, precision=8)  write the dataset and change the writing precision to 8 digits",
"func":1
},
{
"ref":"partycls.workflow.Workflow.disable_output",
"url":16,
"doc":"Disable all outputs. Returns    - None.",
"func":1
},
{
"ref":"partycls.workflow.Workflow.write_trajectory",
"url":16,
"doc":"Write the trajectory file with cluster labels (default) and other additional fields (if any). Parameters      filename : str, optional Filename of the output trajectory. Uses a default naming convention if not specified. The default is None. fmt : str, optional Output trajectory format. The default is 'xyz'. backend : str, optional Name of the backend to use to write the trajectory. Must be either 'atooms' or 'mdtraj'. The default is None. additional_fields : list, optional Additional fields (i.e. particle properties) to write in the output trajectory. Note that all the  Particle objects should have the specified properties as attributes. The default is []. precision : int, optional Number of decimals when writing the output trajectory. The default is 6. Returns    - None. Examples     >>> wf = Workflow('trajectory.xyz') >>> wf.write_trajectory(fmt='rumd') >>> wf.write_trajectory(additional_field=['particle.mass'])   Particle must have the attribute  mass . >>> wf.write_trajectory(filename='my_custom_name', precision=8)",
"func":1
},
{
"ref":"partycls.workflow.Workflow.write_log",
"url":16,
"doc":"Write a log file with all relevant information about the workflow. The log file can be written only if the workflow has been run at least once with the method  Workflow.run . Parameters      filename : str, optional Filename of the log file. Uses a default naming convention if not specified. The default is None. precision : int, optional Number of decimals when writing the log file. The default is 6.  kwargs : TYPE DESCRIPTION. Returns    - None.",
"func":1
},
{
"ref":"partycls.workflow.Workflow.write_centroids",
"url":16,
"doc":"Write the coordinates of the clusters' centroids using the raw features from the descriptor (i.e. nor scaled or reduced). Parameters      filename : str, optional Filename of the centroids file. Uses a default naming convention if not specified. The default is None. precision : int, optional Number of decimals when writing the centroids file. The default is 6. Returns    - None.",
"func":1
},
{
"ref":"partycls.workflow.Workflow.write_labels",
"url":16,
"doc":"Write the clusters' labels only. Parameters      filename : str, optional Filename of the labels file. Uses a default naming convention if not specified. The default is None. Returns    - None.",
"func":1
},
{
"ref":"partycls.workflow.Workflow.write_dataset",
"url":16,
"doc":"Write the full raw dataset from the descriptor as an array (i.e. all the individual raw features of each particle). Parameters      filename : str, optional Filename of the dataset file. Uses a default naming convention if not specified. The default is None. precision : int, optional Number of decimals when writing the dataset file. The default is 6. Returns    - None.",
"func":1
},
{
"ref":"partycls.cell",
"url":17,
"doc":"Simulation cell. This class is inspired by the  atooms framework authored by Daniele Coslovich See https: framagit.org/atooms/atooms"
},
{
"ref":"partycls.cell.Cell",
"url":17,
"doc":"Orthorhombic cell. Parameters      side : list of float or float array List of lengths for the sides of the cell. Attributes      side : float array List of lengths for the sides of the cell. The default is None (will set  True in each direction). periodic : bool array Periodicity of the cell on each axis. Examples     >>> c = Cell([2.0, 2.0, 2.0]) >>> c.volume 8.0"
},
{
"ref":"partycls.cell.Cell.volume",
"url":17,
"doc":"Volume of the cell."
},
{
"ref":"partycls.dim_reduction",
"url":18,
"doc":"Dimensionality reduction techniques (linear and non-linear), to be performed on a dataset stored in a numpy array."
},
{
"ref":"partycls.dim_reduction.PCA",
"url":18,
"doc":"Principal component analysis (PCA). Linear dimensionality reduction using Singular Value Decomposition of the data to project it to a lower dimensional space. The input data is centered but not scaled for each feature before applying the SVD. It uses the LAPACK implementation of the full SVD or a randomized truncated SVD by the method of Halko et al. 2009, depending on the shape of the input data and the number of components to extract. It can also use the scipy.sparse.linalg ARPACK implementation of the truncated SVD. Notice that this class does not support sparse input. See :class: TruncatedSVD for an alternative with sparse data. Read more in the :ref: User Guide   . Parameters      n_components : int, float or 'mle', default=None Number of components to keep. if n_components is not set all components are kept n_components  min(n_samples, n_features) If  n_components  'mle' and  svd_solver  'full' , Minka's MLE is used to guess the dimension. Use of  n_components  'mle' will interpret  svd_solver  'auto' as  svd_solver  'full' . If  0   .  versionadded 0.18.0 Attributes      components_ : ndarray of shape (n_components, n_features) Principal axes in feature space, representing the directions of maximum variance in the data. The components are sorted by  explained_variance_ . explained_variance_ : ndarray of shape (n_components,) The amount of variance explained by each of the selected components. Equal to n_components largest eigenvalues of the covariance matrix of X.  versionadded 0.18 explained_variance_ratio_ : ndarray of shape (n_components,) Percentage of variance explained by each of the selected components. If  n_components is not set then all components are stored and the sum of the ratios is equal to 1.0. singular_values_ : ndarray of shape (n_components,) The singular values corresponding to each of the selected components. The singular values are equal to the 2-norms of the  n_components variables in the lower-dimensional space.  versionadded 0.19 mean_ : ndarray of shape (n_features,) Per-feature empirical mean, estimated from the training set. Equal to  X.mean(axis=0) . n_components_ : int The estimated number of components. When n_components is set to 'mle' or a number between 0 and 1 (with svd_solver  'full') this number is estimated from input data. Otherwise it equals the parameter n_components, or the lesser value of n_features and n_samples if n_components is None. n_features_ : int Number of features in the training data. n_samples_ : int Number of samples in the training data. noise_variance_ : float The estimated noise covariance following the Probabilistic PCA model from Tipping and Bishop 1999. See \"Pattern Recognition and Machine Learning\" by C. Bishop, 12.2.1 p. 574 or http: www.miketipping.com/papers/met-mppca.pdf. It is required to compute the estimated data covariance and score samples. Equal to the average of (min(n_features, n_samples) - n_components) smallest eigenvalues of the covariance matrix of X. See Also     KernelPCA : Kernel Principal Component Analysis. SparsePCA : Sparse Principal Component Analysis. TruncatedSVD : Dimensionality reduction using truncated SVD. IncrementalPCA : Incremental Principal Component Analysis. References      For n_components  'mle', this class uses the method of  Minka, T. P. \"Automatic choice of dimensionality for PCA\". In NIPS, pp. 598-604 Implements the probabilistic PCA model from: Tipping, M. E., and Bishop, C. M. (1999). \"Probabilistic principal component analysis\". Journal of the Royal Statistical Society: Series B (Statistical Methodology), 61(3), 611-622. via the score and score_samples methods. See http: www.miketipping.com/papers/met-mppca.pdf For svd_solver  'arpack', refer to  scipy.sparse.linalg.svds . For svd_solver  'randomized', see:  Halko, N., Martinsson, P. G., and Tropp, J. A. (2011). \"Finding structure with randomness: Probabilistic algorithms for constructing approximate matrix decompositions\". SIAM review, 53(2), 217-288. and also  Martinsson, P. G., Rokhlin, V., and Tygert, M. (2011). \"A randomized algorithm for the decomposition of matrices\". Applied and Computational Harmonic Analysis, 30(1), 47-68. Examples     >>> import numpy as np >>> from sklearn.decomposition import PCA >>> X = np.array( -1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2 ) >>> pca = PCA(n_components=2) >>> pca.fit(X) PCA(n_components=2) >>> print(pca.explained_variance_ratio_) [0.9924 . 0.0075 .] >>> print(pca.singular_values_) [6.30061 . 0.54980 .] >>> pca = PCA(n_components=2, svd_solver='full') >>> pca.fit(X) PCA(n_components=2, svd_solver='full') >>> print(pca.explained_variance_ratio_) [0.9924 . 0.00755 .] >>> print(pca.singular_values_) [6.30061 . 0.54980 .] >>> pca = PCA(n_components=1, svd_solver='arpack') >>> pca.fit(X) PCA(n_components=1, svd_solver='arpack') >>> print(pca.explained_variance_ratio_) [0.99244 .] >>> print(pca.singular_values_) [6.30061 .]"
},
{
"ref":"partycls.dim_reduction.PCA.symbol",
"url":18,
"doc":""
},
{
"ref":"partycls.dim_reduction.PCA.full_name",
"url":18,
"doc":""
},
{
"ref":"partycls.dim_reduction.PCA.reduce",
"url":18,
"doc":"Project the input features onto a reduced space using principal component analysis. Parameters      X : numpy.ndarray Features in the original space. Returns    - numpy.ndarray Features in the reduced space.",
"func":1
},
{
"ref":"partycls.dim_reduction.TSNE",
"url":18,
"doc":"t-distributed Stochastic Neighbor Embedding. t-SNE [1] is a tool to visualize high-dimensional data. It converts similarities between data points to joint probabilities and tries to minimize the Kullback-Leibler divergence between the joint probabilities of the low-dimensional embedding and the high-dimensional data. t-SNE has a cost function that is not convex, i.e. with different initializations we can get different results. It is highly recommended to use another dimensionality reduction method (e.g. PCA for dense data or TruncatedSVD for sparse data) to reduce the number of dimensions to a reasonable amount (e.g. 50) if the number of features is very high. This will suppress some noise and speed up the computation of pairwise distances between samples. For more tips see Laurens van der Maaten's FAQ [2]. Read more in the :ref: User Guide   . Parameters      n_components : int, default=2 Dimension of the embedded space. perplexity : float, default=30.0 The perplexity is related to the number of nearest neighbors that is used in other manifold learning algorithms. Larger datasets usually require a larger perplexity. Consider selecting a value between 5 and 50. Different values can result in significantly different results. early_exaggeration : float, default=12.0 Controls how tight natural clusters in the original space are in the embedded space and how much space will be between them. For larger values, the space between natural clusters will be larger in the embedded space. Again, the choice of this parameter is not very critical. If the cost function increases during initial optimization, the early exaggeration factor or the learning rate might be too high. learning_rate : float, default=200.0 The learning rate for t-SNE is usually in the range [10.0, 1000.0]. If the learning rate is too high, the data may look like a 'ball' with any point approximately equidistant from its nearest neighbours. If the learning rate is too low, most points may look compressed in a dense cloud with few outliers. If the cost function gets stuck in a bad local minimum increasing the learning rate may help. n_iter : int, default=1000 Maximum number of iterations for the optimization. Should be at least 250. n_iter_without_progress : int, default=300 Maximum number of iterations without progress before we abort the optimization, used after 250 initial iterations with early exaggeration. Note that progress is only checked every 50 iterations so this value is rounded to the next multiple of 50.  versionadded 0.17 parameter  n_iter_without_progress to control stopping criteria. min_grad_norm : float, default=1e-7 If the gradient norm is below this threshold, the optimization will be stopped. metric : str or callable, default='euclidean' The metric to use when calculating distance between instances in a feature array. If metric is a string, it must be one of the options allowed by scipy.spatial.distance.pdist for its metric parameter, or a metric listed in pairwise.PAIRWISE_DISTANCE_FUNCTIONS. If metric is \"precomputed\", X is assumed to be a distance matrix. Alternatively, if metric is a callable function, it is called on each pair of instances (rows) and the resulting value recorded. The callable should take two arrays from X as input and return a value indicating the distance between them. The default is \"euclidean\" which is interpreted as squared euclidean distance. init : {'random', 'pca'} or ndarray of shape (n_samples, n_components), default='random' Initialization of embedding. Possible options are 'random', 'pca', and a numpy array of shape (n_samples, n_components). PCA initialization cannot be used with precomputed distances and is usually more globally stable than random initialization. verbose : int, default=0 Verbosity level. random_state : int, RandomState instance or None, default=None Determines the random number generator. Pass an int for reproducible results across multiple function calls. Note that different initializations might result in different local minima of the cost function. See :term:  Glossary   . method : str, default='barnes_hut' By default the gradient calculation algorithm uses Barnes-Hut approximation running in O(NlogN) time. method='exact' will run on the slower, but exact, algorithm in O(N^2) time. The exact algorithm should be used when nearest-neighbor errors need to be better than 3%. However, the exact method cannot scale to millions of examples.  versionadded 0.17 Approximate optimization  method via the Barnes-Hut. angle : float, default=0.5 Only used if method='barnes_hut' This is the trade-off between speed and accuracy for Barnes-Hut T-SNE. 'angle' is the angular size (referred to as theta in [3]) of a distant node as measured from a point. If this size is below 'angle' then it is used as a summary node of all points contained within it. This method is not very sensitive to changes in this parameter in the range of 0.2 - 0.8. Angle less than 0.2 has quickly increasing computation time and angle greater 0.8 has quickly increasing error. n_jobs : int, default=None The number of parallel jobs to run for neighbors search. This parameter has no impact when  metric=\"precomputed\" or ( metric=\"euclidean\" and  method=\"exact\" ).  None means 1 unless in a :obj: joblib.parallel_backend context.  -1 means using all processors. See :term: Glossary   for more details.  versionadded 0.22 square_distances : True or 'legacy', default='legacy' Whether TSNE should square the distance values.  'legacy' means that distance values are squared only when  metric=\"euclidean\" .  True means that distance values are squared for all metrics.  versionadded 0.24 Added to provide backward compatibility during deprecation of legacy squaring behavior.  deprecated 0.24 Legacy squaring behavior was deprecated in 0.24. The  'legacy' value will be removed in 1.1 (renaming of 0.26), at which point the default value will change to  True . Attributes      embedding_ : array-like of shape (n_samples, n_components) Stores the embedding vectors. kl_divergence_ : float Kullback-Leibler divergence after optimization. n_iter_ : int Number of iterations run. Examples     >>> import numpy as np >>> from sklearn.manifold import TSNE >>> X = np.array( 0, 0, 0], [0, 1, 1], [1, 0, 1], [1, 1, 1 ) >>> X_embedded = TSNE(n_components=2).fit_transform(X) >>> X_embedded.shape (4, 2) References      [1] van der Maaten, L.J.P.; Hinton, G.E. Visualizing High-Dimensional Data Using t-SNE. Journal of Machine Learning Research 9:2579-2605, 2008. [2] van der Maaten, L.J.P. t-Distributed Stochastic Neighbor Embedding https: lvdmaaten.github.io/tsne/ [3] L.J.P. van der Maaten. Accelerating t-SNE using Tree-Based Algorithms. Journal of Machine Learning Research 15(Oct):3221-3245, 2014. https: lvdmaaten.github.io/publications/papers/JMLR_2014.pdf"
},
{
"ref":"partycls.dim_reduction.TSNE.symbol",
"url":18,
"doc":""
},
{
"ref":"partycls.dim_reduction.TSNE.full_name",
"url":18,
"doc":""
},
{
"ref":"partycls.dim_reduction.TSNE.reduce",
"url":18,
"doc":"Project the input features onto a reduced space using t-distributed stochastic neighbor embedding. Parameters      X : numpy.ndarray Features in the original space. Returns    - numpy.ndarray Features in the reduced space.",
"func":1
},
{
"ref":"partycls.dim_reduction.LocallyLinearEmbedding",
"url":18,
"doc":"Locally Linear Embedding Read more in the :ref: User Guide   . Parameters      n_neighbors : int, default=5 number of neighbors to consider for each point. n_components : int, default=2 number of coordinates for the manifold reg : float, default=1e-3 regularization constant, multiplies the trace of the local covariance matrix of the distances. eigen_solver : {'auto', 'arpack', 'dense'}, default='auto' auto : algorithm will attempt to choose the best method for input data arpack : use arnoldi iteration in shift-invert mode. For this method, M may be a dense matrix, sparse matrix, or general linear operator. Warning: ARPACK can be unstable for some problems. It is best to try several random seeds in order to check results. dense : use standard dense matrix operations for the eigenvalue decomposition. For this method, M must be an array or matrix type. This method should be avoided for large problems. tol : float, default=1e-6 Tolerance for 'arpack' method Not used if eigen_solver 'dense'. max_iter : int, default=100 maximum number of iterations for the arpack solver. Not used if eigen_solver 'dense'. method : {'standard', 'hessian', 'modified', 'ltsa'}, default='standard' standard : use the standard locally linear embedding algorithm. see reference [1] hessian : use the Hessian eigenmap method. This method requires  n_neighbors > n_components  (1 + (n_components + 1) / 2 see reference [2] modified : use the modified locally linear embedding algorithm. see reference [3] ltsa : use local tangent space alignment algorithm see reference [4] hessian_tol : float, default=1e-4 Tolerance for Hessian eigenmapping method. Only used if  method  'hessian' modified_tol : float, default=1e-12 Tolerance for modified LLE method. Only used if  method  'modified' neighbors_algorithm : {'auto', 'brute', 'kd_tree', 'ball_tree'}, default='auto' algorithm to use for nearest neighbors search, passed to neighbors.NearestNeighbors instance random_state : int, RandomState instance, default=None Determines the random number generator when  eigen_solver  'arpack'. Pass an int for reproducible results across multiple function calls. See :term:  Glossary   . n_jobs : int or None, default=None The number of parallel jobs to run.  None means 1 unless in a :obj: joblib.parallel_backend context.  -1 means using all processors. See :term: Glossary   for more details. Attributes      embedding_ : array-like, shape [n_samples, n_components] Stores the embedding vectors reconstruction_error_ : float Reconstruction error associated with  embedding_ nbrs_ : NearestNeighbors object Stores nearest neighbors instance, including BallTree or KDtree if applicable. Examples     >>> from sklearn.datasets import load_digits >>> from sklearn.manifold import LocallyLinearEmbedding >>> X, _ = load_digits(return_X_y=True) >>> X.shape (1797, 64) >>> embedding = LocallyLinearEmbedding(n_components=2) >>> X_transformed = embedding.fit_transform(X[:100]) >>> X_transformed.shape (100, 2) References       [1] Roweis, S. & Saul, L. Nonlinear dimensionality reduction by locally linear embedding. Science 290:2323 (2000).  [2] Donoho, D. & Grimes, C. Hessian eigenmaps: Locally linear embedding techniques for high-dimensional data. Proc Natl Acad Sci U S A. 100:5591 (2003).  [3] Zhang, Z. & Wang, J. MLLE: Modified Locally Linear Embedding Using Multiple Weights. http: citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.70.382  [4] Zhang, Z. & Zha, H. Principal manifolds and nonlinear dimensionality reduction via tangent space alignment. Journal of Shanghai Univ. 8:406 (2004)"
},
{
"ref":"partycls.dim_reduction.LocallyLinearEmbedding.symbol",
"url":18,
"doc":""
},
{
"ref":"partycls.dim_reduction.LocallyLinearEmbedding.full_name",
"url":18,
"doc":""
},
{
"ref":"partycls.dim_reduction.LocallyLinearEmbedding.reduce",
"url":18,
"doc":"Project the input features onto a reduced space using locally linear embedding. Parameters      X : numpy.ndarray Features in the original space. Returns    - numpy.ndarray Features in the reduced space.",
"func":1
},
{
"ref":"partycls.dim_reduction.AutoEncoder",
"url":18,
"doc":"Multi-layer Perceptron regressor. This model optimizes the squared-loss using LBFGS or stochastic gradient descent.  versionadded 0.18 Parameters      hidden_layer_sizes : tuple, length = n_layers - 2, default=(100,) The ith element represents the number of neurons in the ith hidden layer. activation : {'identity', 'logistic', 'tanh', 'relu'}, default='relu' Activation function for the hidden layer. - 'identity', no-op activation, useful to implement linear bottleneck, returns f(x) = x - 'logistic', the logistic sigmoid function, returns f(x) = 1 / (1 + exp(-x . - 'tanh', the hyperbolic tan function, returns f(x) = tanh(x). - 'relu', the rectified linear unit function, returns f(x) = max(0, x) solver : {'lbfgs', 'sgd', 'adam'}, default='adam' The solver for weight optimization. - 'lbfgs' is an optimizer in the family of quasi-Newton methods. - 'sgd' refers to stochastic gradient descent. - 'adam' refers to a stochastic gradient-based optimizer proposed by Kingma, Diederik, and Jimmy Ba Note: The default solver 'adam' works pretty well on relatively large datasets (with thousands of training samples or more) in terms of both training time and validation score. For small datasets, however, 'lbfgs' can converge faster and perform better. alpha : float, default=0.0001 L2 penalty (regularization term) parameter. batch_size : int, default='auto' Size of minibatches for stochastic optimizers. If the solver is 'lbfgs', the classifier will not use minibatch. When set to \"auto\",  batch_size=min(200, n_samples) learning_rate : {'constant', 'invscaling', 'adaptive'}, default='constant' Learning rate schedule for weight updates. - 'constant' is a constant learning rate given by 'learning_rate_init'. - 'invscaling' gradually decreases the learning rate  learning_rate_ at each time step 't' using an inverse scaling exponent of 'power_t'. effective_learning_rate = learning_rate_init / pow(t, power_t) - 'adaptive' keeps the learning rate constant to 'learning_rate_init' as long as training loss keeps decreasing. Each time two consecutive epochs fail to decrease training loss by at least tol, or fail to increase validation score by at least tol if 'early_stopping' is on, the current learning rate is divided by 5. Only used when solver='sgd'. learning_rate_init : double, default=0.001 The initial learning rate used. It controls the step-size in updating the weights. Only used when solver='sgd' or 'adam'. power_t : double, default=0.5 The exponent for inverse scaling learning rate. It is used in updating effective learning rate when the learning_rate is set to 'invscaling'. Only used when solver='sgd'. max_iter : int, default=200 Maximum number of iterations. The solver iterates until convergence (determined by 'tol') or this number of iterations. For stochastic solvers ('sgd', 'adam'), note that this determines the number of epochs (how many times each data point will be used), not the number of gradient steps. shuffle : bool, default=True Whether to shuffle samples in each iteration. Only used when solver='sgd' or 'adam'. random_state : int, RandomState instance, default=None Determines random number generation for weights and bias initialization, train-test split if early stopping is used, and batch sampling when solver='sgd' or 'adam'. Pass an int for reproducible results across multiple function calls. See :term: Glossary   . tol : float, default=1e-4 Tolerance for the optimization. When the loss or score is not improving by at least  tol for  n_iter_no_change consecutive iterations, unless  learning_rate is set to 'adaptive', convergence is considered to be reached and training stops. verbose : bool, default=False Whether to print progress messages to stdout. warm_start : bool, default=False When set to True, reuse the solution of the previous call to fit as initialization, otherwise, just erase the previous solution. See :term: the Glossary   . momentum : float, default=0.9 Momentum for gradient descent update. Should be between 0 and 1. Only used when solver='sgd'. nesterovs_momentum : bool, default=True Whether to use Nesterov's momentum. Only used when solver='sgd' and momentum > 0. early_stopping : bool, default=False Whether to use early stopping to terminate training when validation score is not improving. If set to true, it will automatically set aside 10% of training data as validation and terminate training when validation score is not improving by at least  tol for  n_iter_no_change consecutive epochs. Only effective when solver='sgd' or 'adam' validation_fraction : float, default=0.1 The proportion of training data to set aside as validation set for early stopping. Must be between 0 and 1. Only used if early_stopping is True beta_1 : float, default=0.9 Exponential decay rate for estimates of first moment vector in adam, should be in [0, 1). Only used when solver='adam' beta_2 : float, default=0.999 Exponential decay rate for estimates of second moment vector in adam, should be in [0, 1). Only used when solver='adam' epsilon : float, default=1e-8 Value for numerical stability in adam. Only used when solver='adam' n_iter_no_change : int, default=10 Maximum number of epochs to not meet  tol improvement. Only effective when solver='sgd' or 'adam'  versionadded 0.20 max_fun : int, default=15000 Only used when solver='lbfgs'. Maximum number of function calls. The solver iterates until convergence (determined by 'tol'), number of iterations reaches max_iter, or this number of function calls. Note that number of function calls will be greater than or equal to the number of iterations for the MLPRegressor.  versionadded 0.22 Attributes      loss_ : float The current loss computed with the loss function. best_loss_ : float The minimum loss reached by the solver throughout fitting. loss_curve_ : list of shape ( n_iter_ ,) The ith element in the list represents the loss at the ith iteration. t_ : int The number of training samples seen by the solver during fitting. coefs_ : list of shape (n_layers - 1,) The ith element in the list represents the weight matrix corresponding to layer i. intercepts_ : list of shape (n_layers - 1,) The ith element in the list represents the bias vector corresponding to layer i + 1. n_iter_ : int The number of iterations the solver has ran. n_layers_ : int Number of layers. n_outputs_ : int Number of outputs. out_activation_ : str Name of the output activation function. loss_curve_ : list of shape (n_iters,) Loss value evaluated at the end of each training step. t_ : int Mathematically equals  n_iters  X.shape[0] , it means  time_step and it is used by optimizer's learning rate scheduler. Examples     >>> from sklearn.neural_network import MLPRegressor >>> from sklearn.datasets import make_regression >>> from sklearn.model_selection import train_test_split >>> X, y = make_regression(n_samples=200, random_state=1) >>> X_train, X_test, y_train, y_test = train_test_split(X, y,  . random_state=1) >>> regr = MLPRegressor(random_state=1, max_iter=500).fit(X_train, y_train) >>> regr.predict(X_test[:2]) array([-0.9 ., -7.1 .]) >>> regr.score(X_test, y_test) 0.4 . Notes   - MLPRegressor trains iteratively since at each time step the partial derivatives of the loss function with respect to the model parameters are computed to update the parameters. It can also have a regularization term added to the loss function that shrinks model parameters to prevent overfitting. This implementation works with data represented as dense and sparse numpy arrays of floating point values. References      Hinton, Geoffrey E. \"Connectionist learning procedures.\" Artificial intelligence 40.1 (1989): 185-234. Glorot, Xavier, and Yoshua Bengio. \"Understanding the difficulty of training deep feedforward neural networks.\" International Conference on Artificial Intelligence and Statistics. 2010. He, Kaiming, et al. \"Delving deep into rectifiers: Surpassing human-level performance on imagenet classification.\" arXiv preprint arXiv:1502.01852 (2015). Kingma, Diederik, and Jimmy Ba. \"Adam: A method for stochastic optimization.\" arXiv preprint arXiv:1412.6980 (2014)."
},
{
"ref":"partycls.dim_reduction.AutoEncoder.symbol",
"url":18,
"doc":""
},
{
"ref":"partycls.dim_reduction.AutoEncoder.full_name",
"url":18,
"doc":""
},
{
"ref":"partycls.dim_reduction.AutoEncoder.n_components",
"url":18,
"doc":""
},
{
"ref":"partycls.dim_reduction.AutoEncoder.reduce",
"url":18,
"doc":"Project the input features onto a reduced space using a neural network autoencoder. The dimension of the reduced space is the number of nodes in the bottleneck layer. Parameters      X : numpy.ndarray Features in the original space. Returns    - numpy.ndarray Features in the reduced space.",
"func":1
}
]