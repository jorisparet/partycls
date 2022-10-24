"""
Clustering algorithms.
"""

import numpy
from sklearn.cluster import KMeans as _KMeans
from sklearn.mixture import GaussianMixture as _GaussianMixture
from partycls.descriptors import StructuralDescriptor, DummyDescriptor, BondOrientationalDescriptor

__all__ = ['Clustering', 'KMeans', 'GaussianMixture', 'CommunityInference']


class Clustering:
    """
    Base class for clustering methods.

    Attributes
    ----------
    n_clusters : int
        Number of clusters.
    
    n_init : int
        Number of times the clustering is run.
    
    labels : list
        Cluster labels. The default is ``None``. Initialized after the ``fit``
        method is called.
    """

    def __init__(self, n_clusters=2, n_init=1, backend=None):
        """
        If a scikit-learn compatible ``backend`` is available, it will be 
        used within Strategy.

        Parameters
        ----------
        n_clusters : int, default: 2
            Requested number of clusters.
        
        n_init : int, default: 1
            Number of times the clustering will be run with different seeds.
                
        backend : scikit-learn compatible backend, default: None
            Backend used for the clustering method. If provided, it must
            be an object implementing an scikit-learn compatible interface,
            with a ``fit`` method and a ``labels_`` attribute. Duck typing
            is assumed.
        """
        self.n_clusters = n_clusters
        self.n_init = n_init
        self.backend = backend
        self.labels = None

    def fit(self, X):
        """
        Run a scikit-learn compatible clustering backend (if available) on ``X``.

        Subclasses implementing a specific clustering algorithm must
        override this method.

        Parameters
        ----------
        X : numpy.ndarray
            Dataset matrix for which to compute the clusters.

        Returns
        -------
        None
        """
        if self.backend is not None:
            if hasattr(X, 'features'):
                self.backend.fit(X.features)
            else:
                self.backend.fit(X)
            self.labels = self.backend.labels_

    @property
    def fractions(self):
        """
        ``numpy.ndarray`` with the fractions of particles in each cluster.
        """
        if self.labels is not None:
            f_k = numpy.empty(self.n_clusters, dtype=numpy.float64)
            for k in range(self.n_clusters):
                f_k[k] = numpy.sum(self.labels == k)
            return f_k / len(self.labels)

    @property
    def populations(self):
        """
        ``numpy.ndarray`` with the number of particles in each cluster.
        """
        if self.labels is not None:
            n_k = numpy.empty(self.n_clusters, dtype=numpy.int64)
            for k in range(self.n_clusters):
                n_k[k] = numpy.sum(self.labels == k)
            return n_k

    def centroids(self, X):
        """
        Central feature vector of each cluster.
        
        Each object in the dataset over which the clustering was performed is 
        assigned a discrete label. This label represents the index of the 
        nearest cluster center to which this object belongs. The centroid (*i.e.* 
        the cluster center), is thus the average feature vector of all the 
        objects in the cluster.
        
        Cluster memberships of the objects are stored in the ``labels``
        attribute. Coordinates of the centroids can then be calculated for an
        arbitrary dataset ``X``, provided it has the same shape as the original 
        dataset used for the clustering.

        Parameters
        ----------
        X : numpy.ndarray
            Array of features (dataset) for which to compute the centroids.

        Returns
        -------
        C_k : numpy.ndarray
            Cluster centroids. ``C_k[n]`` is the coordinates of the n-th cluster 
            center.
        """
        n_features = X.shape[1]
        C_k = numpy.zeros((self.n_clusters, n_features), dtype=numpy.float64)
        n_k = self.populations
        for k in range(self.n_clusters):
            for n in range(len(self.labels)):
                if self.labels[n] == k:
                    C_k[k] += X[n]
            C_k[k] = C_k[k] / n_k[k]
        return C_k

    def __str__(self):
        rep = 'Clustering(method="{}", n_clusters={}, n_init={})'
        return rep.format(self.full_name, self.n_clusters, self.n_init)

    def __repr__(self):
        return self.__str__()

class KMeans(Clustering):
    """
    KMeans clustering.
    
    This class relies on the class ``KMeans`` from the machine learning package 
    scikit-learn. An instance of ``sklearn.cluster.KMeans`` is created when 
    calling the ``fit`` method, and is then accessible through the ``backend``
    attribute for later use. See scikit-learn's documentation for more information on
    the original class.
    """

    def __init__(self, n_clusters=2, n_init=1):
        """
        Parameters
        ----------
        n_clusters : int, default: 2
            Requested number of clusters.
        
        n_init : int, default: 1
            Number of times the clustering will be run with different seeds.
        """
        self.symbol = 'kmeans'
        self.full_name = 'K-Means'
        Clustering.__init__(self, n_clusters=n_clusters, n_init=n_init)

    def fit(self, X):
        """
        Run the K-Means algorithm on ``X``.
        The predicted labels are updated in the attribute ``labels`` of 
        the current instance of ``KMeans``.

        Parameters
        ----------
        X : numpy.ndarray
            Dataset matrix for which to compute the clusters.

        Returns
        -------
        None
        """
        self.backend = _KMeans(n_clusters=self.n_clusters,
                               n_init=self.n_init)
        if hasattr(X, 'features'):
            self.backend.fit(X.features)
        else:
            self.backend.fit(X)
        self.labels = self.backend.labels_


class GaussianMixture(Clustering):
    """
    Gaussian Mixture.
    
    This class relies on the class ``GaussianMixture`` from the machine learning 
    package scikit-learn. An instance of ``sklearn.mixture.GaussianMixture`` is 
    created when calling the ``fit`` method, and is then accessible through the 
    ``backend`` attribute for later use. See scikit-learn's documentation for more 
    information on the original class.
    """

    def __init__(self, n_clusters=2, n_init=1):
        """
        Parameters
        ----------
        n_clusters : int, default: 2
            Requested number of clusters.
        
        n_init : int, default: 1
            Number of times the clustering will be run with different seeds.       
        """
        self.symbol = 'gmm'
        self.full_name = 'Gaussian Mixture'
        Clustering.__init__(self, n_clusters=n_clusters, n_init=n_init)

    def fit(self, X):
        """
        Run the expectation-maximization algorithm on ``X`` using a mixture of 
        Gaussians. The predicted labels are updated in the attribute ``labels`` of the 
        current instance of ``GaussianMixture``.

        Parameters
        ----------
        X : numpy.ndarray
            Dataset matrix for which to compute the clusters.

        Returns
        -------
        None
        """
        self.backend = _GaussianMixture(n_components=self.n_clusters,
                                        n_init=self.n_init)
        if hasattr(X, 'features'):
            self.backend.fit(X.features)
        else:
            self.backend.fit(X)
        self.labels = self.backend.predict(X)


class CommunityInference(Clustering):
    """
    Community Inference is a hard clustering method based on information 
    theory. See https://doi.org/10.1063/5.0004732 (Paret et. al) for more 
    details.
    """

    def __init__(self, n_clusters=2, n_init=1):
        """
        Parameters
        ----------
        n_clusters : int, default: 2
            Requested number of clusters.
        
        n_init : int, default: 1
            Number of times the clustering will be run with different seeds.           
        """
        self.symbol = 'cinf'
        self.full_name = 'Community Inference'
        Clustering.__init__(self, n_clusters=n_clusters, n_init=n_init)
        self.mutual_information = None

    def fit(self, X):
        """
        Run the community inference algorithm on ``X``, where ``X``
        is an instance of ``StructuralDescriptor`` with a ``normalize``
        method. Otherwise ``X`` is converted to a dummy descriptor.

        Parameters
        ----------
        X : StructuralDescriptor
            Descriptor on which the community algorithm inference will be run.

        Returns
        -------
        None
        """
        if isinstance(X, StructuralDescriptor):
            descriptor = X
        else:
            descriptor = DummyDescriptor()
            features = numpy.empty_like(X)
            bins = X.shape[1]
            for n in range(features.shape[0]):
                features[n] = numpy.histogram(X[n], bins=bins)[0]
            descriptor.features = features

        MI_previous, labels_previous = self._inference_loop(descriptor)
        for n in range(self.n_init - 1):
            MI_next, labels_next = self._inference_loop(descriptor)
            # optimization `n` is worse than the previous one
            if MI_next < MI_previous:
                self.mutual_information = MI_previous
                self.labels = labels_previous
            # optimization `n` is better than the previous one
            # it becomes the new standard
            else:
                MI_previous = MI_next
                labels_previous = labels_next

    def _inference_loop(self, descriptor):

        import random

        # shortcuts
        N = descriptor.n_samples
        K = self.n_clusters
        Km1 = K - 1  # loop invariant

        # randomly set the labels
        self.labels = numpy.array([random.randint(0, Km1) for n in range(N)])
        # populations and fractions
        n_k = self.populations
        f_k = self.fractions
        # community histograms
        dtype = descriptor.features.dtype
        H_k = numpy.zeros((K, descriptor.n_features), dtype=dtype)
        for i in range(N):
            k_i = self.labels[i]
            H_k[k_i] += descriptor.features[i]
        # community distributions
        P_k = numpy.empty_like(H_k, dtype=numpy.float64)
        for k in range(K):
            P_k[k] = descriptor.normalize(H_k[k] / n_k[k])
        # average distribution
        P_average = f_k @ P_k

        # community information
        if isinstance(descriptor, BondOrientationalDescriptor):
            # in case of non-successive values of l in the grid
            dx = 1.0
        else:
            dx = descriptor.grid[1] - descriptor.grid[0]

        def _mutual_information(P_average, P_k, f_k, dx):
            MI = 0.0
            f_x = numpy.empty_like(P_average)
            for k in range(len(f_k)):
                for m in range(len(P_average)):
                    if P_k[k, m] != 0.0 and P_average[m] != 0.0:
                        f_x[m] = P_k[k, m] * numpy.log2(P_k[k, m] / P_average[m])
                    else:
                        f_x[m] = 0.0

#                #TODO: is it more efficient?
#                f_x = P_k[k] * (numpy.log2(P_k[k],
#                         out=numpy.zeros_like(P_average),
#                         where=(P_k[k] != 0.0))
#                         - numpy.log2(P_average,
#                         out=numpy.zeros_like(P_average),
#                         where=(P_average != 0.0)))

                f_x = f_x * f_k[k]
                MI = MI + numpy.sum(f_x) * dx
            return MI
        self.mutual_information = _mutual_information(P_average, P_k, f_k, dx)

        # Begin loop
        no_change = 0
        while no_change < N:
            # sample particles
            sample = random.sample(range(N), N)
            for i in sample:
                k_i = self.labels[i]
                k_others = [k for k in range(K) if k != k_i]
                # lists to store attemps for relocation of `i`
                labels_list = [numpy.empty_like(self.labels) for k in range(Km1)]
                n_k_list = [numpy.empty_like(n_k) for k in range(Km1)]
                H_k_list = [numpy.empty_like(H_k) for k in range(Km1)]
                P_k_list = [numpy.empty_like(P_k) for k in range(Km1)]
                I_list = [0.0 for k in range(Km1)]
                # loop over all other available communnities
                for kn, k_j in enumerate(k_others):
                    labels_new = self.labels.copy()
                    n_k_new = n_k.copy()
                    H_k_new = H_k.copy()
                    P_k_new = P_k.copy()
                    # k_i --> k_j
                    labels_new[i] = k_j
                    labels_list[kn] = labels_new
                    # changes in populations and fractions
                    n_k_new[k_i] -= 1
                    n_k_new[k_j] += 1
                    n_k_list[kn] = n_k_new
                    f_k_new = n_k_new / N
                    # changes in histograms
                    H_k_new[k_i] -= descriptor.features[i]
                    H_k_new[k_j] += descriptor.features[i]
                    H_k_list[kn] = H_k_new
                    # changes in distributions
                    P_k_new[k_i] = descriptor.normalize(H_k_new[k_i] / n_k_new[k_i])
                    P_k_new[k_j] = descriptor.normalize(H_k_new[k_j] / n_k_new[k_j])
                    P_k_list[kn] = P_k_new
                    # change in community information
                    I_new = _mutual_information(P_average, P_k_new, f_k_new, dx)
                    I_list[kn] = I_new
                # find the most appropriate community for `i`
                best = numpy.argmax(I_list)
                if I_list[best] > self.mutual_information:
                    self.labels = labels_list[best]
                    self.mutual_information = I_list[best]
                    n_k = n_k_list[best]
                    H_k = H_k_list[best]
                    P_k = P_k_list[best]
                    no_change = 0
                else:
                    no_change += 1
        # final mutual information
        return self.mutual_information, self.labels.copy()
