from sklearn.cluster import KMeans as _KMeans
from sklearn.mixture import GaussianMixture as _GaussianMixture
from pysc.descriptor import StructuralDescriptor, DummyDescriptor
import numpy

class Clustering:
    
    def __init__(self, n_clusters=2, method='kmeans', n_init=1):
        self.n_clusters = n_clusters
        self.method = method
        self.n_init = n_init
        self.labels = None
    
    def __repr__(self):
        # `self.symbol` is defined in child classes
       return self.symbol
    
    def fit(X):
        pass

    @property
    def fractions(self):
        if self.labels is not None:
            f_k = numpy.empty(self.n_clusters, dtype=numpy.float64)
            for k in range(self.n_clusters):
                f_k[k] = numpy.sum(self.labels == k)
            return f_k / len(self.labels)

    @property
    def populations(self):
        if self.labels is not None:
            n_k = numpy.empty(self.n_clusters, dtype=numpy.int64)
            for k in range(self.n_clusters):
                n_k[k] = numpy.sum(self.labels == k)
            return n_k
    
    def centroids(self, X):
        """
        Returns the clusters' centroids of dataset `X` using the labels.
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
    
class KMeans(Clustering):
    
    symbol = 'kmeans'
    
    def __init__(self, n_clusters=2, n_init=1):
        self.symbol = 'kmeans'
        self.full_name = 'K-Means'
        Clustering.__init__(self, n_clusters=n_clusters, method=self.symbol, n_init=n_init)
        
    def fit(self, X):
        """
        Run the K-Means algorithm on `X`.
        The predicted labels are updated in the attribute `labels` of 
        the current instance of `KMeans`.
        """
        kmeans = _KMeans(n_clusters=self.n_clusters,
                         n_init=self.n_init)
        if isinstance(X, StructuralDescriptor):
            kmeans.fit(X.features)
        else:
            kmeans.fit(X)
        self.labels = kmeans.labels_
        
class GaussianMixture(Clustering):
    
    def __init__(self, n_clusters=2, n_init=1):
        self.symbol = 'gmm'
        self.full_name = 'Gaussian Mixture'
        Clustering.__init__(self, n_clusters=n_clusters, method=self.symbol, n_init=n_init)
        
    def fit(self, X):
        """
        Run the EM algorithm on `X` using a mixture of Gaussians.
        The predicted labels are updated in the attribute `labels` of the 
        current instance of `GaussianMixture`.
        """
        gmm = _GaussianMixture(n_components=self.n_clusters,
                               n_init=self.n_init)
        if isinstance(X, StructuralDescriptor):
            gmm.fit(X.features)
        else:
            gmm.fit(X)
        self.labels = gmm.predict(X)
        
class CommunityInference(Clustering):
    
    def __init__(self, n_clusters=2, n_init=1):
        self.symbol = 'cinf'
        self.full_name = 'Community Inference'
        Clustering.__init__(self, n_clusters=n_clusters, method=self.symbol, n_init=n_init)
        self.mutual_information = None
        
    def fit(self, X):
        """
        Community inference algorithm.
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
        for n in range(self.n_init-1):
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
        N = descriptor.size
        K = self.n_clusters
        Km1 = K - 1 # loop invariant
        
        # randomly set the labels
        self.labels = numpy.array([random.randint(0,Km1) for n in range(N)])
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
        dx = descriptor.grid[1] - descriptor.grid[0]
        def _mutual_information(P_average, P_k, f_k, dx):
            MI = 0.0
            f_x = numpy.empty_like(P_average)
            for k in range(len(f_k)):
                for m in range(len(P_average)):
                    if P_k[k,m] != 0.0 and P_average[m] != 0.0:
                        f_x[m] = P_k[k,m] * numpy.log2( P_k[k,m] / P_average[m] )
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