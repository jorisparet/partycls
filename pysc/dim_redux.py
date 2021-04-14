from sklearn.decomposition import PCA as _PCA
from sklearn.manifold import TSNE as _TSNE
from sklearn.manifold import LocallyLinearEmbedding as _LocallyLinearEmbedding
from sklearn.neural_network import MLPRegressor
import numpy

class PCA(_PCA):
    
    symbol = 'pca'
    full_name = 'Principal Component Analysis (PCA)'
    
    def reduce(self, X):
        return self.fit_transform(X)

class TSNE(_TSNE):
    
    symbol = 'tsne'
    full_name = 't-distributed Stochastic Neighbor Embedding (t-SNE)'

    def reduce(self, X):
        return self.fit_transform(X)
    
class LocallyLinearEmbedding(_LocallyLinearEmbedding):
    
    symbol = 'lle'
    full_name = 'Locally Linear Embedding (LLE)'
    
    def reduce(self, X):
        return self.fit_transform(X)

class AutoEncoder(MLPRegressor):
    
    symbol = 'ae'
    full_name = 'Neural-Network Auto-Encoder (AE)'
    
    def __init__(self, layers=(2,), activation='relu'):
        MLPRegressor.__init__(hidden_layer_sizes=layers,
                              activation=activation)
    
    def reduce(self, X):
        """
        Neural-network-based auto-encoder.
        """
        
        # Train the network to reproduce its input as output
        self.fit(X, X)
        
        # Mean absolute error
        Y_pred = self.predict(X)
        err = numpy.abs(Y_pred - X).mean()
        self.mean_abs_error = err
        
        # Weights and biases
        W = self.coefs_
        biases = self.intercepts_
        
        # Keep the encoder part only
        n_components = min(self.hidden_layer_sizes)
        bottleneck_index = self.hidden_layer_sizes.index(n_components)
        encoder_weights = W[0:bottleneck_index+1]
        encoder_biases = biases[0:bottleneck_index+1]
        
        # Encode data
        X_red = X
        for index, (w, b) in enumerate(zip(encoder_weights, encoder_biases)):
            if index+1 == len(encoder_weights):
                X_red = X_red @ w + b 
            else:
                # Use the right activation function here
                if self.activation == 'relu':
                    X_red = numpy.maximum(0, X_red @ w + b)
                if self.activation == 'tanh':
                    X_red = numpy.tanh(X_red @ w + b)
                if self.activation == 'identity':
                    raise NotImplementedError
                if self.activation == 'logistic':
                    raise NotImplementedError
                    
        # Return the dataset in low dimension
        return X_red