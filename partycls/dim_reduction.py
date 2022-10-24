"""
Dimensionality reduction techniques (linear and non-linear), to be performed
on a dataset stored in a numpy array.
"""

import numpy
from sklearn.decomposition import PCA as _PCA
from sklearn.manifold import TSNE as _TSNE
from sklearn.manifold import LocallyLinearEmbedding as _LocallyLinearEmbedding
from sklearn.neural_network import MLPRegressor

__all__ = ['PCA', 'TSNE', 'LocallyLinearEmbedding', 'AutoEncoder']


class PCA(_PCA):

    symbol = 'pca'
    full_name = 'Principal Component Analysis (PCA)'

    def reduce(self, X):
        """
        Project the input features onto a reduced space using principal
        component analysis.

        Parameters
        ----------
        X : numpy.ndarray
            Features in the original space.

        Returns
        -------
        numpy.ndarray
            Features in the reduced space.
        """
        return self.fit_transform(X)


class TSNE(_TSNE):

    symbol = 'tsne'
    full_name = 't-distributed Stochastic Neighbor Embedding (t-SNE)'

    def reduce(self, X):
        """
        Project the input features onto a reduced space using t-distributed 
        stochastic neighbor embedding.

        Parameters
        ----------
        X : numpy.ndarray
            Features in the original space.

        Returns
        -------
        numpy.ndarray
            Features in the reduced space.
        """
        return self.fit_transform(X)


class LocallyLinearEmbedding(_LocallyLinearEmbedding):

    symbol = 'lle'
    full_name = 'Locally Linear Embedding (LLE)'

    def reduce(self, X):
        """
        Project the input features onto a reduced space using locally
        linear embedding.

        Parameters
        ----------
        X : numpy.ndarray
            Features in the original space.

        Returns
        -------
        numpy.ndarray
            Features in the reduced space.
        """
        return self.fit_transform(X)


class AutoEncoder(MLPRegressor):

    symbol = 'ae'
    full_name = 'Neural-Network Auto-Encoder (AE)'

    def __init__(self, layers=(100, 2, 100), activation='relu', solver='adam', alpha=1e-4):
        MLPRegressor.__init__(self, hidden_layer_sizes=layers,
                              activation=activation, solver=solver,
                              alpha=alpha)

    @property
    def n_components(self):
        """
        Number of nodes at the level of the bottleneck layer (*i.e.* dimension
        after reduction).
        """
        return min(self.hidden_layer_sizes)

    def reduce(self, X):
        """
        Project the input features onto a reduced space using a neural network
        autoencoder. The dimension of the reduced space is the number of 
        nodes in the bottleneck layer.

        Parameters
        ----------
        X : numpy.ndarray
            Features in the original space.

        Returns
        -------
        numpy.ndarray
            Features in the reduced space.
        """

        # Train the network to reproduce its input as output
        self.fit(X, X)

        # Mean absolute error
        Y_pred = self.predict(X)
        # MAE
        MAE = numpy.abs(Y_pred - X).mean()
        self.mean_absolute_error = MAE
        # MSE / MSD
        MSE = 0.0
        MSD = 0.0
        Xmean = numpy.mean(X, axis=0)
        for i in range(X.shape[0]):
            MSE += numpy.sum((X[i] - self.predict(X[i].reshape(1, -1)))**2)
            MSD += numpy.sum((X[i] - Xmean)**2)
        MSE /= X.shape[0]
        MSD /= X.shape[0]
        self.mean_squared_error = MSE
        self.mean_squared_deviation = MSD

        # Weights and biases
        W = self.coefs_
        biases = self.intercepts_

        # Keep the encoder part only
        bottleneck_index = self.hidden_layer_sizes.index(self.n_components)
        encoder_weights = W[0:bottleneck_index + 1]
        encoder_biases = biases[0:bottleneck_index + 1]

        # Encode data
        X_red = X
        for index, (w, b) in enumerate(zip(encoder_weights, encoder_biases)):
            if index + 1 == len(encoder_weights):
                X_red = X_red @ w + b
            else:
                # Use the right activation function here
                if self.activation == 'relu':
                    X_red = numpy.maximum(0, X_red @ w + b)
                if self.activation == 'tanh':
                    X_red = numpy.tanh(X_red @ w + b)
                if self.activation == 'identity':
                    X_red = X_red @ w + b
                if self.activation == 'logistic':
                    X_red = 1.0 / (1.0 + numpy.exp(-(X_red @ w + b)))

        # Return the dataset in low dimension
        return X_red
