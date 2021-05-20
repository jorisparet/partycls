from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler

__all__ = ['ZScore', 'MinMax']

class ZScore(StandardScaler):
    
    symbol = 'zscore'
    full_name = 'Z-Score'
    
    def scale(self, X):
        """
        Rescale the input features using Z-scores (or "standard scores").

        Parameters
        ----------
        X : numpy.ndarray
            Original features.

        Returns
        -------
        numpy.ndarray
            Scaled features.

        """
        return self.fit_transform(X)
    
class MinMax(MinMaxScaler):
    
    symbol = 'minmax'
    full_name = 'Min-Max'
    
    def scale(self, X):
        """
        Rescale the input features using min-max renormalization.

        Parameters
        ----------
        X : numpy.ndarray
            Original features.

        Returns
        -------
        numpy.ndarray
            Scaled features.

        """
        return self.fit_transform(X)