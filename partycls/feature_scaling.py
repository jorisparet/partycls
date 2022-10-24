"""
Feature scaling techniques, to be performed on a dataset stored in a numpy
array.
"""

from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import MaxAbsScaler
from sklearn.preprocessing import RobustScaler

__all__ = ['ZScore', 'MinMax', 'MaxAbs', 'Robust']


class ZScore(StandardScaler):

    symbol = 'zscore'
    full_name = 'Z-Score'

    def scale(self, X):
        """
        Standardize features by removing the mean and scaling to unit variance.

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
        Transform features by scaling each feature to a given range (default 
        is :math:`[0,1]`).

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


class MaxAbs(MaxAbsScaler):

    symbol = 'maxabs'
    full_name = 'Max-Abs'

    def scale(self, X):
        """
        Scale each feature by its maximum absolute value.

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


class Robust(RobustScaler):

    symbol = 'robust'
    full_name = 'Robust'

    def scale(self, X):
        """
        Scale features using statistics that are robust to outliers.

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
