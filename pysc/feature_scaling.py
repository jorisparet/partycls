from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler

__all__ = ['ZScore', 'MinMax']

class ZScore(StandardScaler):
    
    symbol = 'zscore'
    full_name = 'Z-Score'
    
    def scale(self, X):
        return self.fit_transform(X)
    
class MinMax(MinMaxScaler):
    
    symbol = 'minmax'
    full_name = 'Min-Max'
    
    def scale(self, X):
        return self.fit_transform(X)