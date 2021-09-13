class FeatureSelector:
    
    def __init__(self):
        self.selections = []
    
    def show_feature_space(self, X_1, X_2, method='heatmap', show_selections=True, **kwargs):
        import matplotlib.pyplot as plt
        
        if method == 'heatmap':
            plt.hist2d(X_1, X_2, **kwargs)
        elif method == 'scatter':
            plt.scatter(X_1, X_2, **kwargs)
        else:
            print('Unknown method.')
            
        if show_selections:
            for S in self.selections:
                from matplotlib.patches import Rectangle as RectanglePatch
                if S.shape == 'rectangle':
                    shape = RectanglePatch((S.x, S.y), S.width, S.height, **S.kwargs)
                else:
                    print('Unknown shape.')
                plt.gca().add_patch(shape)
        plt.show()
        
    def add_rectangular_selection(self, x, y, width, height, **kwargs):
        R = Rectangle(x, y, width, height, **kwargs)
        self.selections.append(R)

    def clear_selections(self):
        self.selections = []
        
    def selected(self, X_1, X_2):
        selected_list = [False for i in range(X_1.size)]
        for n, (x, y) in enumerate(zip(X_1,X_2)):
            for m, selection in enumerate(self.selections):
                selected_list[n] = selection.is_inside(x, y)
        return selected_list
        
class Rectangle:
    
    def __init__(self, x, y, width, height, label=None, **kwargs):
        self.x = x
        self.y = y
        self.width = width
        self.height = height
        self.label = label
        self.kwargs = kwargs
        self.shape = 'rectangle'
        
    def is_inside(self, x, y):
        return self.x < x < self.x+self.width and self.y < y < self.y+self.height
        
    def __repr__(self):
        return 'Rectangle=(origin={}, width={}, height={}, label={})'.format((self.x,self.y),
                                                                             self.width,
                                                                             self.height,
                                                                             self.label)