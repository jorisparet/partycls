class FeatureSelector:
    
    def __init__(self, descriptor, X_1, X_2):
        self.descriptor = descriptor
        self.X_1 =  X_1
        self.X_2 = X_2
        self.active_selections = []
        
    def add_rectangular_selection(self, x0, y0, x1, y1, label=None, **kwargs):
        '''
        Rectangular selection.
        
        +----------(x1,y1)
        |                |
        |                |
        |                |
        (x0,y0)----------+
        
        '''
        # Default name for label
        if label is None:
            label = 'S_{}'.format(len(self.selections))
            
        # Set particle property trajectory-wide to 0
        if not hasattr(self.descriptor.trajectory[0].particle[0], label):
            self.descriptor.trajectory.set_property(label, 0)
        # Set property to 1 if particle is in the current selection
        R = Rectangle(x0, y0, x1, y1, label=label, **kwargs)
        count = 0
        for frame in self.descriptor.groups[0]:
            for particle in frame:
                is_inside = int(R.is_inside(self.X_1[count], self.X_2[count]))
                setattr(particle, label, is_inside)
                count += 1
            
        # Add selection
        self.active_selections.append(R)

    def clear_selection(self, label):
        '''
        Clear all selections labeled as `label` and clear the associated
        particle property from the trajectory.
        '''
        # Delete particle property trajectory-wide           
        for system in self.descriptor.trajectory:
            for particle in system.particle:
                delattr(particle, label)
        # Delete selection(s) from list of active selections
        self.active_selections = [s for s in self.active_selections if s.label != label]

    def clear_all_selections(self):
        '''
        Clear all selections their associated
        particle properties from the trajectory.
        '''
        for selection in self.active_selections:
            self.clear_selection(selection.label)

    def show_feature_space(self, method='heatmap', show=True, show_selections=True, **kwargs):
        '''
        Show the feature space (X_1,X_2) in the form of a scatter plot or a
        heatmap. Active selections can be superposed.
        '''
        import matplotlib.pyplot as plt
        
        fig = plt.figure()
        if method == 'heatmap':
            plt.hist2d(self.X_1, self.X_2, **kwargs)
        elif method == 'scatter':
            plt.scatter(self.X_1, self.X_2, **kwargs)
        else:
            print('Unknown method.')
            
        if show_selections:
            for S in self.active_selections:
                from matplotlib.patches import Rectangle as RectanglePatch
                if isinstance(S, Rectangle):
                    anchor = (S.x0, S.y0)
                    width = S.x1 - S.x0
                    height = S.y1 - S.y0
                    shape = RectanglePatch(anchor, width, height, **S.kwargs)
                else:
                    print('Unknown shape.')
                plt.gca().add_patch(shape)
        if show:
            plt.show()
        return fig
        
       
class Rectangle:
    
    def __init__(self, x0, y0, x1, y1, label, **kwargs):
        self.x0 = x0
        self.y0 = y0
        self.x1 = x1
        self.y1 = y1
        self.label = label
        self.kwargs = kwargs
        
    def is_inside(self, x, y):
        return self.x0 < x < self.x1 and self.y0 < y < self.y1
        
    def __repr__(self):
        return 'Rectangle=(label="{}", x0={}, y0={}, x1={}, y1={})'.format(self.label,
                                                                           self.x0,
                                                                           self.y0,
                                                                           self.x1,
                                                                           self.y1)