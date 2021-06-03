import numpy

__all__ = ['show_matplotlib',
           'show_3dmol',
           'shannon_entropy',
           'merge_clusters',
           'sort_clusters']

def show_matplotlib(system, color, cmap='viridis', outfile=None, linewidth=0.5, alpha=1.0, show=False):
    """
    Make a snapshot of the `system` using matplotlib.
    The figure is returned for further
    customization or visualization in jupyter notebooks.
    """
    import matplotlib.pyplot as plt
    from .trajectory import tipify
    from matplotlib.cm import cmaps_listed

    discrete_colors = ['r', 'b', 'g', 'y']
    colormap = cmaps_listed[cmap]
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    ax.set_xlim((-system.cell.side[0]/2, system.cell.side[0]/2))
    ax.set_ylim((-system.cell.side[1]/2, system.cell.side[1]/2))
    # check for radius property
    # set default value if non-existant
    if not hasattr(system.particle[0], 'radius'):
        system.set_property('radius', 0.5)
    # color according to a specific property
    property_vals = system.dump('particle.{}'.format(color))
    # discrete property
    discrete = isinstance(tipify(str(property_vals[0])), (str, int))
    if discrete:
        property_set = list(set(property_vals))
        property_set.sort()
        color_db = discrete_colors
    else:
        color_db = colormap(property_vals)
    for pn, p in enumerate(system.particle):
        if discrete:
            p_color = color_db[property_set.index(p.__getattribute__(color))]
        else:
            p_color = color_db[pn]
        c = plt.Circle(p.position[: 2], p.radius,
                       facecolor=p_color, edgecolor='k',
                       alpha=alpha, linewidth=linewidth,
                       zorder=p.position[2])
        ax.add_artist(c)
    if outfile is not None:
        fig.savefig(outfile, bbox_inches='tight')
    if show:
        plt.show()
    return fig

def show_3dmol(system, color, radius=1.0, palette=None):
    """
    Visualize the `system` in cell using 3dmol http://3dmol.csb.pitt.edu/
    """
    import py3Dmol
    from .trajectory import tipify
    
    if palette is None:
        palette = ["#50514f", "#f25f5c", "#ffe066", "#247ba0", "#70c1b3",
                   "#0cce6b", "#c200fb", "#e2a0ff", "#6622cc", "#119822"]
    view = py3Dmol.view()
    view.setBackgroundColor('white')
    # check for radius property
    # set default value if non-existant
    if not hasattr(system.particle[0], 'radius'):
        system.set_property('radius', 0.5)
    # color according to a specific property
    property_vals = system.dump('particle.{}'.format(color))
    property_set = list(set(property_vals))
    property_set.sort()
    if not isinstance(tipify(str(property_vals[0])), (str, int)):
        raise ValueError('cannot color particle according to a float property')
    for p in system.particle:
        p_color = palette[property_set.index(p.__getattribute__(color))]
        view.addSphere({'center': {'x': p.position[0], 'y': p.position[1], 'z': p.position[2]},
                        'radius': radius * p.radius, 'color': p_color})
    if system.cell is not None:
        view.addBox({'center': {'x': 0.0,
                                'y': 0.0,
                                'z': 0.0},
                     'dimensions': {'w': system.cell.side[0],
                                    'h': system.cell.side[1],
                                    'd': system.cell.side[2]},
                     'wireframe': True, 'color': "#000000"})
    return view

def shannon_entropy(px, dx=1.0):
    """
    Shannon entropy of distribution p(x).

    Parameters
    ----------
    px : list or numpy.array
        Distribution p(x).
    dx : float, optional
        Differential of x. The default is 1.0.

    Returns
    -------
    S : float
        Shannon entropy.
    """
    S = 0.0
    P = px * dx
    for p in P:
        if p != 0.0:
            S += p * numpy.log(p)
    return -S

def merge_clusters(weights, n_clusters_min=2, epsilon_=1e-15):
    """
    Merge clusters into `n_clusters_min` new clusters based on the
    probabilities that particles initially belong to each of the original
    clusters with a certain probability and using an entropy criterion.
    
    See https://doi.org/10.1198/jcgs.2010.08111 (Baudry et al.)

    Parameters
    ----------
    weights : list or numpy.ndarray
        Probabilities that each particle belongs to each cluster.
        If there are N particles, then the length of the list (or first
        dimension of the array) must be N. If there are K original clusters,
        each element of `weights` (or the first dimension of the array) must
        be K. `weights[i][j]` (list) or `weights[i,k]` (array) is the 
        probability that particle `i` belongs to cluster `k` before merging.
        For each particle, sum(weights[i]) = 1.
    n_clusters_min : int, optional
        Final number of clusters after merging. The default is 2.
    epsilon_ : float
        Small number (close to zero). This is needed as a replacement for zero
        when computing a logarithm to avoid errors. The default is 1e-15.

    Returns
    -------
    new_weights : numpy.ndarray
        New weights after merging. Same shape and interpretation as the
        `weights` input parameter.
    new_labels : list
        New discrete labels based on the weights after merging.
    """
    new_weights = numpy.asarray(weights).copy()

    # number of clusters (to be changed during the merge)
    n_clusters = weights[0].size

    # print("# i  j  delta_entropy")
    while n_clusters > n_clusters_min:
        # loop over all pairs of clusters and consider what happens if we merge
        d_evals = []
        labs = []
        for i in range(new_weights[0].size):  # i loops over clusters
            for j in range(i):
                delta_ent = _compute_delta_ent(i, j, new_weights)
                d_evals.append(delta_ent)
                labs.append([i,j])
                # print entropy change for (i,j)
                # print('  {i}  {j}  {:.4f}'.format(i,j,delta_ent))

        best_merge_index = d_evals.index(max(d_evals))
        # print("# best merge: {}, gain={:.4f}".format(labs[best_merge_index],
        #                                              d_evals[best_merge_index]))

        # do the merge...
        for w in new_weights:
            w[ labs[best_merge_index][1] ] += w[ labs[best_merge_index][0] ]
            w[ labs[best_merge_index][0] ] = epsilon_

        # print("# ICL entropy before (total) : {:.4f}".format(_compute_ICL_ent(weights, epsilon_)))
        # print("# ICL entropy after  (total) : {:.4f}".format(_compute_ICL_ent(new_weights, epsilon_)))

        n_clusters -= 1
        
    # new discrete labels based on the new weights
    new_labels = [0 for n in range(new_weights.shape[0])]
    for i, wi in enumerate(new_weights):
        new_labels[i] = numpy.argmax(wi)
        
    return new_weights, new_labels

def sort_clusters(labels, centroids, func=shannon_entropy):
    """
    Make a consistent labeling of the clusters based on their centroids by
    computing an associated numerical value as sorting criterion. By default, 
    the labeling is based on the Shannon entropy of each cluster.

    Parameters
    ----------
    labels : list
        Original labels.
    centroids : numpy.ndarray
        Cluster centroids.
    func : function, optional
        Function used to associate a numerical value to each cluster, to be
        used as sorting criterion. This function must accept a list or a 
        one dimensional array as parameter (this parameter being the 
        coordinates of a given centroid). The default is shannon_entropy.

    Returns
    -------
    new_labels : list
        New labels based on centroid entropies.
    new_centroids : numpy.ndarray
        Centroids arranged in order of descending entropies.

    """
    n_clusters = centroids.shape[0]
    new_centroids = numpy.zeros_like(centroids)
    new_labels = labels.copy()
    entropies = numpy.empty(n_clusters)
    # Shannon entropy of each cluster based on centroids
    for k in range(n_clusters):
        entropies[k] = func(centroids[k])
    # rank distributions according to their entropies
    entropies *= -1
    ranks = numpy.argsort(numpy.argsort(entropies))
    # Change particles
    for i in range(len(labels)):
        k_i = labels[i]
        k_new = ranks[k_i]
        new_labels[i] = k_new
    # change fractions, hist and dist
    for k in range(n_clusters):
        k_new = ranks[k]
        new_centroids[k_new] = centroids[k]
    return new_labels, new_centroids

def _compute_delta_ent(i, j, weights):
    """
    Entropy change on merging two clusters (following Baudry) 
    """
    delta_ent = 0.0
    for w in weights:  # each w is a vector of weights (this is a loop over particles)
        w_merge = w[i] + w[j]  # for this particle, add the weights for cluster i and j
        delta_ent += w_merge * numpy.log(w_merge)
        delta_ent -= w[i] * numpy.log(w[i])
        delta_ent -= w[j] * numpy.log(w[j])
    return delta_ent

def _compute_ICL_ent(weights, epsilon_):
    """
    ICL entropy from list of weights
    """
    ICL_ent = 0.0
    for w in weights:
        wts = numpy.maximum(w, epsilon_)
        ICL_ent -= numpy.sum( wts * numpy.log(wts) )
    return ICL_ent