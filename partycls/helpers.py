"""
Various helper functions for visualization, cluster analysis, etc.
"""

import numpy
from sklearn.metrics import adjusted_mutual_info_score as AMI
from sklearn.metrics import adjusted_rand_score as ARI
from sklearn.metrics import silhouette_samples, silhouette_score

__all__ = ['AMI',
           'ARI',
           'silhouette_samples',
           'silhouette_score',
           'show_matplotlib',
           'show_ovito',
           'show_3dmol',
           'shannon_entropy',
           'merge_clusters',
           'sort_clusters']

_palette = ["#50514f", "#f25f5c", "#ffe066", "#247ba0", "#70c1b3",
            "#0cce6b", "#c200fb", "#e2a0ff", "#6622cc", "#119822"]


def hex_to_rgb(h):
    h = h.lstrip('#')
    return tuple(int(h[i:i + 2], 16) / 256 for i in (0, 2, 4))


def show_matplotlib(system, color, view='top', palette=None, cmap='viridis',
                    outfile=None, linewidth=0.5, alpha=1.0, show=False):
    """
    Make a snapshot of the ``system`` using matplotlib.
    The figure is returned for further customization or visualization 
    in jupyter notebooks.    

    Parameters
    ----------
    system : System
        The system to visualize.

    color : str
        Particle property to use for color coding, *e.g.* ``"species"``, ``"label"``.

    view : str, default: "top"
        View type, *i.e.* face of the box to show. Only works for a 3D system.

    palette : list, default: None
        List of colors when coloring particles according to a discrete property,
        such as ``"species"`` or ``"label"``. A default palette will be used if 
        not specified.

    cmap : str, default: "viridis"
        Name of a matplotlib colormap to use when coloring particles according
        to a continuous property such as ``"velocity"`` or ``"energy"``. List of 
        available colormaps can be found in ``matplotlib.cm.cmaps_listed``.

    outfile : str, default: None
        Output filename to save the snapshot. Default is to not save.

    linewidth : float, default: 0.5
        Line width.

    alpha : float, default: 1.0
        Transparency parameter.

    show : bool, default: False
        Show the snapshot when calling the function.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure of the snapshot.
    """
    try:
        import matplotlib.pyplot as plt
        from matplotlib.cm import cmaps_listed
    except ModuleNotFoundError:
        print('install matplotlib to display the particles')
        return
    from .core.utils import tipify
    from numpy import array, sign, argsort

    views = {'top': [1, 2, 3],
             'bottom': [1, -2, -3],
             'front': [1, 3, -2],
             'back': [-1, 3, 2],
             'left': [-2, 3, -1],
             'right': [2, 3, 1]}

    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    ax.set_xlim((-system.cell.side[0] / 2, system.cell.side[0] / 2))
    ax.set_ylim((-system.cell.side[1] / 2, system.cell.side[1] / 2))
    # scale marker size relative to box size
    M = ax.transData.get_matrix()
    scale = M[0, 0]
    # color according to a specific property
    property_vals = system.get_property('particle.{}'.format(color))

    # discrete property?
    discrete = isinstance(tipify(str(property_vals[0])), (str, int))
    # corresponding color system
    discrete_colors = _palette if palette is None else palette
    colormap = cmaps_listed[cmap]
    if discrete:
        property_set = list(set(property_vals))
        property_set.sort()
        color_db = discrete_colors
    else:
        color_db = colormap(list(property_vals))

    # list of individual colors
    colors = []
    for pn, p in enumerate(system.particle):
        if discrete:
            p_color = color_db[property_set.index(p.__getattribute__(color))]
            colors.append(p_color)
        else:
            colors.append(color_db[pn])
    colors = array(colors)

    # positions and radii
    pos = system.get_property('position')
    R = system.get_property('radius')

    # plot 3D
    if system.n_dimensions == 3:
        xi, yi, zi = views[view]
        X = sign(xi) * pos[:, abs(xi) - 1]
        Y = sign(yi) * pos[:, abs(yi) - 1]
        Z = sign(zi) * pos[:, abs(zi) - 1]
        order = argsort(Z)
        sizes = numpy.array((scale * R[order])**2, dtype=numpy.int64)
        ax.scatter(X[order], Y[order], c=colors[order],
                   marker='o', ec='k', s=sizes,
                   linewidths=linewidth, alpha=alpha)
    # plot 2D
    if system.n_dimensions == 2:
        X = pos[:, 0]
        Y = pos[:, 1]
        ax.scatter(X, Y, c=colors, marker='o', ec='k', s=(scale * R)**2,
                   linewidths=linewidth, alpha=alpha)

    if outfile is not None:
        fig.savefig(outfile, bbox_inches='tight')
    if show:
        plt.show()
    return fig


def show_ovito(system, color, view='top', palette=None, cmap='viridis',
               outfile=None, size=(640, 480), zoom=True):
    """
    Make a snapshot of the ``system`` using Ovito.
    The image is returned for further customization or visualization 
    in jupyter notebooks.        

    Parameters
    ----------
    system : System
        The system to visualize.

    color : str
        Particle property to use for color coding, e.g. ``"species"``, ``"label"``.

    view : str, default: "top"
        View type, *i.e.* face of the box to show. Only works for a 3D system.

    palette : list, default: "viridis"
        List of colors when coloring particles according to a discrete property,
        such as ``"species"`` or ``"label"``. Colors must be expressed in RGB format through
        tuples (*e.g.* ``palette=[(0,0,1), (1,0,0)]``). A default palette will be used
        if not specified.
   
    cmap : str, optional
        Name of a matplotlib colormap to use when coloring particles according
        to a continuous property such as ``"velocity"`` or "``energy"``. List of 
        available colormap can be found in ``matplotlib.cm.cmaps_listed``.

    outfile : str, default: None
        Output filename to save the snapshot. The default is to not save.

    size : tuple, default: (640, 480)
        Size of the image to render.

    zoom : bool, default: True
        Zoom on the simulation box.

    Returns
    -------
    Image
        Rendered image.
    """
    try:
        from ovito.io import import_file
        from ovito.vis import Viewport, TachyonRenderer
    except ImportError:
        print('Warning: install `ovito` to display the particles.')
        return
    import os
    import tempfile
    from .core.utils import tipify
    from matplotlib.cm import cmaps_listed

    # discrete property?
    property_vals = system.get_property('particle.{}'.format(color))
    discrete = isinstance(tipify(str(property_vals[0])), (str, int))
    # corresponding color system
    discrete_colors = _palette if palette is None else palette
    colormap = cmaps_listed[cmap]

    # discrete or continuous color field
    if palette is None:
        discrete_colors = _palette
        discrete_colors = [hex_to_rgb(c) for c in discrete_colors]
    else:
        discrete_colors = palette
    # corresponding palette / colormap
    if discrete:
        property_set = list(set(property_vals))
        property_set.sort()
        color_db = discrete_colors
    else:
        color_db = colormap(list(property_vals))
    
    # individual particle color
    for pn, p in enumerate(system.particle):
        if discrete:
            p_color = color_db[property_set.index(p.__getattribute__(color))]
            p.color = p_color
        else:
            p.color = color_db[pn]
    
    # Get a temporary file to write the sample
    fh = tempfile.NamedTemporaryFile('w', suffix='.xyz', delete=False)
    tmp_file = fh.name
    # Self-contained EXYZ dump (it is not clean to use trajectories here)
    fh.write('{}\n'.format(len(system.particle)))
    fh.write('Properties=species:S:1:pos:R:3:radius:R:1:color:R:3 Lattice="{},0.,0.,0.,{},0.,0.,0.,{}"\n'.format(
        *system.cell.side))
    for p in system.particle:
        fh.write('{} {} {} {} {} {} {} {}\n'.format(p.species, *p.position, p.radius, *p.color))
    fh.close()

    # Ovito stuff. Can be customized by client code.
    pipeline = import_file(tmp_file)
    # Ovito seems to ignore the lattice info of exyz file
    # so we forcibly set the cell info here
    pipeline.source.data.cell_[0, 0] = system.cell.side[0]
    pipeline.source.data.cell_[1, 1] = system.cell.side[1]
    pipeline.source.data.cell_[2, 2] = system.cell.side[2]
    pipeline.source.data.cell_[:, 3] = -system.cell.side / 2
    pipeline.add_to_scene()

    views = {'top':    Viewport.Type.Top,
             'bottom': Viewport.Type.Bottom,
             'front':  Viewport.Type.Front,
             'back':   Viewport.Type.Back,
             'left':   Viewport.Type.Left,
             'right':  Viewport.Type.Right}
    vp_type = views[view]
    vp = Viewport(type=vp_type)

    # Render
    if zoom:
        vp.zoom_all()
    if outfile is None:
        outfile = tmp_file + '.png'
        
    vp.render_image(filename=outfile, 
                    size=size, 
                    renderer=TachyonRenderer())

    # Scene is a singleton, so we must clear it
    pipeline.remove_from_scene()
    
    # remove temporary exyz file
    os.remove(tmp_file)

    # Try to display the image (e.g. in a jupyter notebook)
    try:
        from IPython.display import Image
        return Image(outfile)
    except ImportError:
        return outfile


def show_3dmol(system, color, palette=None):
    """
    Visualize the ``system`` using 3dmol.
    The py3Dmol view is returned for further customization or visualization 
    in jupyter notebooks.

    Parameters
    ----------
    system : System
        The system to visualize.

    color : str
        Particle property to use for color coding, *e.g.* ``"species"``, ``"label"``.
        This property must be a string or an integer.

    palette : list, default: None
        List of colors when coloring particles according to a discrete property,
        such as ``"species"`` or ``"label"``. A default palette will be used if not 
        specified.

    Raises
    ------
    ValueError
        If the ``color`` parameter refers to a ``float`` particle property.

    Returns
    -------
    view : py3Dmol.view
        py3Dmol view.
    """
    try:
        import py3Dmol
    except ImportError:
        print('install py3Dmol to display the particles')
        return
    from .trajectory import tipify

    if palette is None:
        palette = _palette
    view = py3Dmol.view()
    view.setBackgroundColor('white')
    # color according to a specific property
    property_vals = system.get_property('particle.{}'.format(color))
    property_set = list(set(property_vals))
    property_set.sort()
    if not isinstance(tipify(str(property_vals[0])), (str, int)):
        raise ValueError('cannot color particle according to a float property')
    # plot particles
    for p in system.particle:
        p_color = palette[property_set.index(p.__getattribute__(color))]
        view.addSphere({'center': {'x': p.position[0],
                                   'y': p.position[1],
                                   'z': p.position[2]},
                        'radius': p.radius,
                        'color': p_color})
    # plot cell
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
    Shannon entropy of distribution :math:`p(x)`.

    Parameters
    ----------
    px : numpy.array
        Distribution :math:`p(x)`.

    dx : float, default: 1.0
        Differential of x.

    Returns
    -------
    S : float
        Shannon entropy.
    """
    S = 0.0
    P = px * dx
    for p in P:
        if p > 0.0:
            S += p * numpy.log(p)
    return -S


def merge_clusters(weights, n_clusters_min=2, epsilon_=1e-15):
    """
    Merge clusters into ``n_clusters_min`` new clusters based on the
    probabilities that particles initially belong to each of the original
    clusters with a certain probability and using an entropy criterion.
    
    See https://doi.org/10.1198/jcgs.2010.08111 (Baudry et al.).

    Parameters
    ----------
    weights : list
        Probabilities that each particle belongs to each cluster.
        If there are :math`N` particles, then the length of the list (or first
        dimension of the array) must be :math:`N`. If there are math:`K` original
        clusters, each element of ``weights`` (or the first dimension of the array)
        must be :math:`K`. ``weights[i][j]`` (list) or ``weights[i,k]`` (array) is the 
        probability that particle ``i`` belongs to cluster ``k`` before merging.
        For each particle, ``sum(weights[i])`` is equal to 1.

    n_clusters_min : int, default: 2
        Final number of clusters after merging.

    epsilon_ : float
        Small number (close to zero). This is needed as a replacement for zero
        when computing a logarithm to avoid errors.

    Returns
    -------
    new_weights : numpy.ndarray
        New weights after merging. Same shape and interpretation as the
        ``weights`` input parameter.

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
                labs.append([i, j])
                # print entropy change for (i,j)
                # print('  {i}  {j}  {:.4f}'.format(i,j,delta_ent))

        best_merge_index = d_evals.index(max(d_evals))
        # print("# best merge: {}, gain={:.4f}".format(labs[best_merge_index],
        #                                              d_evals[best_merge_index]))

        # do the merge...
        for w in new_weights:
            w[labs[best_merge_index][1]] += w[labs[best_merge_index][0]]
            w[labs[best_merge_index][0]] = epsilon_

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

    func : function, default: ``shannon_entropy``
        Function used to associate a numerical value to each cluster, to be
        used as sorting criterion. This function must accept a list or a 
        one dimensional array as parameter (this parameter being the 
        coordinates of a given centroid).

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
    Entropy change on merging two clusters (following Baudry).
    """
    delta_ent = 0.0
    for w in weights:  # each w is a vector of weights (this is a loop over particles)
        w_merge = w[i] + w[j]  # for this particle, add the weights for cluster i and j
        delta_ent += w_merge * numpy.log(w_merge)
        if w[i] > 0.0:
            delta_ent -= w[i] * numpy.log(w[i])
        if w[j] > 0.0:
            delta_ent -= w[j] * numpy.log(w[j])
    return delta_ent


def _compute_ICL_ent(weights, epsilon_):
    """
    ICL entropy from list of weights.
    """
    ICL_ent = 0.0
    for w in weights:
        wts = numpy.maximum(w, epsilon_)
        ICL_ent -= numpy.sum(wts * numpy.log(wts))
    return ICL_ent
