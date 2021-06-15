partycls
========

**partycls** is a Python package for spatio-temporal cluster analysis of interacting particles. It provides descriptors suitable for applications in condensed matter physics and integrates the necessary tools of unsupervised learning into a streamlined workflow. Thanks to a flexible system of filters, it makes it easy to restrict the analysis to a given subset of particles based on arbitrary particle properties.

Quick start
-----------

The goal of *partycls* is to provide a coherent interface to the basic objects of structural clustering. After opening a trajectory file, the process of clustering and analyzing the results comes in the form of a simple and tunable workflow that facilitates the study of the various parameters as well as comparisons with other clusterings.

In this simple example, we read the XYZ trajectory file of a tridimensional system composed of crystalline domains separated by grains boundaries. Using a built-in visualization tool, we show this system:

```python
from partycls import Trajectory

traj = Trajectory('grains.xyz')
traj.show()
```

![](data/snapshots/grains_species.png)

Using the local distribution of bond angles around each particle in this system as numerical fingerprint, we perform a clustering using the [K-Means](https://en.wikipedia.org/wiki/K-means_clustering) algorithm. We show the system again, this time using the resulting cluster labels for coloring:

```python
from partycles import Workflow

wf = Workflow(traj, descriptor='ba', clustering='kmeans')
wf.run()
traj.show(color='label')
```

![](data/snapshots/grains_labels.png)

Using simple angular correlations, we were able to identify the grain boundaries. This execution will also write a set of files relative the execution of the designed workflow, such as a labeled trajectory file or additional information on the clustering.

We can also choose to restrict the analysis to a specific subset of particles by adding a filter on any particle property. Say we have a binary system composed of particles with types A and B, and are only interested in the angular correlations of B particles in the left side of the box (with respect to x-axis):

```python
from partycls import Trajectory
from partycls.descriptor import BondAngleDescriptor

traj = Trajectory('trajectory.xyz')
D = BondAngleDescriptor(traj)
D.add_filter("species == 'B'")
D.add_filter("x < 0.0")
D.compute()

# List of active filters on particles' properties
print(D.active_filters)
# Angular correlations for the selected particles
print(D.features)
```

We can then perform a clustering based on these structural features, asking for *e.g.* 3 clusters:

```python
from partycls import KMeans

clustering = KMeans(n_clusters=3)
clustering.fit(D.features)
print('Cluster membership of the particles', clustering.labels)
```

In addition to the built-in trajectory reader, *partycls* is designed to accept a large variety of formats for trajectory files by relying on third-party packages such as [MDTraj](www.mdtraj.org) and [atooms](https://framagit.org/atooms/atooms). It also supports additional descriptors thanks to an adapter for the [DScribe](https://singroup.github.io/dscribe) package.

Requirements
------------

* [numpy](https://pypi.org/project/numpy/)
* [scikit-learn](https://scikit-learn.org)
* [optional] [mdtraj](https://www.mdtraj.org) (additional trajectory formats)
* [optional] [atooms](https://framagit.org/atooms/atooms) (additional trajectory formats)
* [optional] [dscribe](https://singroup.github.io/dscribe) (additional descriptors)

Documentation
-------------

See the [tutorial](https://github.com/jorisparet/partycls/tree/master/tutorial) in the form of Jupyter notebooks for a step-by-step introduction to *partycls* objects and main features.

Installation
------------

From pypi:

```bash
pip install partycls
```

From the code repository:

```bash
git clone https://github.com/jorisparet/partycls.git
cd partycls
make install
```

Authors
-------

Joris Paret

Daniele Coslovich: http://www-dft.ts.infn.it/~coslovich/
