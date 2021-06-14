partycls
========

**partycls** is a Python framework for the clustering of condensed matter systems using simple few-body spatial correlations as structural descriptors. It allows to easily create a customized workflow, from reading the input trajectory to the clustering results, through feature scaling and dimensionality reduction methods. Thanks to a flexible system of filters, it makes it easy to restrict the analysis to a given subset of particles based on arbitrary particle properties.

Quick start
-----------

The goal of partycls is to provide a coherent interface to the basic objects of structural clustering. After opening a trajectory file, the process of clustering and analyzing the results comes in the form of a simple and tunable workflow to facilitate the study of the various parameters and comparisons with other clusterings.

In this simple example, we read a trajectory file in XYZ format and perform a clustering based on the radial correlations of the particle, using the radial distribution of particles around a central particle as numerical fingerprint and the K-Means clustering algorithm to form the clusters:

```python
from partycls import Workflow

wf = Workflow('trajectory.xyz', descriptor='gr', clustering='kmeans')
wf.run()
print('Cluster membership of the particles', wf.labels)
```

This will also write a set of files relative the execution of the workflow. Among those, a trajectory file similar to the input trajectory, with an additional row for the cluster membership of each particle.

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

We can then perform a clustering based on these structural features, asking for e.g. 3 clusters:

```python
from partycls import KMeans

clustering = KMeans(n_clusters=3)
clustering.fit(D.features)
print('Cluster membership of the particles', clustering.labels)
```

Documentation
-------------

See the [tutorial](https://github.com/jorisparet/partycls/tree/master/tutorial) in the form of Jupyter notebooks for a step-by-step introduction to partycls objects and main features. 

Installation
------------

From the code reposity

```bash
git clone https://github.com/jorisparet/partycls.git
cd partycls
make install
```

Authors
-------

Joris Paret

Daniele Coslovich: http://www-dft.ts.infn.it/~coslovich/