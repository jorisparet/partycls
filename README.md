partycls
========

[![pypi](https://img.shields.io/pypi/v/partycle.svg)](https://pypi.python.org/pypi/partycls/)
[![version](https://img.shields.io/pypi/pyversions/partycls.svg)](https://pypi.python.org/pypi/partycls/)
[![license](https://img.shields.io/pypi/l/partycls.svg)](https://en.wikipedia.org/wiki/GNU_General_Public_License)
[![pipeline](https://framagit.org/coslo/partycls/badges/develop/pipeline.svg)](https://framagit.org/atooms/atooms/badges/develop/pipeline.svg)
[![coverage report](https://framagit.org/coslo/partycls/badges/develop/coverage.svg)](https://framagit.org/atooms/atooms/-/commits/develop)

**partycls** is a Python package for spatio-temporal cluster analysis of interacting particles. It provides descriptors suitable for applications in condensed matter physics and integrates the necessary tools of unsupervised learning into a streamlined workflow. Thanks to a flexible system of filters, it makes it easy to restrict the analysis to a given subset of particles based on arbitrary particle properties.

Quick start
-----------

Here is a simple example that shows how to use partycls to identify grain boundaries in a polycrystalline system. The system configuration is stored in a trajectory file with a single frame. We use the local distribution of bond angles around each particle as a structural fingerprint and perform a clustering using the [K-Means](https://en.wikipedia.org/wiki/K-means_clustering) algorithm. We show the system coloring the particles according to the cluster they belong to.

```python
from partycls import Trajectory, Workflow

traj = Trajectory('grains.xyz')
wf = Workflow(traj, descriptor='ba', clustering='kmeans')
wf.run()
traj[0].show(color='label', backend='ovito')
```

![](https://raw.githubusercontent.com/jorisparet/partycls/master/data/snapshots/grains_labels.png)

The results are also written to a set of files including a labeled trajectory file and additional information on the clustering results. The whole workflow can be easily tuned and customized, check out the [tutorials](https://github.com/jorisparet/partycls/tree/master/tutorial) to see how and for further examples.

We can restrict the analysis to specific a subset of particles by adding filters. Say we have a binary mixture composed of particles with types A and B, and are only interested in the angular correlations of B particles in the left side of the box (with respect to x-axis):

```python
from partycls import Trajectory
from partycls.descriptor import BondAngleDescriptor

traj = Trajectory('trajectory.xyz')
D = BondAngleDescriptor(traj)
D.add_filter("species == 'B'")
D.add_filter("x < 0.0")
D.compute()

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

Features
--------

partycls is designed to accept a large variety of trajectory formats (including custom ones!) either through its built-in trajectory reader or via third-party packages, such as [MDTraj](www.mdtraj.org) and [atooms](https://framagit.org/atooms/atooms). It relies on the [scikit-learn](https://scikit-learn.org) package to perform feature scaling as well as a number of dimensionality reduction and clustering methods. In addition to its native descriptors, partycls supports additional structural descriptors via [DScribe](https://singroup.github.io/dscribe).

Requirements
------------

* [numpy](https://pypi.org/project/numpy/)
* [scikit-learn](https://scikit-learn.org)
* [optional] [mdtraj](https://www.mdtraj.org) (additional trajectory formats)
* [optional] [atooms](https://framagit.org/atooms/atooms) (additional trajectory formats)
* [optional] [dscribe](https://singroup.github.io/dscribe) (additional descriptors)
* [optional] [matplotlib](https://matplotlib.org/) (visualization)
* [optional] [ovito](https://ovito.org/) (visualization)
* [optional] [py3Dmol](https://github.com/avirshup/py3dmol) (interactive 3D visualization)

Documentation
-------------

See the [tutorials](https://github.com/jorisparet/partycls/tree/master/tutorial) (Jupyter notebooks) for a step-by-step introduction to the main features of partycls and some of its applications.

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

Support and contribution
------------------------

If you wish to contribute or report an issue, feel free to [contact us](mailto:joris.paret@umontpellier.fr) or to use the [issue tracker](https://github.com/jorisparet/partycls/issues) and [pull requests](https://github.com/jorisparet/partycls/pulls) from the code repository. 

Authors
-------

Joris Paret

Daniele Coslovich: http://www-dft.ts.infn.it/~coslovich/
