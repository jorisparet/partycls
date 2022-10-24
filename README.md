<p align="center">
<a href="https://github.com/jorisparet/partycls"><img src="https://github.com/jorisparet/partycls/blob/master/logo/logo.svg" width="250"></a>
</p>

[![pypi](https://img.shields.io/pypi/v/partycls.svg)](https://pypi.python.org/pypi/partycls/)
[![version](https://img.shields.io/badge/python-3.6+-blue.svg)](https://pypi.python.org/pypi/partycls/)
[![license](https://img.shields.io/pypi/l/partycls.svg)](https://en.wikipedia.org/wiki/GNU_General_Public_License)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03723/status.svg)](https://doi.org/10.21105/joss.03723)
[![build](https://github.com/jorisparet/partycls/actions/workflows/build-test.yml/badge.svg)](https://github.com/jorisparet/partycls/actions/workflows/build-test.yml)
![coverage](https://img.shields.io/badge/coverage-89%25-green)
  
**partycls** is a Python package for cluster analysis of systems of interacting particles. By grouping particles that share similar structural or dynamical properties, partycls enables rapid and unsupervised exploration of the system's relevant features. It provides descriptors suitable for applications in condensed matter physics, such as structural analysis of disordered or partially ordered materials, and integrates the necessary tools of unsupervised learning into a streamlined workflow.

Homepage
--------

For more details and tutorials, visit the homepage at: https://www.jorisparet.com/partycls

Quick start
-----------

This quick example shows how to use partycls to identify grain boundaries in a polycrystalline system. The system configuration is stored in a [XYZ](https://en.wikipedia.org/wiki/XYZ_file_format) trajectory file with a single frame. We use the local distribution of bond angles around each particle as a structural descriptor and perform a clustering using the [K-Means](https://en.wikipedia.org/wiki/K-means_clustering) algorithm.

```python
from partycls import Trajectory, Workflow

traj = Trajectory('grains.xyz')
wf = Workflow(traj, descriptor='ba', clustering='kmeans')
wf.run()
traj[0].show(color='label', backend='ovito')
```

![](https://raw.githubusercontent.com/jorisparet/partycls/master/data/snapshots/grains_labels.png)

The results are also written to a set of files including a labeled trajectory file and additional information on the clustering results. The whole workflow can be tuned and customized, check out the [tutorials](https://www.jorisparet.com/partycls/tutorials) to see how and for further examples.

Thanks to a flexible system of filters, partycls makes it easy to restrict the analysis to a given subset of particles based on arbitrary particle properties. Say we have a binary mixture composed of particles with types A and B, and we are only interested in analyzing the bond angles of B particles in a vertical slice:

```python
from partycls import Trajectory
from partycls.descriptors import BondAngleDescriptor

traj = Trajectory('trajectory.xyz')
D = BondAngleDescriptor(traj)
D.add_filter("species == 'B'")
D.add_filter("x > 0.0")
D.add_filter("x < 1.0")
D.compute()

# Angular correlations for the selected particles
print(D.features)
```

We can then perform a clustering based on these structural features and ask for 3 clusters:

```python
from partycls import KMeans

clustering = KMeans(n_clusters=3)
clustering.fit(D.features)
print('Cluster membership of the particles', clustering.labels)
```

Main features
-------------

### Trajectory formats

partycls accepts several trajectory formats (including custom ones) either through its built-in trajectory reader or via third-party packages, such as [MDTraj](www.mdtraj.org) and [atooms](https://framagit.org/atooms/atooms). The code is currently optimized for small and medium system sizes (of order 10⁴ particles). Multiple trajectory frames can be analyzed to extend the structural dataset.

### Structural descriptors

partycls implements various structural descriptors: 

* [Radial descriptor](https://www.jorisparet.com/partycls/tutorials/descriptors/gr.html)
* [Tetrahedral descriptor](https://www.jorisparet.com/partycls/tutorials/descriptors/tetra.html)
* [Bond-angle descriptor](https://www.jorisparet.com/partycls/tutorials/descriptors/ba.html)
* [Smoothed bond-angle descriptor](https://www.jorisparet.com/partycls/tutorials/descriptors/sba.html)
* [Bond-orientational descriptor](https://www.jorisparet.com/partycls/tutorials/descriptors/bo.html)
* [Smoothed bond-orientational descriptor](https://www.jorisparet.com/partycls/tutorials/descriptors/sbo.html)
* [Locally averaged bond-orientational descriptor](https://www.jorisparet.com/partycls/tutorials/descriptors/labo.html)
* [Radial bond-orientational descriptor](https://www.jorisparet.com/partycls/tutorials/descriptors/rbo.html)
* [Compactness descriptor](https://www.jorisparet.com/partycls/tutorials/descriptors/compact.html)
* [Coordination descriptor](https://www.jorisparet.com/partycls/tutorials/descriptors/coord.html)

### Machine learning

partycls performs feature scaling, dimensionality reduction and cluster analysis using the [scikit-learn](https://scikit-learn.org) package and additional built-in algorithms.

Dependencies
------------

partycls relies on several external packages, most of which only provide additional features and are not necessarily required.

### Required

* Fortran compiler (*e.g.* [gfortran](https://gcc.gnu.org/wiki/GFortran))
* [NumPy](https://pypi.org/project/numpy/)
* [scikit-learn](https://scikit-learn.org)

### Optional

* [MDTraj](https://www.mdtraj.org) (additional trajectory formats)
* [atooms](https://framagit.org/atooms/atooms) (additional trajectory formats)
* [DScribe](https://singroup.github.io/dscribe) (additional descriptors)
* [Matplotlib](https://matplotlib.org/) (visualization)
* [OVITO](https://ovito.org/) < 3.7.0 (visualization)
* [Py3DMol](https://github.com/avirshup/py3dmol) (interactive 3D visualization)
* [pyvoro](https://github.com/joe-jordan/pyvoro) or its [memory-optimized fork](https://framagit.org/coslo/pyvoro) for large systems (Voronoi neighbors and tessellation)
* [tqdm](https://tqdm.github.io/) (progress bars)

Documentation
-------------

Check the [tutorials](https://www.jorisparet.com/partycls/tutorials) to see various examples and detailed instructions on how to run the code, as well as an in-depth presentation of the built-in structural descriptors.

For a more detailed documentation, you can check the [API](https://www.jorisparet.com/partycls/api).

Installation
------------

### From PyPI

The latest stable release is available on [PyPI](https://pypi.org/project/partycls/). Install it with `pip`:

```bash
pip install partycls
```

### From source 

To install the latest development version from source, clone the source code from the official [GitHub repository](https://github.com/jorisparet/partycls) and install it with:

```bash
git clone https://github.com/jorisparet/partycls.git
cd partycls
make install
```

Run the tests using:

```bash
make test
```

or manually compile the Fortran sources and run the tests:

```bash
cd partycls/
f2py -c -m neighbors_wrap neighbors.f90
cd descriptor/
f2py -c -m realspace_wrap realspace.f90
cd ../../
pytest tests/
```

Support and contribution
------------------------

If you wish to contribute or report an issue, feel free to [contact us](mailto:joris.paret@gmail.com) or to use the [issue tracker](https://github.com/jorisparet/partycls/issues) and [pull requests](https://github.com/jorisparet/partycls/pulls) from the [code repository](https://github.com/jorisparet/partycls).

We largely follow the [GitHub flow](https://guides.github.com/introduction/flow/) to integrate community contributions. In essence:
1. Fork the repository.
2. Create a feature branch from `master`.
3. Unleash your creativity.
4. Run the tests.
5. Open a pull request.

We also welcome contributions from other platforms, such as GitLab instances. Just let us know where to find your feature branch.

Citing partycls
---------------

If you use partycls in a scientific publication, please consider citing the following article:

*[partycls: A Python package for structural clustering](https://joss.theoj.org/papers/10.21105/joss.03723). Paret et al., (2021). Journal of Open Source Software, 6(67), 3723*

Bibtex entry:
```
@article{Paret2021,
  doi = {10.21105/joss.03723},
  url = {https://doi.org/10.21105/joss.03723},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {67},
  pages = {3723},
  author = {Joris Paret and Daniele Coslovich},
  title = {partycls: A Python package for structural clustering},
  journal = {Journal of Open Source Software}
}
```

Authors
-------

[Joris Paret](https://www.jorisparet.com/)

[Daniele Coslovich](https://www.units.it/daniele.coslovich/)
