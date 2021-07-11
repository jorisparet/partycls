---
title: 'partycls: A Python package for structural clustering'
tags:
  - Python
  - Fortran
  - clustering
  - structure
  - condensed matter
  - physics
authors:
  - name: Joris Paret
    orcid: 0000-0002-1300-1227
    affiliation: 1
  - name: Daniele Coslovich
    affiliation: 2
affiliations:
 - name: Laboratoire Charles Coulomb (L2C), Université de Montpellier, CNRS, Montpellier, France
   index: 1
 - name: Dipartimento di Fisica, Università di Trieste, Italy
   index: 2
date: 2021
bibliography: paper.bib

---

# Summary

partycls is a Python framework to perform a clustering of structural features on systems composed of interacting particles. It provides descriptors suitable for applications in condensed matter physics and integrates the necessary tools of unsupervised learning, such as dimensionality reduction, into a streamlined workflow. Through a simple and expressive interface, partycls allows one to open a trajectory file, perform a clustering based on the selected structural descriptor, and analyze and save the results with only a few lines of code.

# Statement of need

Analysis of the local arrangements of atoms and molecules in dense liquids and solids is crucial to understand their emergent physical properties.
This is particularly important in systems whose local structure is heterogeneous, which include polycrystalline materials and partially ordered systems, like semi-crystalline polymers [@Ganda_2020] or metastable liquids during crystal nucleation [@Russo_Tanaka_2016].
Even more challenging is the case of glass-forming liquids and glasses [@Royall_Williams_2015], as these systems may display locally favored structures, whose symmetry and chemical concentration differ in a subtle way from the bulk.

Recently, unsupervised learning has emerged as a novel approach to structural analysis of disordered materials [@reinhart_machine_2017;@boattini_unsupervised_2019].
In particular, clustering methods based on simple observables, such as radial distribution functions, bond angle distributions, and bond orientational parameters (BOP), provide useful insight into the structural heterogenity of glassy systems [@boattini_autonomously_2020;@paret_assessing_2020].
By grouping the particles according to the similarity of their local structure, these methods avoid the profileration of distinct structural signatures that affects conventional methods, like Voronoi-based analysis [@tanemura_geometrical_1977] or common neighbor analysis (CNA) [@honeycutt_molecular_1987], as well as topological classification approaches [@malins_identification_2013;@Lazar_Han_Srolovitz_2015].

With the present code, we aim to provide a coherent numerical framework for unsupervised learning of structural and dynamical features of condensed matter systems.
To the best of our knowledge, there is currently no publicly available code that provides and integrates all the necessary tools needed for this kind of analysis.
Through a variety of structural descriptors, dimensionality reduction methods, clustering algorithms and filtering options, partycls makes it possible to discover the key structural features of a system and to assess the robustness of the results.
The code has already been used in the context of a recent publication [@paret_assessing_2020] and can be easily extended.
In particular, future versions will implement clustering in space *and* time, to learn about the dynamics of the system as well.

# Design

partycls is mostly written in Python, with a few parts coded in Fortran 90 for efficiency.
It provides a simple and configurable workflow, from reading the input trajectory, through the pre-processing steps, to the final clustering results.
It accepts a large variety of trajectory file formats, by relying on optional third-party packages such as `MDTraj` [@McGibbon2015MDTraj], which supports several well-known trajectory formats, and `atooms` [@coslovich_daniele_2018_1183302], which makes it easy to interface custom formats often used by in-house simulation codes. Thanks to a flexible system of filters, it is possible to compute the structural descriptors or perform the clustering on restricted subsets of particles of the system, based on arbitrary particle properties. In addition to its native descriptors, partycls also supports additional structural descriptors via DScribe [@dscribe].

Some parts of the code act as a wrapper around functions of the machine learning package `scikit-learn` [@scikit-learn]. This allows non-experienced users to rely on the simplicity of partycls's interface without any prior knowledge of this external package, while experienced users can take full advantage of the many options provided by `scikit-learn`. In addition, the code integrates the relevant tools for distributional clustering, such as a community inference method tailored to amorphous materials [@paret_assessing_2020], and several helper functions, e.g. for merging mixture models [@baudry_combining_2010] and consistent centroid-based cluster labeling. A simple diagram of the different steps and combinations to create a custom workflow is shown in \autoref{fig:workflow}. A collection of notebooks, with various examples and detailed instructions on how to run the code, is available in the [partycls's repository](https://github.com/jorisparet/partycls).

![The different steps to perform a structural clustering. The input is a file written in any of the trajectory formats supported by partycls. After selecting the structural descriptor and optional filters, two key steps for pre-processing the data are possible: feature scaling and dimensionality reduction. Finally, a clustering is performed using the selected algorithm. Several output files are produced for further analysis. \label{fig:workflow}](figures/workflow.pdf)

To maintain a consistent API as the code base evolves, partycls will rigorously follow [semantic versioning](https://semver.org/) as the code is designed for maintainable extension. Future work will focus on the underlooked case of dynamical clustering by implementing time-dependent descriptors for individual particle trajectory.

# Examples

As a simple example, we consider the detection of the grain boundaries in a polycrystal formed by differently oriented FCC crystallites. This is easily achieved even with a simple radial descriptor, since that the average radial distribution of particles at the boundaries is different than the one of the crystal in the bulk. The following short piece of code opens the input trajectory stored in the file `grains.xyz`, computes the local radial distribution functions of the particles, applies a standard Z-Score normalization on the data, and finally performs a clustering using the Gaussian mixture model (GMM) with $K=2$ clusters (default):

```python
from partycls import Workflow

wf = Workflow('grains.xyz',
              descriptor='gr',
              scaling='zscore',
              clustering='gmm')
wf.run()
```

Each of these steps is easily tunable, so as to change the workflow with little effort. The labels are available as a simple attribute of the `Workflow` instance. Optionally, a set of output files can be produced for further analysis, including a trajectory file with the cluster labels. Quick visualization of the clusters, as in \autoref{fig:grains}, is possible within partycls through optional visualization backends.

![(a) A polycrystalline material with differently oriented FCC crystallites. (b) Using the individual radial distributions of the particles as structural descriptor, the algorithm identifies the crystalline domains (blue, $k=0$) and the grain boundaries (red, $k=1$). (c) The radial distribution functions restricted to these two clusters display a marked difference, with higher peaks for the crystals. All 3D visualizations were performed with OVITO [@ovito]. \label{fig:grains}](figures/grains_figure.pdf)

The local structure of a glass-forming liquid provides a more challenging bench-case, since the system is amorphous overall, but subtle structural features emerge at low temperature. Here, we consider a binary metallic alloy Cu$_{64}$Zr$_{36}$, which shows a tendency for local icosahedral arrangements around copper atoms [@soklaski_locally_2016]. The fraction of atoms that form such locally favored structures increases markedly when the system is cooled at low temperature. We use LAMMPS [@plimpton_fast_1995] to perform a molecular dynamics simulation using an embedded atom potential. After a rapid quench from high temperature, the supercooled liquid is annealed at $T=900$K. In the following piece of code, we open a LAMMPS trajectory using `atooms` as backend, we restrict the analysis to the copper atoms and use bond-angle correlations and the K-Means algorithm to form the clusters: 

```python
from partycls import Trajectory, Workflow
from partycls.descriptor import BondAngleDescriptor

trajectory = Trajectory('cuzr_900K.dat', fmt='lammps', backend='atooms')
descriptor = BondAngleDescriptor(trajectory)
descriptor.add_filter("species == 'Cu'")

wf = Workflow(trajectory,
              descriptor=descriptor,
              scaling='zscore',
              clustering='kmeans')
wf.run()
```

Here, we directly access classes for the trajectory and the structural descriptor, and then pass them to the `Workflow` instance. Every step of the workflow can also be performed manually by directly instantiating the desired classes, without creating an instance of `Workflow`.

In \autoref{fig:cuzr}, we see that the distribution of the cluster $k=1$ is similar to what is expected for icosahedral structural environments, whereas that of the cluster $k=0$ is flatter and thus more disordered. This provides evidence of the local structural heterogeneity of the system. Similar results have been obtained using related clustering algorithms for simpler models of glass-forming liquids based on Lennard-Jones interactions [@boattini_autonomously_2020;@paret_assessing_2020].

![(a) A glassy copper-zirconium alloy at $T=900$K. Copper and zirconium atoms are colored in orange and grey, respectively. We focus on the bond-angle distribution around the copper atoms only. (b) Copper atoms are now colored blue ($k=0$) and red ($k=1$) based on their cluster membership. Zirconium atoms (grey) are discarded from the analysis. (c) Bond-angle distributions of the clusters. \label{fig:cuzr}](figures/cuzr_figure.pdf)

# Acknowledgements

We thank Robert Jack for his contribution to the community inference method, the merging of mixture models, and for his helpful comments.

# References
