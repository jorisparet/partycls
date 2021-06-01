---
title: 'PySCL: A Python package for structural clustering'
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

# TODO

- mention freud?
- mention possible applications? (liquid-liquid, molecules...)

# Summary

PySCL is a Python framework to perform a structural clustering of a condensed matter system, that is, grouping particles that share similar local structural environments. Applications include local structure discovery in heterogeneous materials such as polycrystalline and partially ordered solids, as well as in supercooled liquids and glasses. The code provides a coherent workflow for the calculation of structural descriptors and for the common tasks of unsupervised machine learning. Through a simple and expressive interface, PySCL allows one to open a trajectory file, perform a clustering based on the selected structural descriptor, and analyze and/or save the results with only few lines of code. Additional pre-processing steps such as feature scaling and dimensionality reduction are organically integrated into the workflow and make it easy to assess the robustness of the results.

# Statement of need & Design

Analysis of the local arrangements of atoms and molecules in dense liquids and solids is crucial to understand their emergent physical properties.
This is particularly important in systems whose local structure is heterogeneous, which include polycrystalline materials and partially ordered systems, like semi-crystalline polymers [@McIlroy_2017] or metastable liquids during crystal nucleation [@Russo_Tanaka_2016].
Even more challenging is the case of glass-forming liquids and glasses [@Royall_Williams_2015], which often display locally stable arrangements, known as locally favored structures, whose symmetry and local chemical concentration differ in a subtle way from the bulk.
In the glassy systems, structure-property relationships have been long sought, but are difficult to identify in general [@Hocky_Coslovich_Ikeda_Reichman_2014] and require tailored local structural descriptors [@Richard_Ozawa_2020].

Several methods are available to classify the particles according to local arrangements of their neighbors.
Traditional methods include the Voronoi tessellation [@tanemura_geometrical_1977] and common neighbor analysis (CNA) [@honeycutt_molecular_1987], while more recent approaches provide detailed insight into the topology of the particles' arrangements [@malins_identification_2013,@Lazar_Han_Srolovitz_2015]. Many of these methods are implemented in open source code and can be directly applied to trajectories produced by computer simulations, but also to experimental data of colloidal suspensions analyzed using confocal microscopes [@Royall_Williams_2015].
One of the shortcomings of these approaches, however, is that they tend to produce a very large number of distinct signatures, especially in disordered systems.
Moreover, small distortions of the local environments can substantially affect the structural fingerprint of the particles.

Recently, unsupervised learning has emerged as an alternative approach to characterize the local structure of disordered materials [@reinhart_machine_2017;@boattini_unsupervised_2019].
In particular, clustering methods based on simple observables, such as radial distribution functions, bond angle distributions, and bond orientational parameters (BOP), can provide useful insight into the structural heterogenity of glassy systems [@boattini_autonomously_2020;@paret_assessing_2020].
With the present code, we aim to provide a coherent framework to facilitate unsupervised learning of local structure in condensed matter systems.
The idea is to differentiate the particles' structural environments through the prism of clustering, i.e. by grouping the particles according to the similarity of their local structure.
Through a variety of structural descriptors, dimensionality reduction methods, clustering algorithms and filtering options, PySCL makes it possible to customize these steps to study specific aspects of the structure and to assess the robustness of the results.

PySCL provides a simple and configurable workflow, from reading the input trajectory, through the pre-processing steps, to the final clustering results.
It is designed to accept a large variety of formats for trajectory files, by relying on third-party packages such as `MDTraj` [@McGibbon2015MDTraj], which supports several well-known trajectory formats, and `atooms` [@coslovich_daniele_2018_1183302], which makes it easy to interface custom formats often used by in-house simulation codes. Thanks to a flexible system of filters, it is possible to compute the structural descriptors or perform the clustering on restricted subsets of particles of the system, based on arbitrary particle properties. A substantial fraction of the code acts as a wrapper around functions of the machine learning package `scikit-learn` [@scikit-learn]. This allows non-experienced users to rely on the simplicity of PySCL's interface without any prior knowledge of this external package, while experienced users can take full advantage of the many options provided by `scikit-learn`. In addition, the code also integrates a statistical inference method tailored to amorphous materials [@paret_assessing_2020] and several additional helper functions such as cluster merging for mixture models [@baudry_combining_2010] and consistent centroid-based cluster labeling. A simple diagram of the different steps and combinations to create a custom workflow is shown in \autoref{fig:workflow}. A collection of notebooks, with various examples and detailed instructions on how to run the code, is available on [PySCL's repository](https://gitlab.etu.umontpellier.fr/jorisparet/pysc).

![The different steps to perform a structural clustering. The input must be trajectory a file with a supported format. After selecting the type of structural descriptor (and optional filters) to use for the clustering, optional steps for pre-processing the data are possible: feature scaling and dimensionality reduction. Finally, a clustering is performed using the selected algorithm. Output files are written (unless disabled by the user), such as a labeled trajectory file (*i.e.* containing a row with cluster labels, to facilitate visualization) or the dataset used by the clustering algorithm. \label{fig:workflow}](figures/workflow.pdf)

# Examples

As a simple first example, we consider the detection of the grain boundaries in a polycrystal formed by differently oriented FCC crystallites. This can be done, for instance, using a simple radial descriptor, since that the average radial distribution of particles at the boundaries should be different to that of the crystalline ones. Once the grain boundaries and the crystalline domains are identified as two distinct groups of particles, clustering in real-space allows one to identify each individual grain. The following short piece of code opens the input trajectory stored in the file `grains.xyz`, computes the local radial distribution functions of the particles, applies a standard Z-Score normalization on the data, and finally performs a clustering using the Gaussian mixture model (GMM) with $K=2$ clusters (default):

```python
from pysc import Optimization

opt = Optimization('grains.xyz',
                   descriptor='gr',
                   scaling='zscore',
                   clustering='gmm')
opt.run()
```

Each of these steps is easily tunable, so as to change the workflow with little effort. The labels are available as a simple attribute of the optimization instance. Optionally, a set of output files can be produced for further analysis, including a trajectory file with the cluster labels. All this allows one to quickly visualize the nature of the clusters, as shown in \autoref{fig:grains}.

![(a) A polycrystalline material with differently oriented FCC crystallites. (b) Using the individual radial distributions of the particle, we can distinguish between the crystalline particles (blue, $k=0$) and particles at the boundaries (red, $k=1$). (c) The radial distribution functions restricted to the clusters show a clear difference between the two local environments, with higher peaks for the crystals. All 3D visualizations were rendered in OVITO [@ovito]. \label{fig:grains}](figures/grains_figure.pdf)

The local structure of a glass-forming liquid provides a more challenging bench-case, since the system is amorphous overall, but subtle structural features emerge at low temperature. Here, we consider a binary metallic alloy Cu$_{64}$Zr$_{36}$, which shows a tendency for local icosahedral arrangements around copper atoms [@soklaski_locally_2016]. The fraction of atoms that form such locally favored structures increases markedly when the system is cooled at low temperature. We use LAMMPS [@plimpton_fast_1995] to perform a molecular dynamics simulation using an embedded atom potential. After a rapid quench from high temperature, the supercooled liquid is annealed at $T=900$K. In the following piece of code, we open a LAMMPS trajectory using `atooms` as backend, we restrict the analysis to the copper atoms and use bond-angle correlations and the K-Means algorithm to form the clusters: 

```python
from pysc import Trajectory, Optimization
from pysc.descriptor import BondAngleDescriptor

trajectory = Trajectory('cuzr_900K.dat', fmt='lammps', backend='atooms')
descriptor = BondAngleDescriptor(trajectory)
descriptor.add_filter("species == 'Cu'")

opt = Optimization(trajectory,
                   descriptor=descriptor,
                   scaling='zscore',
                   clustering='kmeans')
opt.run()
```

Note that here, we directly access classes for the trajectory and the structural descriptor, and then pass them to the `Optimization` instance. Every step realized during the optimization can thus be done manually by directly instantiating the desired classes, without even the need for an `Optimization` instance.

In \autoref{fig:cuzr}, we see that the distribution of the cluster $k=1$ is similar to one expected for icosahedra, whereas that of the cluster $k=0$ is flatter and thus more disordered. This provides evidence of local structural heterogeneity in the system. Similar results have been obtained using related clustering algorithms for simpler models of glass-forming liquids, based on Lennard-Jones interactions [@boattini_autonomously_2020;@paret_assessing_2020].

![(a) Sample of a copper-zirconium mixture at $T=900$K. Copper atoms are colored orange and zirconium atoms are colored grey. We look at the angular correlations around the copper atoms only (orange). (b) Copper atoms are now colored blue ($k=0$) and red ($k=1$) based on their cluster membership. Zirconium atoms (grey) are discarded from the analysis. (c) Bond angle distributions of the clusters. \label{fig:cuzr}](figures/cuzr_figure.pdf)

# Acknowledgements

Thank you Mum and Dad.

# References
