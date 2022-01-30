# partycls

[![pypi](https://img.shields.io/pypi/v/partycls.svg)](https://pypi.python.org/pypi/partycls/)
[![version](https://img.shields.io/badge/python-3.6+-blue.svg)](https://pypi.python.org/pypi/partycls/)
[![license](https://img.shields.io/pypi/l/partycls.svg)](https://en.wikipedia.org/wiki/GNU_General_Public_License)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03723/status.svg)](https://doi.org/10.21105/joss.03723)

**partycls** is a Python package for cluster analysis of systems of interacting particles. By grouping particles that share similar structural or dynamical properties, partycls enables rapid and unsupervised exploration of the system's relevant features. It provides descriptors suitable for applications in condensed matter physics, such as structural analysis of disordered or partially ordered materials, and integrates the necessary tools of unsupervised learning into a streamlined workflow.

## Tutorials

We provide several tutorials through a collection of Jupyter notebooks, with examples and detailed instructions on how to run the code:

1. [Trajectories](1_trajectory)
2. [Workflow](2_workflow)
3. [Descriptors](3_descriptors)
4. [Going further](4_going_further)

These notebooks can be launched into an interactive session on *Binder*. Currently, all the dependencies (including the optional ones) will be installed before launching the notebooks. The session can thus take a few minutes to start.

```{warning}
If you download the notebooks, be sure to download the content of the [data](https://github.com/jorisparet/partycls/tree/master/data) folder as well before executing them. The data folder should be placed one level above the folder containing the notebooks (`../data/`) for the paths to remain correct.
```

:::{seealso}
See the full [API documentation](https://jorisparet.github.io/partycls/docs/API/) for a more detailed description of the code.
:::