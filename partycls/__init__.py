"""
partycls is a Python package for cluster analysis of systems of interacting particles.
By grouping particles that share similar structural or dynamical properties, partycls 
enables rapid and unsupervised exploration of the system's relevant features.
It provides descriptors suitable for applications in condensed matter physics, 
such as structural analysis of disordered or partially ordered materials, and 
integrates the necessary tools of unsupervised learning into a streamlined workflow.
"""

from .core._version import __version__

from .workflow import Workflow
from .trajectory import Trajectory
from .clustering import *
from .dim_reduction import *
from .feature_scaling import *
from .helpers import *
from .particle import aliases as particle_aliases