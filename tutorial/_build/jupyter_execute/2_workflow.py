#!/usr/bin/env python
# coding: utf-8

# # Workflow
# 
# In this notebook, we will see how to perform a cluster analysis of the system using the `Workflow` class.
# 
# ## The simplest workflow
# 
# The `Workflow` class provides the most straightforward way to perform a spatial clustering of a system.
# 
# To be instantiated, this class only requires one positional argument: the path to a trajectory file in a supported format or an instance of `partycls.trajectory.Trajectory`. See the **Trajectories** notebook for details about trajectories.
# 
# Here we analyze a trajectory obtained from a simulation of a Lennard-Jones binary mixture, known as a Kob-Andersen mixture:

# In[1]:


from partycls import Workflow

wf = Workflow('../data/kalj_N150.xyz')
wf.run()


# This is the most minimalist script: we open a trajectory file through the `Workflow` class, and use the default parameters to perform a clustering.
# 
# Other parameters can be set when creating a workflow:
# * `descriptor`: the type of structural descriptor to use to perform the clustering. This must be a short string or an instance of `partycls.descriptor.StructuralDescriptor`. The default value is set to `"gr"`, which is the symbol of the `RadialDescriptor` class. All compatible string aliases are stored in the class attribute `descriptor_db`.
# * `scaling`: type of feature scaling to apply on the data. Default is `None`, but a short string (e.g. `"zscore"` or `"minmax"`) or an instance of the associated classes is also possible. All compatible string aliases are stored in the class attribute `scaling_db`.
# * `dim_reduction`: dimensionality reduction method. Default value is `None`, but a short string (e.g. `"pca"` or `"tsne"`) or an instance of the associated classes is also possible. All compatible string aliases are stored in the class attribute `dim_redeux_db`.
# * `clustering`: clustering algorithm. Default is `"kmeans"`, but other short strings (`"gmm"` or `"cinf"`) or an instance of the associated classes is also possible. All compatible string aliases are stored in the class attribute `clustering_db`.
# 
# With the three lines of code above, we used the radial distribution of the particles to form clusters using the K-Means algorithm (no feature scaling or dimensionality reduction by default). Let's have a look at the bulk *unnormalized* radial distribution, and the distributions restricted to the clusters:

# In[2]:


import matplotlib.pyplot as plt

# grid and average distribution
r = wf.descriptor.grid
p_r = wf.descriptor.average

# dataset: all the individual radial distributions of the particles
data = wf.descriptor.features
# distributions of the clusters
pk_r = wf.clustering.centroids(data)

# plot
plt.plot(r, p_r, c='k', ls='--', label='bulk')
plt.plot(r, pk_r[0], c='b', label='k=0')
plt.plot(r, pk_r[1], c='r', label='k=1')
plt.xlabel('r')
plt.ylabel('Unnormalized radial distributions')
plt.legend(frameon=False)
plt.show()


# **N.B.** :
# * `descriptor.average` refers to the average feature vector of all the objects in the dataset. In this example, this is the average radial distribution of the particles.
# * `clustering.centroids` are the coordinates of the cluster centers for the dataset `data`. This is the average feature vector of all the objects that belong to this cluster. In this example, the average radial distribution of the clusters.
# 
# We see that the bulk distribution (dashed line) has been split into two quite different distributions. Let's have a look at a more common quantity: the radial distribution function (RDF), $g(r)$.
# 
# The class `RadialDescriptor`, used here to compute the radial correlations, implements a normalization function `normalize` that allows for two different normalization specified by a parameter `method`:
# * `method='r2'` returns $r^2g(r)$ (default) ;
# * `method='gr'` returns the standard $g(r)$ ;
# 
# Let us look at the standard $g(r)$:

# In[3]:


# normalized g(r)
g_r = wf.descriptor.normalize(p_r, method='gr')
g0_r = wf.descriptor.normalize(pk_r[0], method='gr')
g1_r = wf.descriptor.normalize(pk_r[1], method='gr')

# plot normalized g(r)
plt.plot(r, g_r, c='k', ls='--', label='bulk')
plt.plot(r, g0_r, c='b', label='k=0')
plt.plot(r, g1_r, c='r', label='k=1')
plt.xlabel('r')
plt.ylabel('RDF')
plt.legend(frameon=False)
plt.show()


# The numbers, fractions and labels of the clusters are accessible through the following attributes:

# In[4]:


print('Populations:', wf.populations)
print('Fractions;', wf.fractions)
print('Labels:', wf.labels)


# ## Output files
# 
# Unless disabled with `wf.disable_output()`, a set output files is written once the cluster analysis is finished. Some are written by default while some need to be enabled. The output files include:
# * `"trajectory"`: the optimized trajectory file with an additional column for the clusters' labels (if the output format allows it) ;
# * `"log"`: a log file with all the relevant information about the optimization ;
# * `"centroids"`: a data file with the centroids of the clusters, using the raw features from the descriptor ;
# * `"labels"`: a text file with the clusters' labels only ;
# * `"dataset"`: a data file with the raw features from the descriptor ;
# 
# All the information regarding output files is stored in the `output_metadata` attribute (a dictionary) :

# In[5]:


# all possible output files
print(wf.output_metadata.keys())


# Each of these files has its specific options, for example, here are the options for the output trajectory:

# In[6]:


# options for the output trajectory
print(wf.output_metadata['trajectory'].keys())

# default values of some options
print('Enable output:', wf.output_metadata['trajectory']['enable'])
print('Output format:', wf.output_metadata['trajectory']['fmt'])


# This can be changed through the method `set_output_metadata`, for example:

# In[7]:


wf.set_output_metadata('trajectory',
                       filename='my_traj.xyz.gz',
                       fmt='rumd')


# If no `filename` is provided for the output files, a default naming convention will be used. 
# 
# The default naming convention of every file is defined by the attribute `wf.naming_convention`:

# In[8]:


print(wf.naming_convention)


# Each tag in the previous string will be replaced by its value in the current instance of `Worfklow`. This naming convention can be changed by using any combination of the following tags: `{filename}`, `{code}`, `{descriptor}`, `{scaling}`, `{dim_reduction}`, `{clustering}`. For instance, setting
# 
# `wf.naming_convention = "{filename}_{clustering}_{descriptor}"`
# 
# would create output files with names : `kalj_N150.xyz_kmeans_gr.log`, `kalj_N150.xyz_kmeans_gr.dataset`, etc.

# ## A customized workflow
# 
# Let us look at a different example, where the `Workflow` needs to be set differently.
# 
# We consider a cubic crystal in which some particles have been dislocated from their original lattice positions. We use angular correlations, a standard feature scaling method (Z-score) and a linear dimensionality reduction technique (PCA) to find the clusters:

# In[9]:


from partycls.descriptor import BondAngleDescriptor
from partycls import Trajectory, PCA

# trajectory
xtal_traj = Trajectory('../data/dislocation_N8000.xyz')

# bond-angle descriptor
D_ba = BondAngleDescriptor(xtal_traj)
X = D_ba.compute()

# dimensionality reduction method
redux = PCA(n_components=2)

# optimization
xtal_wf = Workflow(xtal_traj,
                   descriptor=D_ba,
                   scaling='zscore',
                   dim_reduction=redux,  
                   clustering='kmeans')
xtal_wf.clustering.n_init = 100
xtal_wf.run()


# **N.B.** the raw features, rescaled features, and features in the reduced space are accesible through the attributes `features`, `scaled_features` and `reduced_features` respectively.
# 
# Let's look at the clusters in the 2D latent space corresponding to the two principal components with largest eigenvalues:

# In[10]:


import numpy as np

# features in latent space and labels
X_red = xtal_wf.reduced_features
labels = xtal_wf.labels

# scatter plot of latent space with clusters' labels
c = np.array(['b', 'r'])
plt.scatter(X_red[:,0], X_red[:,1], c=c[labels], alpha=0.3)
plt.show()


# If we now look at the centroids of the clusters (bond angle distribution functions):

# In[11]:


# grid and centroids
theta = D_ba.grid
C_k = xtal_wf.clustering.centroids(X)

# make centroids normalized bond angle distributions
q_0 = D_ba.normalize(C_k[0], method='pdf')
q_1 = D_ba.normalize(C_k[1], method='pdf')

plt.plot(theta, q_0, c='b', label='k=0')
plt.plot(theta, q_1, c='r', label='k=1')
plt.legend(frameon=False)
plt.show()


# This is what is expected for particles on the lattice and for dislocated particles. It is straightforward to see in real space if the clusters are what they should be:

# In[12]:


# show the system blank
fig = xtal_traj[0].show(show=True)

# show the system with labels
fig = xtal_traj[0].show(color='label', palette=['b', 'r'], show=True)


# Dislocated particles are colored differently from the crystal, as expected.
