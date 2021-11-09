#!/usr/bin/env python
# coding: utf-8

# # Descriptors
# 
# In this notebook, we will take a look at the different structural descriptors and learn how to parametrize them and add filters.
# 
# ## Common properties
# 
# A structural descriptor $S(x)$ is a collection of $N$ individual empirical distribution functions or arbitrary features vectors $\{s_i(\vec{x})\}$, defined over a grid $\{x_j\}$ of $M$ features. This data set is stored in the `features` attribute of each descriptor class as a matrix:
# 
# $$
# \begin{pmatrix}
# s_0(x_0) & s_0(x_1) & ... & s_0(x_M) \\
# s_1(x_0) & s_1(x_1) & ... & s_1(x_M) \\ 
# ...      & ...      & ... & ... \\
# s_N(x_0) & s_N(x_1) & ... & s_N(x_M)
# \end{pmatrix}
# $$
# 
# `features` is None by default and is defined only when the `compute()` method is called.
# 
# The features can be calculated between two arbitrary subsets of particles called "groups":
# - group 0 is the main group, *i.e.* particles for which the features are being calculated ;
# - group 1 is the secondary group, *i.e.* particles that are being considered when calculating the features ;
# These groups are formed by adding filters on particles' properties (species, radius, position, etc.).
# 
# All descriptor classes inherit from the abstract class `StructuralDescriptor`, thus having several common properties and attributes such as:
# * `trajectory` : a descriptor is defined for a given `Trajectory` (see the **Trajectories** notebook) ;
# * `dimension` : the spatial dimension of the descriptor. Indeed, some descriptors can only be computed for 2D or 3D trajectories ;
# * `grid` : grid over which the structural features will be computed. For instance, the values of interparticle distances $r$ (*i.e.* bins) when using a radial descriptor ;
# * `active_filters` : list of all active filters on the both groups of particles ;
# * `features` : array of all the computed structural features for particles in the main group ($i.e.$ group 0) ;
# * `size` : number of individual correlation functions $N$ ;
# * `n_features` : the number of features $M$ for the descriptor (equivalent to the length of `grid` and `features.shape[1]`) ;
# * `average` : average feature vector of the descriptor, *i.e.* $\sum_i s_i(\vec{x})$ ;
# 
# Let us consider a pratical example by creating a structural descriptor on a trajectory file:

# In[1]:


from partycls import Trajectory
from partycls.descriptor import StructuralDescriptor

# Open a trajectory
traj = Trajectory('../data/kalj_N150.xyz', last=20)

# Create a descriptor
D = StructuralDescriptor(traj)


# Note that the class `StructuralDescriptor` is abstract and does not perform any computation. This is only an example to show the common properties of the subclasses. Actual examples will be shown later in this notebook.
# 
# By default, all particles are included in groups 0 and 1. The composition of these groups is accessible through the attribute `groups`, a tuple with two elements: 
# - `groups[0]` is the list of particles in the main group (group #0)
# - `groups[1]` the particles in the secondary group (group #1)

# In[2]:


# A tuple with two elements
print(len(D.groups))


# Each element of this tuple is a list whose length is equal to the number of frames in the trajectory:

# In[3]:


# Number of frames in the trajectory
print('Frames :', len(traj))
# Length of the groups
print('Groups :', len(D.groups[0]), len(D.groups[1]))


# Each of these sublists contains all the particles that belong to the group:

# In[4]:


# The first 3 particles in the first frame of group 0
print(D.groups[0][0][0:3])

# All particles are in both groups by default.
# Number of particles in the groups of the first frame
print('Group sizes (first frame) :', len(D.groups[0][0]), len(D.groups[1][0]))


# ### Filters
# 
# To restrict the computation of the descriptor to a given subset, we apply a filter with the method `add_filter`. This is done by specifying the filter with a string of the form:
# 
# `"particle.<attribute> _operator_ <value>"`
# 
# or simply
# 
# `"<attribute> _operator_ <value>"`
# 
# for common attributes. Ex. `"species == 'A'"`, `"particle.radius < 0.5"`.
# 
# We must also specify the `group` over which this filter is applied.
# 
# The current input trajectory file is a binary mixture of particles labeled `"A"` and `"B"` with chemical fractions of 0.8 and 0.2 respectively. Let us restrict the computation of the structural features to B particles only:

# In[5]:


# Add a filter on group 0
D.add_filter("species == 'B'", group=0)

# Active filters
print(D.active_filters)


# Group #0 now contains only B particles, and group #1 still contains all the particles:

# In[6]:


# Group sizes in the first frame
print('Size of group #0 :', len(D.groups[0][0]), '(first frame)')
print('Size of group #1 :', len(D.groups[1][0]), '(first frame)')

# Group sizes over all frames
print('\nSize of group #0 :', D.group_size(0), '(total)')
print('Size of group #1 :', D.group_size(1), '(total)')

# Fractions of particles in the groups
#  (over the whole trajectory)
print('\nGroup #0 fraction :', D.group_fraction(0), '(total)')
print('Group #1 fraction :', D.group_fraction(1), '(total)')


# This means that the features will be computed for B particles only, but by considering all particles (both A and B). Now, let us make this a partial correlation by focusing the analysis on A particles:

# In[7]:


D.add_filter("species == 'A'", group=1)
print(D.active_filters)

# New fractions
print('Group #0 fraction :', D.group_fraction(0))
print('Group #1 fraction :', D.group_fraction(1))


# If the descriptor was describing radial correlations, this would be equivalent to computing the partial RDF between B and A particles, $g_{BA}(r)$.
# 
# **N.B.** We can cumulate filters by adding them successively (as long as there are particles left in the specified group).
# 
# Filters can be removed with the methods `clear_filters` (clears all filters on the specified group) and `clear_all_filters` (clears all filters on both groups):

# In[8]:


# Remove filter on group #0
D.clear_filters(0)
print('Group #0 fraction :', D.group_fraction(0))

# Remove all filters (both groups #0 and #1)
D.clear_all_filters()
print('Group #1 fraction :', D.group_fraction(1))


# All particles are included again in both groups.
# 
# ### Accessing a group property
# 
# Much like the method `Trajectory.get_property` (see the **Trajectories** notebook), the properties of the particles in each group can be accessed through the method `get_group_property` (also called `dump` for short):

# In[9]:


# Species of the particles in group #0
species = D.get_group_property('species', 0)
# 10 first particles in the first frame
print(species[0][0:10])

# Keep only B particles in group #0
D.add_filter("species == 'B'")
# 10 first particles in the first frame
species = D.get_group_property('species', 0)
print(species[0][0:10])


# A few more examples:

# In[10]:


# Positions of particles in group #0
pos_0 = D.dump('pos', 0)

# Radii of particles in group #1
rad_1 = D.dump('radius', 1)

# Position x of particles in group #0
x_0 = D.dump('x', 0)


# **N.B** note that, unlike `Trajectory.get_property`, this method only applies to particles. Moreover, the groups being a data structure natively related to filters, it *does not* accept a `subset` parameter.

# --------------------------------------------------------------------
# 
# **<span style="color:red"> IMPORTANT:</style>** mastering the notion of groups (their data structure, etc.) is not essential to use the code. This is only useful for specific usage, as will become clear with the rest of this notebook and the next.
# 
# Now, in order to actually compute the structural features, the method `compute` must be called: `D.compute()`. This will set the `features` attribute. Let us move to more practical examples with actual descriptors.
# 
# ## Radial descriptor
# 
# The class `RadialDescriptor` computes the radial correlations of a central particle with all its surrounding particles up to a given distance $R$. In practice, for a particle $i$, it computes a histogram for particles between $r$ and $r + \mathrm{d}r$ between the specified range:

# In[11]:


from partycls.descriptor import RadialDescriptor

# radial descriptor between r=0 and r=2 with bin size dr=0.1
D_r = RadialDescriptor(traj, bounds=(0,2), dr=0.05)
# compute the features (also sets D_r.features)
X_r = D_r.compute()

# grid
print('grid = ', D_r.grid)
# radial histogram of particle #0
print('radial histogram of particle #0 :', X_r[0])


# The average of the descriptor is the average over all histograms:

# In[12]:


import matplotlib.pyplot as plt

plt.plot(D_r.grid, D_r.average)
plt.show()


# Note that the bounds for the grid can be set automatically through the integer parameter `n_shells` (instead of specifying `bounds` manually). This will set the upper bound automatically by counting the number of coordination shells (default value is 3).
# 
# An example with a filter:

# In[13]:


# A-A correlation
D_r_AA = RadialDescriptor(traj, bounds=(0,2), dr=0.05)
D_r_AA.add_filter("species == 'A'", 0)
D_r_AA.add_filter("species == 'A'", 1)
_ = D_r_AA.compute()

# A-B correlation
D_r_AB = RadialDescriptor(traj, bounds=(0,2), dr=0.05)
D_r_AB.add_filter("species == 'A'", 0)
D_r_AB.add_filter("species == 'B'", 1)
_ = D_r_AB.compute()


# Let us normalize the distributions in order to get true $g(r)$.
# 
# As mentioned in the **Workflow** notebook, the class `RadialDescriptor` implements a normalization function that allows for two types of normalization through a `method` parameter:
# * `method='r2'` returns $r^2g(r)$ (default) ;
# * `method='gr'` returns the standard $g(r)$ ;
# 
# In our case, we are interested in standard $g(r)$ and $g_{\alpha\beta}(r)$:

# In[14]:


# /!\ normalization depends on the descriptor in this case
g    = D_r.normalize(D_r.average, method='gr')
g_AA = D_r_AA.normalize(D_r_AA.average, method='gr')
g_AB = D_r_AB.normalize(D_r_AB.average, method='gr')

# plot
plt.plot(D_r.grid, g, label='total')
plt.plot(D_r_AA.grid, g_AA, label='A-A')
plt.plot(D_r_AB.grid, g_AB, label='A-B')
plt.legend(frameon=False)
plt.show()


# ## Bond angle descriptor
# 
# Similarly, the class `BondAngleDescriptor` counts the number of particles $(j,k)$ around a central particle $i$ such that the angle $\widehat{jik}$ is in the range $[\theta, \theta + \mathrm{d} \theta$]. Particles $j$ and $k$ must be in the first coordination shell of particle particle $i$.
# 
# This class inherits from a class `AngularStructuralDescriptor` (that inherits from `StructuralDescriptor`). It has several additional attributes:
# 
# * `nearest_neighbors_method` : must be `"FC"` (fixed cutoff) or `"SANN"` (solid angle based nearest-neighbors). See the documentation for more details ;
# * `cutoffs` : cutoffs can be set manually when using `nearest_neighbors_method = "FC"`, or through the method `set_cutoff`. If not specified, they are computed automatically ;
# * `neighbors` : list of nearest neighbors for each particle (in group 0). This is set once the method `nearest_neighbors` has been called ;
# 
# **N.B.** When using `"FC"` as nearest neighbors method and letting the code find the cutoffs automatically, it may happen that some cutoffs do not coincide with the actual cutoff $r_{cut}^{\alpha\beta}$ for a pair $(\alpha,\beta)$ of species. Indeed, if the shape of $g_{\alpha\beta(r)}$ is peculiar (*e.g.* very ordered systems or splitting of the first peak), the identification of the first coordination shell may be wrong. In this case, it is recommended to check the values of `cutoffs`.

# In[15]:


from partycls.descriptor import BondAngleDescriptor

# bond angle descriptor with bin size dtheta=4 (in degrees)
D_ba = BondAngleDescriptor(traj, dtheta=4.0)
# compute the features
X_ba = D_ba.compute()

# cutoffs used to find neighbors
print('all pairs of species :', traj[0].pairs_of_species)
print('associated cutoffs :', D_ba.cutoffs)

# grid (from 0 to 180 deg.)
print('\ngrid = ', D_ba.grid)
# bond angles histogram of particle #0
print('radial histogram of particle #0 :', X_ba[0])
# neighbors of particles #0
print('neighbors of particle #0 :', D_ba.neighbors[0][0])


# Average of the descriptor:

# In[16]:


plt.plot(D_ba.grid, D_ba.average)
plt.show()


# The class `BondAngleDescriptor` also implements a normalization function with a `method` parameter. :
# * `method='sin'`: by construction, the probability density of $\theta$ has a sinusoidal shape in 3D for uniformly distributed points on a sphere (default) ;
# * `method='pdf'` : gives a flat probability density for uniformly distributed points on a sphere ;

# In[17]:


theta = D_ba.grid
q_sin = D_ba.normalize(D_ba.average, method='sin')
q_flat = D_ba.normalize(D_ba.average, method='pdf')

plt.plot(theta, q_sin, label='sin')
plt.plot(theta, q_flat, label='pdf')
plt.legend(frameon=False)
plt.show()


# ## Bond orientational parameters
# 
# The class `BondOrientationalDescriptor` computes bond orientational parameters $q_l$ (see [Steinhardt et al.](http://doi.org/10.1103/PhysRevB.28.784)) around particle $i$ surrounded by nearest neighbors $\{j\}$.
# 
# This class also inherits from the class `AngularStructuralDescriptor` (see the section **Bond angle descriptor** or the documentation for more details).
# 
# The minimal and maximal orders of $l$ must be specified through the parameters `lmin` and `lmax`. Otherwise, specific values of $l$ can be computed by using the list parameter `orders`:

# In[18]:


from partycls.descriptor import BondOrientationalDescriptor

# BOP descriptor
D_bo = BondOrientationalDescriptor(traj, lmin=1, lmax=8)
print('range of values for l :', D_bo.grid)

# compute specific values instead
D_bo.orders = [4,6]
print('specific values of l :', D_bo.grid)

# change nearest neighbors method and compute
D_bo.nearest_neighbors_method = 'SANN'
X_bo = D_bo.compute()


# Distributions of $q_4$ and $q_6$:

# In[19]:


# q_4 (/!\ FIRST element of the grid)
plt.figure(1)
plt.hist(X_bo[:,0], ec='k', bins=20, density=True)
plt.title(r'$q_4$')
plt.show()

# q_6 (/!\ SECOND element of the grid)
plt.figure(1)
plt.hist(X_bo[:,1], ec='k', bins=20, density=True)
plt.title(r'$q_6$')
plt.show()


# ## Lechner-Dellago bond orientational parameters
# 
# The class `LechnerDellagoDescriptor` computes a variant $\bar{q}_l$ (as introduced by [Lechner & Dellago](http://doi.org/10.1063/1.2977970)) of the previously mentioned bond orientational parameters, $q_l$. This class inherits from `BondOrientationalDescriptor` and thus has the same parameters and attributes.
# 
# Note that this descriptor is more computationally expensive than its parent class.

# In[20]:


from partycls.descriptor import LechnerDellagoDescriptor

# Lechner-Dellago descriptor
D_ld = LechnerDellagoDescriptor(traj, orders=[4,6])
# compute specific of values of l
print('grid of l :', D_ld.grid)

# change nearest neighbors method and compute
D_ld.nearest_neighbors_method = 'SANN'
X_ld = D_ld.compute()


# In[21]:


# q_4 (/!\ FIRST element of the grid)
plt.figure(1)
plt.hist(X_ld[:,0], ec='k', bins=20, density=True)
plt.title(r'$\bar{q}_4$')
plt.show()

# q_6 (/!\ SECOND element of the grid)
plt.figure(1)
plt.hist(X_ld[:,1], ec='k', bins=20, density=True)
plt.title(r'$\bar{q}_6$')
plt.show()


# ## Interfacing with DScribe
# 
# Several common machine learning descriptors, such as SOAP (smooth overlap of atomic positions) or ACSF (atom-centered symmetry functions), can be computed using an adapter for the [DScribe](https://singroup.github.io/dscribe/latest/) package.
# 
# This is done using two classes:
# - `DscribeDescriptor` : does not account for the chemical information. Essentially, all the particles are considered as hydrogen atoms ;
# - `DscribeChemicalDescriptor` : accounts for the chemical information. Be sure that the species in the input trajectory are labeled with atomic symbols (`"C"`, `"H"`, `"O"`, `"Si"`) instead of labels such as `"1"`, `"2"` and `"A"`, `"B"`.
# 
# As a proof of principle, let us perform two distinct clusterings using these classes and the SOAP descriptor. We will use a sample of a Wahnström polydisperse binary mixture: 50% large particles ("A"), 50% small particles ("B").

# In[22]:


from partycls import Trajectory, Workflow
from partycls.descriptor import DscribeDescriptor, DscribeChemicalDescriptor
from dscribe.descriptors import SOAP

# Trajectory
traj = Trajectory('../data/wahn_N1000.xyz')

# Change the radius sizes (for visualization purpose)
traj.set_property('radius', 0.5, "species == 'A'")
traj.set_property('radius', 0.4, "species == 'B'")


# ### Without chemical information

# In[23]:


# SOAP descriptor
D = DscribeDescriptor(traj, SOAP, sigma=0.1, rcut=3.0, lmax=7, nmax=7, rbf='gto')

# Workflow
# Z-score + PCA with 2 components (default) + K-Means
wf = Workflow(traj, D, scaling='z-score', 
              dim_reduction='PCA', 
              clustering='K-means')
wf.n_init = 10
wf.run()


# We expect the clustering to partition the particles based on their species (A and B), purely based on spatial correlations since the `DscribeDescriptor` does not account for the chemistry.
# 
# Let us look at two snapshots of the sample: in the first we shall color based on the species, in the second based on cluster labels:

# In[24]:


# First snapshot: species
_ = traj[0].show(color='species', show=True)

# Second snapshot: labels
_ = traj[0].show(color='label', show=True)


# We see that, for the vast majority of the particles, the clustering was able to recover the species solely based on spatial correlation.
# 
# **N.B.** the colors in the second snapshot depend on the execution (labels 0 and 1 can be reversed).

# ### With chemical information
# 
# Let us now repeat the previous clustering using chemical information. As mentioned before, the species of the particles should be atomic symbols. Let us first change the name of the species labels, say `"A"`-->`"C"` and `"B"`-->`"H"` :

# In[25]:


# Change the name of the species
traj.set_property('species', 'C', "species == 'A'")
traj.set_property('species', 'H', "species == 'B'")


# Let us now repeat the clustering using `DscribeChemicalDescriptor`:

# In[26]:


# Descriptor
D = DscribeChemicalDescriptor(traj, SOAP, sigma=0.1, rcut=3.0, lmax=7, nmax=7, rbf='gto')

# Workflow
# Z-score + PCA with 2 components (default) + K-Means
wf = Workflow(traj, D, scaling='z-score', 
              dim_reduction='PCA', 
              clustering='K-means')
wf.n_init = 10
wf.run()


# This time again, a quick way to visualize the results is to visually compare species and cluster labels:

# In[27]:


# First snapshot: species
_ = traj[0].show(color='species', show=True)

# Second snapshot: labels
_ = traj[0].show(color='label', show=True)


# The first clustering was already able to distinguish between the species, so this comes as no surprise that this is still the case here.
# 
# **N.B.** like the previous clustering, the colors in the second snapshot depend on the execution (labels 0 and 1 can be reversed).
