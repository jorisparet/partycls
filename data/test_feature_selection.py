from partycls import Trajectory, ZScore, PCA, FeatureSelector
from partycls.descriptor import BondAngleDescriptor

# Trajectory
traj = Trajectory('kalj_fcc_N1200.xyz')

# Descriptor
D = BondAngleDescriptor(traj)
D.nearest_neighbors_method = 'SANN'
X = D.compute()

# Feature scaling
scaler = ZScore()
X_s = scaler.scale(X)

# PCA
redux = PCA(n_components=2)
X_red = redux.reduce(X_s)

# THIS MIGHT WORK!
F = FeatureSelector(D, X_red[:,0], X_red[:,1])
F.add_rectangular_selection(-3.6, 1.7, 2.7, 6, fill=False, label='outlier', color='k', ls='--')
F.add_rectangular_selection(-3.6, 6, 2.7, 11, fill=False, label='outlier', color='r', ls='--')
F.add_rectangular_selection(-6.1, -3.2, -2.3, 1.1, fill=False, label='B', color='g', ls='--')
fig = F.show_feature_space(method='scatter', alpha=0.3)

# Output trajectory with labels for selections
# traj[0].show(color='outlier', view='left')