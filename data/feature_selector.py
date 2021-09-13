from partycls import Trajectory, ZScore, PCA, FeatureSelector
from partycls.descriptor import BondAngleDescriptor
import matplotlib.pyplot as plt

# Trajectory
traj = Trajectory('kalj_fcc_N1200.xyz')

# Descriptor
D = BondAngleDescriptor(traj)
X = D.compute()

# Feature scaling
scaler = ZScore()
X_s = scaler.scale(X)

# PCA
redux = PCA(n_components=2)
X_red = redux.reduce(X_s)

# Plot
F = FeatureSelector()
F.add_rectangular_selection(-0.2, -1.2, 18, 4, label='zone_1', fill=False, color='b', ls='--')
F.add_rectangular_selection(-5, -6, 7, 13, label='zone_2', fill=False, color='r', ls='--')
F.show_feature_space(X_red[:,0], X_red[:,1], method='scatter', alpha=0.3)

# THIS MIGHT WORK!
# F = FeatureSelector(D, X_1, X_2)
# F.add_rectangular_selection(-0.2, -1.2, 18, 4, label='zone_1')
# F.add_rectangular_selection(-0.2, -1.2, 18, 4, label='zone_2')