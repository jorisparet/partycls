#!/usr/bin/env python

import unittest
import os

from pysc.trajectory import Trajectory
from pysc.descriptor import AngularDescriptor
from pysc.processing import ZScore, PCA, KMeans, CommunityInference
from pysc.optimization import Optimization


class Test(unittest.TestCase):

    def setUp(self):
        data = os.path.join(os.path.dirname(__file__), '../data/')
        self.traj = Trajectory(os.path.join(data, 'kalj_N150.xyz'), first=0, last=10)
        self.cutoffs = [1.45, 1.25, 1.25, 1.075]

    def test_angular_zscore_pca_kmeans(self):
        D = AngularDescriptor(self.traj)
        D.cutoffs = self.cutoffs
        D.add_filter("species == 'B'", group=0)
        X = D.compute()
        scaler = ZScore()
        X = scaler.scale(X)        
        reducer = PCA(n_components=3)
        clustering = KMeans(n_clusters=2, n_init=100)
        clustering.fit(X)
        clustering.labels

        # Same via optimization
        opt = Optimization(self.traj, descriptor='ba', scaling='zscore', clustering='kmeans')
        opt.run()
        opt.labels
        
    def test_angular_ci(self):
        D = AngularDescriptor(self.traj)
        D.cutoffs = self.cutoffs
        D.add_filter("species == 'B'", group=0)
        X = D.compute()
        clustering = CommunityInference(n_init=1)
        clustering.fit(X)
        clustering.labels
        
        # Same via optimization
        opt = Optimization(self.traj, descriptor='ba', clustering='cinf')
        # TODO: fix warning
        opt.descriptor.add_filter("species == 'B'", group=0)
        opt.run()
        opt.labels
        
if __name__ == '__main__':
    unittest.main()
