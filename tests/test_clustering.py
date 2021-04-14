#!/usr/bin/env python

import unittest
import os

from pysc.trajectory import Trajectory
from pysc.descriptor import BondAngleDescriptor, RadialDescriptor
from pysc import Optimization, ZScore, PCA, KMeans, CommunityInference

class Test(unittest.TestCase):

    def setUp(self):
        data = os.path.join(os.path.dirname(__file__), '../data/')
        self.traj = Trajectory(os.path.join(data, 'dislocation.xyz'))
        self.cutoffs = [1.25]

    def test_angular_zscore_pca_kmeans(self):
        D = BondAngleDescriptor(self.traj)
        D.cutoffs = self.cutoffs
        X = D.compute()
        scaler = ZScore()
        X = scaler.scale(X)        
        reducer = PCA(n_components=3)
        clustering = KMeans(n_clusters=2, n_init=100)
        clustering.fit(X)
        # check if the dislocated particles are in the same cluster
        labels = list(clustering.labels)
        self.assertEqual(labels[4], labels[6], 'not the expected labels')
        # check if only the dislocated particles are in the same cluster
        self.assertEqual(set(clustering.fractions), set([2/27, 25/27]),
                         'not the expected cluster fractions')
        
        # Same via optimization
        opt = Optimization(self.traj, descriptor='ba', scaling='zscore', clustering='kmeans')
        opt.descriptor.cutoffs = self.cutoffs
        opt.clustering.n_init = 100
        opt.disable_output()
        opt.run()
        # check if both methods give the same result
        self.assertEqual(set(opt.fractions), set(clustering.fractions),
                         'different cluster fractions')

    """
    def test_radial_ci(self):
        D = RadialDescriptor(self.traj, rlim=(0.0, 1.25))
        X = D.compute()
        clustering = CommunityInference(n_init=100)
        clustering.fit(X)
        # check if the dislocated particles are in the same cluster
        labels = list(clustering.labels)
        self.assertEqual(labels[4], labels[6], 'not the expected labels')
        # check if only tqhe dislocated particles are in the same cluster
        self.assertEqual(set(opt.fractions), set([2/27, 25/27]),
                         'not the expected cluster fractions')
        
        # Same via optimization
        opt = Optimization(self.traj, descriptor='gr', clustering='cinf')
        opt.descriptor.rlim = (0.0, 1.25)
        opt.clustering.n_init = 100
        opt.disable_output()
        # TODO: fix warning
        opt.run()
        # check if both methods give the same result
        self.assertEqual(set(opt.fractions), set(clustering.fractions),
                         'different cluster fractions')
    """
 
if __name__ == '__main__':
    unittest.main()
