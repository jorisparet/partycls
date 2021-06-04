#!/usr/bin/env python

import unittest
import os

from pysc import Trajectory
from pysc.descriptor import BondAngleDescriptor
from pysc import Workflow, ZScore, PCA, KMeans

class Test(unittest.TestCase):

    def setUp(self):
        data = os.path.join(os.path.dirname(__file__), '../data/')
        self.traj = Trajectory(os.path.join(data, 'dislocation.xyz'), fmt='xyz')
        self.cutoffs = [1.45, 1.25, 1.25, 1.075]

    def test_xyz(self):
        data = os.path.join(os.path.dirname(__file__), '../data/')
        traj = Trajectory(os.path.join(data, 'dislocation.xyz'), fmt='xyz')
        self.assertEqual(len(traj[0].particle), 27)

    def test_angular_zscore_pca_kmeans(self):
        D = BondAngleDescriptor(self.traj)
        D.cutoffs = self.cutoffs
        D.add_filter("species == 'A'", group=0)
        X = D.compute()
        scaler = ZScore()
        X = scaler.scale(X)        
        reducer = PCA(n_components=3)
        Xred = reducer.reduce(X)
        clustering = KMeans(n_clusters=2, n_init=100)
        clustering.fit(Xred)
        print('Fractions :', clustering.fractions, '(clustering alone)')
        
        # Same via workflow
        wf = Workflow(self.traj, descriptor='ba', scaling='zscore',
                      dim_redux='pca', clustering='kmeans')
        wf.descriptor.add_filter('species == "A"', group=0)
        wf.descriptor.cutoffs = self.cutoffs
        wf.dim_redux.n_components = 3
        wf.clustering.n_init = 100
        wf.disable_output()
        wf.run()
        print('Fractions :', clustering.fractions, '(via workflow)')
        self.assertEqual(set(clustering.fractions), set(wf.fractions))
        
if __name__ == '__main__':
    unittest.main()
