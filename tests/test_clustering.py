#!/usr/bin/env python

import unittest
import os

from partycls import Trajectory
from partycls.descriptor import BondAngleDescriptor, RadialDescriptor
from partycls import Workflow, ZScore, PCA, KMeans, CommunityInference

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
        
        # Same via workflow
        wf = Workflow(self.traj, descriptor='ba', scaling='zscore', clustering='kmeans')
        wf.descriptor.cutoffs = self.cutoffs
        wf.clustering.n_init = 100
        wf.disable_output()
        wf.run()
        # check if both methods give the same result
        self.assertEqual(set(wf.fractions), set(clustering.fractions),
                         'different cluster fractions')

    def test_radial_ci(self):
        D = RadialDescriptor(self.traj)
        D.compute()
        clustering = CommunityInference(n_init=100)
        clustering.fit(D)
        # check if the dislocated particles are in the same cluster
        labels = list(clustering.labels)
        self.assertEqual(labels[4], labels[6], 'not the expected labels')
        # check if only the dislocated particles are in the same cluster
        self.assertEqual(set(clustering.fractions), set([2/27, 25/27]),
                         'not the expected cluster fractions')
        
        # Same via workflow
        wf = Workflow(self.traj, descriptor='gr', clustering='cinf')
        wf.clustering.n_init = 100
        wf.disable_output()
        wf.run()
        # check if both methods give the same result
        self.assertEqual(set(wf.fractions), set(clustering.fractions),
                         'different cluster fractions')

        # Dummy test of outputs methods
        tmp = '/tmp/partycls_tests'
        try:
            os.makedirs(tmp)
        except OSError:
            pass
        wf.write_log(filename=os.path.join(tmp, 'traj'))
        wf.write_trajectory(filename=os.path.join(tmp, 'traj'))
        wf.write_centroids(filename=os.path.join(tmp, 'traj'))
        try:
            shutil.rmtree(tmp)
        except:
            pass
 
if __name__ == '__main__':
    unittest.main()
