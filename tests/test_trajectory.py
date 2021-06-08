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

    def test_get_property(self):
        data = os.path.join(os.path.dirname(__file__), '../data/')
        traj = Trajectory(os.path.join(data, 'traj_with_masses.xyz'), 
                          additional_fields=['mass'])
        
        # full dump
        dumps = {'pos[0]': [0.1, 0.0, 0.0],
                 'x': [0.1, 0.0, 0.0],
                 'mass': [0.5, 0.25, 0.25],
                 'particle.mass': [0.5, 0.25, 0.25],
                 'cell.side': [1.0, 1.0, 1.0]}
        for key, value in  dumps.items():
            self.assertEqual(value, list(traj[0].dump(key)))
            
        # subset
        dumps = {'pos[0]': [0.1],
                 'x': [0.1],
                 'mass': [0.5],
                 'particle.mass': [0.5]}
        for key, value in  dumps.items():
            self.assertEqual(value, list(traj[0].dump(key, "species == 'A'")))        
    
    def test_set_property(self):
        data = os.path.join(os.path.dirname(__file__), '../data/')
        traj = Trajectory(os.path.join(data, 'traj_with_masses.xyz'), 
                          additional_fields=['mass'])
        
        # cell
        traj.set_property('cell.side[0]', 2.0)
        self.assertEqual(traj[0].cell.side[0], 2.0)
        
        # particle (full)
        traj.set_property('radius', 0.4)
        self.assertEqual(list(traj[0].dump('radius')), [0.4, 0.4, 0.4])
        traj.set_property('radius', 0.7, "species == 'A'")
        self.assertEqual(list(traj[0].dump('radius', "species == 'A'")), [0.7])

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
