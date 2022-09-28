#!/usr/bin/env python

import unittest
import os

from partycls import Trajectory
from partycls.descriptor import RadialDescriptor, BondAngleDescriptor
from partycls.descriptor import BondOrientationalDescriptor, LechnerDellagoDescriptor
from partycls.descriptor import SmoothedBondOrientationalDescriptor, SmoothedBondAngleDescriptor
from partycls.descriptor import TetrahedralDescriptor
from partycls.descriptor import CompactnessDescriptor

from numpy import float32

class Test(unittest.TestCase):

    def setUp(self):
        self.data_dir = os.path.join(os.path.dirname(__file__), '../data/')
        self.traj = Trajectory(os.path.join(self.data_dir, 'kalj_N150.xyz'), first=0, last=10)
        self.traj.nearest_neighbors_method = 'fixed'
        self.traj.nearest_neighbors_cutoffs = [1.45, 1.25, 1.25, 1.075]

    def _compute(self, D):
        D.add_filter("species == 'B'", group=0)
        D.compute()
        
    def test_radial(self):
        D = RadialDescriptor(self.traj, dr=0.1)
        self._compute(D)
        # check automatically computed bounds
        self.assertEqual(D.bounds, (0.05, 2.45), 'wrong bounds for the radial grid')
        # check average value of g(r) at the first peak
        gr = D.normalize(D.average, method="gr")
        self.assertEqual(float32(gr[8]), float32(2.939288),
                         'wrong average value at the first peak of g(r)')

    def test_angular(self):
        D = BondAngleDescriptor(self.traj, dtheta=3.0)
        self._compute(D)
        q = D.normalize(D.average, method="pdf")
        self.assertEqual(float32(q[22]), float32(0.015544709),
                         'wrong average value at the peak \theta=67.5°')
        
    def test_steinhardt(self):
        D = BondOrientationalDescriptor(self.traj)
        self._compute(D)
        # test value of q_1
        self.assertEqual(float32(D.average[0]), float32(0.09393699),
                         'wrong average value for q_1')
        
    def test_lechner_dellago(self):
        D = LechnerDellagoDescriptor(self.traj)
        # D.cutoffs = self.cutoffs
        self._compute(D)
        self.assertEqual(float32(D.average[0]), float32(0.02164681),
                         'wrong average value for qbar_1')
        
    def test_smoothed_bo(self):
        # test average q_1
        D = SmoothedBondOrientationalDescriptor(self.traj, 
                                                cutoff_enlargement=1.3,
                                                exponent=8)
        self._compute(D)
        self.assertEqual(float32(D.average[0]), float32(0.06444802),
                         'wrong average value for qs_1')        
        # test convergence towards Steinhardt BO
        D = SmoothedBondOrientationalDescriptor(self.traj, 
                                                cutoff_enlargement=1.01,
                                                exponent=9223372036854775807)
        self._compute(D)
        self.assertAlmostEqual(float32(D.average[0]), float32(0.09393699),
                         places=3, msg='wrong average value for qs_1')
        self.assertAlmostEqual(float32(D.average[1]), float32(0.10234044),
                         places=3, msg='wrong average value for qs_2')
        self.assertAlmostEqual(float32(D.average[7]), float32(0.28741154),
                         places=3, msg='wrong average value for qs_7')
        
    
    def test_smoothed_ba(self):
        pass
    
    def test_radial_bo(self):
        pass
    
    def test_tetrahedral(self):
        D = TetrahedralDescriptor(self.traj)
        self._compute(D)
        # test average
        self.assertEqual(float32(D.average[0]), float32(0.48001548248253856),
                         'wrong average value for the tetrahedrality')

    def test_compactness(self):
        traj = Trajectory(os.path.join(self.data_dir, 'kalj_N150.xyz'), last=0)
        traj.compute_nearest_neighbors(method='voronoi')
        traj.set_property("radius", 0.54, subset="species == 'A'")
        traj.set_property("radius", 0.44, subset="species == 'B'")
        D = CompactnessDescriptor(traj)
        D.compute()
        self.assertEqual(float32(D.average[0]), float32(0.12505058),
                         'wrong average value for compactness')

if __name__ == '__main__':
    unittest.main()
