#!/usr/bin/env python

import unittest
import os

from partycls import Trajectory
from partycls.descriptor import RadialDescriptor, BondAngleDescriptor
from partycls.descriptor import BondOrientationalDescriptor, LechnerDellagoDescriptor

from numpy import float32

class Test(unittest.TestCase):

    def setUp(self):
        data = os.path.join(os.path.dirname(__file__), '../data/')
        self.traj = Trajectory(os.path.join(data, 'kalj_N150.xyz'), first=0, last=10)
        self.cutoffs = [1.45, 1.25, 1.25, 1.075]

    def _compute(self, D):
        D.add_filter("species == 'B'", group=0)
        D.compute()
        
    def test_radial(self):
        D = RadialDescriptor(self.traj, dr=0.1)
        self._compute(D)
        # check automatically computed bounds
        self.assertEqual(D.bounds, (0.05, 2.45), 'wrong bounds for the radial grid')
        # check average value of g(r) at the first peak
        gr = D.normalize_gr(D.average)
        self.assertEqual(float32(gr[8]), float32(2.939288),
                         'wrong average value at the first peak of g(r)')

    def test_angular(self):
        D = BondAngleDescriptor(self.traj, dtheta=3.0)
        D.cutoffs = self.cutoffs
        self._compute(D)
        q = D.normalize_sin(D.average)
        self.assertEqual(float32(q[22]), float32(0.015544709),
                         'wrong average value at the peak \theta=67.5°')
        
    def test_steinhardt(self):
        D = BondOrientationalDescriptor(self.traj)
        D.cutoffs = self.cutoffs
        self._compute(D)
        # test value of q_1
        self.assertEqual(float32(D.average[0]), float32(0.09393699),
                         'wrong average value for q_1')
        
    def test_lechner_dellago(self):
        D = LechnerDellagoDescriptor(self.traj)
        D.cutoffs = self.cutoffs
        self._compute(D)
        self.assertEqual(float32(D.average[0]), float32(0.02164681),
                         'wrong average value for qbar_1')

if __name__ == '__main__':
    unittest.main()
