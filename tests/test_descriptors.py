#!/usr/bin/env python

import unittest
import os

from pysc.trajectory import Trajectory
from pysc.descriptor import RadialDescriptor, AngularDescriptor
from pysc.descriptor import BondOrientationalDescriptor, LechnerDellagoDescriptor


class Test(unittest.TestCase):

    def setUp(self):
        data = os.path.join(os.path.dirname(__file__), '../data/')
        self.traj = Trajectory(os.path.join(data, 'kalj_N150.xyz'), first=0, last=10)
        self.cutoffs = [1.45, 1.25, 1.25, 1.075]

    def _compute(self, D):
        D.add_filter("species == 'B'", group=0)
        X = D.compute()
        
    def test_radial(self):
        D = RadialDescriptor(self.traj, rlim=(0.0, 2.5))
        self._compute(D)
        
    def test_angular(self):
        D = AngularDescriptor(self.traj)
        D.cutoffs = self.cutoffs
        self._compute(D)
        
    def test_steinhardt(self):
        D = BondOrientationalDescriptor(self.traj)
        D.cutoffs = self.cutoffs
        self._compute(D)
        
    def test_lechner_dellago(self):
        D = LechnerDellagoDescriptor(self.traj)
        D.cutoffs = self.cutoffs
        self._compute(D)


if __name__ == '__main__':
    unittest.main()
