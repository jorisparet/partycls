#!/usr/bin/env python

import unittest
import os

try:
    import atooms
    HAS_ATOOMS = True
except ModuleNotFoundError:
    HAS_ATOOMS = False

from partycls import Trajectory

class Test(unittest.TestCase):

    def setUp(self):
        self.data_dir = os.path.join(os.path.dirname(__file__), '../data/')

    @unittest.skipIf(not HAS_ATOOMS, 'no atooms module')
    def test_xyz(self):
        fname = os.path.join(self.data_dir, 'kalj_N150.xyz')
        traj = Trajectory(fname, fmt='xyz', backend='atooms')

        # Sanity checks
        self.assertEqual(len(traj), 101)
        self.assertEqual(len(traj[0].particle), 150)
        self.assertEqual(list(traj[0].distinct_species), ['A', 'B'])
        self.assertEqual(list(traj[0].cell.side), [5.0, 5.0, 5.0])

    @unittest.skipIf(not HAS_ATOOMS, 'no atooms module')
    def test_additional_fields(self):
        fname = os.path.join(self.data_dir, 'traj_with_neighbors.xyz')
        traj = Trajectory(fname, fmt='xyz', backend='atooms', additional_fields=['neighbors'])
        # Neighbors
        #  radius should not be read
        self.assertEqual(traj[0].particle[1].radius, 0.5,
                         "radius was not a requested additional field")
        #  neighbors should be read
        self.assertEqual(set(traj[0].particle[0].nearest_neighbors),
                         set([1,2,3]),
                         "neighbors were read incorrectly")
        # Read radii
        traj = Trajectory(fname, fmt='xyz', backend='atooms', additional_fields=['radius'])
        #  radius should be read
        self.assertEqual(traj[0].particle[1].radius, 0.4,
                         "radius was not a requested additional field")
        #  neighbors should not be read
        self.assertEqual(traj[0].particle[0].nearest_neighbors, None,
                         "neighbors was not a requested additional field")


    @unittest.skipIf(not HAS_ATOOMS, 'no atooms module')
    def test_rumd(self):
        fname = os.path.join(self.data_dir, 'kalj_N256_rho1.185_rumd.xyz.gz')
        traj = Trajectory(fname, fmt='rumd', backend='atooms')

        # Sanity checks
        self.assertEqual(len(traj), 1)
        self.assertEqual(len(traj[0].particle), 256)
        self.assertEqual(list(traj[0].distinct_species), ['0', '1'])
        self.assertEqual(list(traj[0].cell.side), [6.0, 6.0, 6.0])

    @unittest.skipIf(not HAS_ATOOMS, 'no atooms module')
    def test_superrumd(self):
        fname = os.path.join(self.data_dir, 'ka-2_N300_rumd/')
        traj = Trajectory(fname, fmt='superrumd', backend='atooms')

        # Sanity checks
        self.assertEqual(len(traj), 19)
        self.assertEqual(len(traj[0].particle), 303)
        self.assertEqual(list(traj[0].distinct_species), ['0', '1', '2'])
        self.assertEqual(list(traj[0].cell.side), [6.32053, 6.32053, 6.32053])
        
    @unittest.skipIf(not HAS_ATOOMS, 'no atooms module')
    def test_lammps(self):
        fname = os.path.join(self.data_dir, 'lj_N256_rho1.0.atom')
        traj = Trajectory(fname, fmt='lammps', backend='atooms')

        # Sanity checks
        self.assertEqual(len(traj), 6)
        self.assertEqual(len(traj[0].particle), 256)
        self.assertEqual(list(traj[0].distinct_species), ['1'])
        self.assertEqual(list(traj[0].cell.side), [6.3496, 6.3496, 6.3496])

    @unittest.skipIf(not HAS_ATOOMS, 'no atooms module')
    def test_folderlammps(self):
        fname = os.path.join(self.data_dir, 'lj_N256_rho1.0_lammps/')
        traj = Trajectory(fname, fmt='folderlammps', backend='atooms')

        # Sanity checks
        self.assertEqual(len(traj), 3)
        self.assertEqual(len(traj[0].particle), 256)
        self.assertEqual(list(traj[0].distinct_species), ['1'])
        self.assertEqual(list(traj[0].cell.side), [6.3496, 6.3496, 6.3496])
        
if __name__ == '__main__':
    unittest.main()
