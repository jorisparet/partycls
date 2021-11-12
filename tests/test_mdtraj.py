#!/usr/bin/env python

import unittest
import os

# MDTraj
try:
    import mdtraj
    HAS_MDTRAJ = True
except ModuleNotFoundError:
    HAS_MDTRAJ = False

# h5py
try:
    import h5py
    import tables
    HAS_H5PY = True
except ModuleNotFoundError:
    HAS_H5PY = False

# networkx (HOOMD)
try:
    import networkx
    HAS_NETWORKX = True
except ModuleNotFoundError:
    HAS_NETWORKX = False
    
from partycls import Trajectory

class Test(unittest.TestCase):

    @unittest.skipIf(not HAS_MDTRAJ, 'no mdtraj module')
    def setUp(self):
        self.data_dir = os.path.join(os.path.dirname(__file__), '../data/')

    def test_gro(self):
        fname = os.path.join(self.data_dir, 'two_residues_same_resnum.gro')
        traj = Trajectory(fname, backend='mdtraj')

        # Sanity checks
        self.assertEqual(len(traj), 1)
        self.assertEqual(len(traj[0].particle), 193)
        self.assertEqual(list(traj[0].distinct_species), ['C', 'H', 'N', 'O'])
        self.assertEqual(list(traj[0].cell.side), [3.34809, 3.37095, 3.38426])

    @unittest.skipIf(not HAS_H5PY, 'no h5py module')
    def test_h5(self):
        fname = os.path.join(self.data_dir, 'frame0.h5')
        traj = Trajectory(fname, backend='mdtraj')

        # Sanity checks
        self.assertEqual(len(traj), 501)
        self.assertEqual(len(traj[0].particle), 22)
        self.assertEqual(set(traj[0].distinct_species), set(['C', 'H', 'N', 'O']))
        self.assertEqual(list(traj[0].cell.side), [1.0, 1.0, 1.0])

    def test_pdb(self):
        fname = os.path.join(self.data_dir, 'frame0.pdb')
        traj = Trajectory(fname, backend='mdtraj')

        # Sanity checks
        self.assertEqual(len(traj), 501)
        self.assertEqual(len(traj[0].particle), 22)
        self.assertEqual(set(traj[0].distinct_species), set(['C', 'H', 'N', 'O']))
        self.assertEqual(list(traj[0].cell.side), [1.0, 1.0, 1.0])

    def test_pdbgz(self):
        fname = os.path.join(self.data_dir, 'frame0.pdb.gz')
        traj = Trajectory(fname, backend='mdtraj')

        # Sanity checks
        self.assertEqual(len(traj), 501)
        self.assertEqual(len(traj[0].particle), 22)
        self.assertEqual(set(traj[0].distinct_species), set(['C', 'H', 'N', 'O']))
        self.assertEqual(list(traj[0].cell.side), [1.0, 1.0, 1.0])
        
    def test_arc(self):
        fname = os.path.join(self.data_dir, 'nitrogen.arc')
        traj = Trajectory(fname, backend='mdtraj')

        # Sanity checks
        self.assertEqual(len(traj), 50)
        self.assertEqual(len(traj[0].particle), 212)
        self.assertEqual(list(traj[0].distinct_species), ['N'])
        self.assertEqual(list(traj[0].cell.side), [1.8273600339889526]*3)

    @unittest.skipIf(not HAS_NETWORKX, 'no networkx module')
    def test_hoomdxml(self):
        fname = os.path.join(self.data_dir, 'water-box.hoomdxml')
        traj = Trajectory(fname, backend='mdtraj')

        # Sanity checks
        self.assertEqual(len(traj), 1)
        self.assertEqual(len(traj[0].particle), 302)
        self.assertEqual(list(traj[0].distinct_species), ['Cl', 'H', 'Na', 'O'])
        self.assertEqual(list(traj[0].cell.side), [52.62300129979849,
                                                   52.07199896220118,
                                                   51.33700057864189])
    
    # TODO: what happened between 1.9.6 and 1.9.7
    @unittest.skip('Stopped working with mdtraj-1.9.7')
    def test_xtc(self):
        fname = os.path.join(self.data_dir, 'frame0.xtc')
        top = os.path.join(self.data_dir, 'frame0.pdb')
        traj = Trajectory(fname, backend='mdtraj', top=top)

        # Sanity checks
        self.assertEqual(len(traj), 501)
        self.assertEqual(len(traj[0].particle), 22)
        self.assertEqual(set(traj[0].distinct_species), set(['C', 'H', 'N', 'O']))
        self.assertEqual(list(traj[0].cell.side), [2.573309898376465, 
                                                   2.573335647583008,
                                                   2.5733468532562256])
        
if __name__ == '__main__':
    unittest.main()
