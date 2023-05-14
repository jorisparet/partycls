#!/usr/bin/env python

from tabnanny import verbose
import unittest
import os

from partycls import Trajectory
from partycls.descriptors.descriptor import StructuralDescriptor
from partycls.descriptors import RadialDescriptor, BondAngleDescriptor
from partycls.descriptors import BondOrientationalDescriptor, LocallyAveragedBondOrientationalDescriptor
from partycls.descriptors import SmoothedBondOrientationalDescriptor, SmoothedBondAngleDescriptor
from partycls.descriptors import RadialBondOrientationalDescriptor
from partycls.descriptors import TetrahedralDescriptor
from partycls.descriptors import CompactnessDescriptor
from partycls.descriptors import CoordinationDescriptor

import numpy
from numpy import float32

try:
    import pyvoro
    HAS_PYVORO = True
except ModuleNotFoundError:
    HAS_PYVORO = False

try:
    from partycls.descriptors import DscribeChemicalDescriptor
    from dscribe.descriptors import SOAP
    HAS_DSCRIBE = True
except ModuleNotFoundError:
    HAS_DSCRIBE = False

class Test(unittest.TestCase):

    def setUp(self):
        self.data_dir = os.path.join(os.path.dirname(__file__), '../data/')
        self.traj = Trajectory(os.path.join(self.data_dir, 'kalj_N150.xyz'), first=0, last=10)
        self.traj.nearest_neighbors_method = 'fixed'
        self.traj.nearest_neighbors_cutoffs = [1.45, 1.25, 1.25, 1.075]
        self.traj.compute_nearest_neighbors()

    def _compute(self, D):
        D.add_filter("species == 'B'", group=0)
        D.compute()
        
    def test_filtered_neighbors(self):
        traj = Trajectory(os.path.join(self.data_dir, 'kalj_N150.xyz'), last=0)
        traj.nearest_neighbors_method = 'fixed'
        traj.nearest_neighbors_cutoffs = [1.45, 1.25, 1.25, 1.075]
        traj.compute_nearest_neighbors()
        D = StructuralDescriptor(self.traj)
        # filters
        filter_g0 = "species == 'B'"
        filter_g1 = "species == 'A'"
        D.add_filter(filter_g0, group=0)
        D.add_filter(filter_g1, group=1)
        # filter neighbors that are not in group=1 in the descriptor
        D._filter_neighbors()
        # compare original neighbors and filtered neighbors
        original_neighbors = traj[0].dump('neighbors', subset=filter_g0)
        filtered_neighbors = D._neighbors[0]
        neigh_numbers = D._neighbors_number[0]
        for i in range(len(D.groups[0][0])):
            nn_i = neigh_numbers[i]
            original_set = set(original_neighbors[i])
            filtered_set = set(filtered_neighbors[i][0:nn_i])
            # neighbors that were filtered out in `filtered_neighbors`
            filtered_out = filtered_set ^ original_set
            # are the two sets identical if we include filtered-out neighbors again?
            self.assertEqual((filtered_set | filtered_out), original_set,
                             'sets are different')
            # are filtered-out particles of type 'B'?
            for j in filtered_out:
                self.assertEqual(traj[0].particle[j].species, 'B',
                                 'particle was not filtered properly')

    def test_filtered_subsidiary_neighbors(self):
        traj = Trajectory(os.path.join(self.data_dir, 'kalj_N150.xyz'), last=0)
        traj.nearest_neighbors_method = 'fixed'
        traj.nearest_neighbors_cutoffs = [1.45, 1.25, 1.25, 1.075]
        traj.compute_nearest_neighbors()
        # no filters
        D = StructuralDescriptor(traj)
        D.add_filter("species == 'A'", group=0)
        D.add_filter("species == 'B'", group=1)
        D._filter_neighbors()
        D._filter_subsidiary_neighbors()
        # neighbors and filtered neighbors
        i = 12
        ni = traj[0].particle[i].nearest_neighbors
        nn_i_filtered = D._neighbors_number[0][i]
        ni_filtered = D._neighbors[0][i][0:nn_i_filtered]
        # neighbors of neighbors (with and without filters)
        n_ni = []
        for j in ni:
            nj = traj[0].particle[j].nearest_neighbors
            n_ni.append(nj)
        n_ni_filtered = D._subsidiary_neighbors[0][i]
        # test filtered neighbors of i
        self.assertEqual(set(ni) & set(ni_filtered), {142, 145, 148},
                         'wrong filtered neighbors')
        # neighbors of neighbors (j=145,148,142)
        self.assertEqual(set(traj[0].particle[ni_filtered[0]].nearest_neighbors) & set(n_ni_filtered[0]),
                         {142}, 'wrong neighbors of neighbors')
        self.assertEqual(set(traj[0].particle[ni_filtered[1]].nearest_neighbors) & set(n_ni_filtered[1]),
                         set(), 'wrong neighbors of neighbors')
        self.assertEqual(set(traj[0].particle[ni_filtered[2]].nearest_neighbors) & set(n_ni_filtered[2]),
                         {145}, 'wrong neighbors of neighbors')

    def test_extended_neighbors(self):
        # trajectory
        traj = Trajectory(os.path.join(self.data_dir, 'SiO2_N2000.xyz'))
        traj.nearest_neighbors_method = 'fixed'
        cutoffs = [1.1800, 0.6300, 0.6300, 1.0200]
        traj.nearest_neighbors_cutoffs = cutoffs
        traj.compute_nearest_neighbors()
        # full descriptor (standard cutoffs)
        D = StructuralDescriptor(traj)
        D._compute_extended_neighbors(cutoffs)
        indices = D.dump('_index', group=0)
        for i_, i in enumerate(indices[0]):
            ni = set(traj[0].particle[i].nearest_neighbors)
            ni_ex = ni_ex = set(D._extended_neighbors[0][i_][0:D._extended_neighbors_number[0][i_]])
            # neighbors must be the same
            self.assertEqual(ni, ni_ex, 'neighbors are different')
        # full descriptor (extended cutoffs)
        extended_cutoffs = 1.3 * numpy.array(cutoffs)
        D._compute_extended_neighbors(extended_cutoffs)
        ni = set(traj[0].particle[0].nearest_neighbors)
        ni_ex = set(D._extended_neighbors[0][0][0:D._extended_neighbors_number[0][0]])
        self.assertEqual(ni ^ ni_ex, {336, 341}, 'wrong extended neighbors')
        # partial descriptor 1-1 (standard cutoffs)
        D.add_filter("species == '1'", group=0)
        D.add_filter("species == '1'", group=1)
        D._compute_extended_neighbors(cutoffs)
        ni_full = set(traj[0].particle[0].nearest_neighbors)
        ni_11_ex = set(D._extended_neighbors[0][0][0:D._extended_neighbors_number[0][0]])
        ni_11_filtered = ni_full ^ ni_11_ex
        self.assertEqual(ni_11_filtered, {670, 743, 957, 1803},
                         'filtered neighbors are different')
        # partial descriptor 1-2 (standard cutoffs)
        D.clear_filters(group=1)
        D.add_filter("species == '2'", group=1)
        D._compute_extended_neighbors(cutoffs)
        ni_full = set(traj[0].particle[0].nearest_neighbors)
        ni_12_ex = set(D._extended_neighbors[0][0][0:D._extended_neighbors_number[0][0]])
        ni_12_filtered = ni_full ^ ni_12_ex
        self.assertEqual(ni_12_filtered, {209, 345, 513, 625},
                         'filtered neighbors are different')
        # sum of filtered neighbors of both partial descriptors
        # must be equal to the original neighbors
        self.assertEqual(ni_11_ex | ni_12_ex, ni_full)

    def test_radial(self):
        D = RadialDescriptor(self.traj, dr=0.1)
        self._compute(D)
        # check automatically computed bounds
        self.assertEqual(D.bounds, (0.05, 2.45),
                         'wrong bounds for the radial grid')
        # check average value of g(r) at the first peak
        gr = D.normalize(D.average, method="gr")
        self.assertEqual(float32(gr[8]), float32(2.939288),
                         'wrong average value at the first peak of g(r)')
        # check g(r) is still correct with non-trivial bounds
        # (max. is at bin=3 instead of bin=8 when rmin=0.5)
        D_shifted = RadialDescriptor(self.traj, bounds=(0.5, 2.5), dr=0.1)
        self._compute(D_shifted)
        g_shifted = D_shifted.normalize(D_shifted.average, method="gr")
        self.assertEqual(float32(gr[8]), float32(g_shifted[3]),
                         'the two g(r) are not equal at the first peak')
        # check dr changes the grid dynamically
        len_old = len(D.grid)
        D.dr = 0.05
        len_new = len(D.grid)
        self.assertGreater(len_new, len_old, 'incorrect grid size')

    def test_angular(self):
        # Average distribution
        D = BondAngleDescriptor(self.traj, dtheta=3.0)
        self._compute(D)
        q = D.normalize(D.average, method="pdf")
        self.assertEqual(float32(q[22]), float32(0.015544709),
                         'wrong average value at the peak theta=67.5°')
        
        # Partial distributions
        #  trajectory
        traj = Trajectory(os.path.join(self.data_dir, 'wahn_N1000.xyz'))
        traj.nearest_neighbors_method = 'fixed'
        traj.nearest_neighbors_cutoffs = [1.4250, 1.3250, 1.3250, 1.2750] #[1.45, 1.175, 1.175, 1.05]
        traj.compute_nearest_neighbors()
        #  distributions
        #   B-all
        D_B = BondAngleDescriptor(traj)
        D_B.add_filter("species == 'B'", group=0)
        _ = D_B.compute()
        #   # B-A
        D_BA = BondAngleDescriptor(traj)
        D_BA.add_filter("species == 'B'", group=0)
        D_BA.add_filter("species == 'A'", group=1)
        _ = D_BA.compute()
        #   B-B
        D_BB = BondAngleDescriptor(traj)
        D_BB.add_filter("species == 'B'", group=0)
        D_BB.add_filter("species == 'B'", group=1)
        _ = D_BB.compute()
        #  check error on the average of the partials
        x_A, x_B = traj[0].chemical_fractions
        q_B = D_B.normalize(D_B.average, method='pdf')
        q_part = D_B.normalize(x_A*D_BA.average + x_B*D_BB.average, method='pdf')
        mean_avg_err = numpy.abs(q_part - q_B).mean()
        self.assertLess(mean_avg_err, 1e-3, 'average error is too large')
        
    def test_steinhardt(self):
        D = BondOrientationalDescriptor(self.traj, lmin=2, lmax=4)
        # grid
        self.assertEqual(set(D.grid), set([2,3,4]), 'wrong grid')
        D.orders = [1,2,3,4,5,6,7,8]
        self.assertEqual(set(D.grid), set([1,2,3,4,5,6,7,8]), 'wrong grid')
        # test value of q_1
        self._compute(D)
        self.assertEqual(float32(D.average[0]), float32(0.09393699),
                         'wrong average value for q_1')
        
    def test_lechner_dellago(self):
        D = LocallyAveragedBondOrientationalDescriptor(self.traj)
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
        self.assertEqual(float32(D.average[0]), float32(0.052115675),
                         'wrong average value for qs_1')        
        # test convergence towards Steinhardt BO
        D = SmoothedBondOrientationalDescriptor(self.traj, 
                                                cutoff_enlargement=1.3,
                                                exponent=1000000000)
        self._compute(D)
        self.assertAlmostEqual(float32(D.average[0]), float32(0.09393699),
                         places=8, msg='wrong average value for qs_1 (not converged)')
        self.assertAlmostEqual(float32(D.average[1]), float32(0.10234044),
                         places=8, msg='wrong average value for qs_2 (not converged)')
        self.assertAlmostEqual(float32(D.average[7]), float32(0.28741154),
                         places=8, msg='wrong average value for qs_7 (not converged)')
        
    def test_smoothed_ba(self):
        # compute the descriptor on particles A for better statistics
        D = SmoothedBondAngleDescriptor(self.traj)
        D.add_filter("species == 'A'", group=0)
        D.compute()
        # value at the peak
        q = D.normalize(D.average, method="pdf")
        self.assertEqual(float32(q[18]), float32(0.01290775), 
                         'wrong average value for first peak')
        # convergence towards non-smoothed descriptor
        #  sba
        D_sba = SmoothedBondAngleDescriptor(self.traj, dtheta=3.0, 
                                            cutoff_enlargement=1.3,
                                            exponent=1000000000)
        D_sba.add_filter("species == 'A'", group=0)
        D_sba.compute()
        q_sba = D_sba.normalize(D_sba.average, method="pdf")
        #  ba
        D_ba = BondAngleDescriptor(self.traj)
        D_ba.add_filter("species == 'A'", group=0)
        D_ba.compute()
        q_ba = D.normalize(D_ba.average, method="pdf")
        #  compare
        self.assertAlmostEqual(max(q_sba), max(q_ba), places=8,
                               msg='wrong average value for first peak (not converged)')
    
    def test_radial_bo(self):
        D = RadialBondOrientationalDescriptor(self.traj, bounds=(1.1, 1.5), dr=0.1)
        # bounds
        self.assertEqual(set(map(float32, D.distance_grid)),
                         set(map(float32, [1.1, 1.2, 1.3, 1.4, 1.5])),
                         'incorrect bounds')
        D.dr = 0.2
        self.assertEqual(D.dr, 0.2, 'wrong value for dr')
        D.bounds = (1.1, 1.9)
        self.assertEqual(set(map(float32, D.distance_grid)),
                         set(map(float32, [1.1, 1.3, 1.5, 1.7, 1.9])),
                         'incorrect bounds')
        # distance grid
        r_grid = [1.1, 1.2, 1.3]
        D.distance_grid = r_grid
        self.assertEqual(set( map(float32, D.distance_grid) ),
                         set( map(float32, r_grid)),
                         'incorrect distance grid')
        self.assertEqual(D.bounds, (1.1, 1.3), 
                         'incorrect bounds associated to distance grid')
        # grid
        self.assertEqual(D.grid[0], (1, 1.1),
                         'incorrect (l,r) grid')
        self.assertEqual(D.grid[-1], (8, 1.3),
                         'incorrect (l,r) grid')
        # test average
        self._compute(D)
        self.assertEqual(float32(D.average[0]), float32(0.27672228),
                         'wrong average value for tuple (l,r)={}'.format(D.grid[0]))
        self.assertEqual(float32(D.average[-1]), float32(0.43605693),
                         'wrong average value for tuple (l,r)={}'.format(D.grid[-1]))
    
    def test_tetrahedral(self):
        D = TetrahedralDescriptor(self.traj)
        self._compute(D)
        # test average
        self.assertEqual(float32(D.average[0]), float32(0.48001548248253856),
                         'wrong average value for the tetrahedrality')

    @unittest.skipIf(not HAS_PYVORO, 'no pyvoro module')
    def test_compactness(self):
        # use only one frame
        traj = Trajectory(os.path.join(self.data_dir, 'wahn_N1000.xyz'))
        traj.compute_nearest_neighbors(method='voronoi')
        # radii based on the first peak of g_aa(r)
        traj.set_property("radius", 0.54, subset="species == 'A'")
        traj.set_property("radius", 0.43, subset="species == 'B'")
        D = CompactnessDescriptor(traj)
        D.compute()
        self.assertEqual(float32(D.average[0]), float32(0.11607557),
                         'wrong average value for compactness')

    def test_coordination(self):
        # trajectory and neighbors
        traj = Trajectory(os.path.join(self.data_dir, 'kalj_N150.xyz'), last=0)
        traj.nearest_neighbors_method = 'fixed'
        traj.nearest_neighbors_cutoffs = [1.45, 1.25, 1.25, 1.075]
        traj.compute_nearest_neighbors()
        # total (filter on group=0)
        D = CoordinationDescriptor(traj, total=True, partial=False)
        D.add_filter("species == 'B'", group=0)
        X_tot = D.compute()
        for i, pi in enumerate(D.groups[0][0]):
            self.assertEqual(len(pi.nearest_neighbors), X_tot[i,0], 
                             'wrong total coordination number (g0)')
        # total (filter on group=0)
        D.add_filter("species == 'B'", group=1)
        X_tot = D.compute()
        n_B = [1, 0, 1, 0, 1, 0, 0, 0, 0, 0]
        for i in range(len(n_B)):
            self.assertEqual(n_B[i], X_tot[i,0], 
                            "wrong total coordination number (g0,g1)")
        # total + partial (filter on group=0)
        D = CoordinationDescriptor(traj, total=True, partial=True)
        D.add_filter("species == 'B'", group=0)
        self.assertEqual(set(D.grid), set(['all', 'A', 'B']),
                         "wrong grid of species (g0)")
        X_all = D.compute()
        for i in range(len(D.groups[0][0])):
            self.assertEqual(X_all[i,0], sum(X_all[i,1:]))
        # change grid
        D.total = False
        self.assertEqual(set(D.grid), set(['A', 'B']))
        D.total = True
        D.partial = False
        self.assertEqual(set(D.grid), set(['all']))

    #@unittest.skipIf(not HAS_DSCRIBE, 'no dscribe module')
    @unittest.skipIf(True, 'deprecation error with numpy==1.24')
    def test_dscribe(self):
        # trajectory and neighbors
        traj = Trajectory(os.path.join(self.data_dir, 'wahn_N1000.xyz'))
        traj.nearest_neighbors_method = 'fixed'
        traj.nearest_neighbors_cutoffs = [1.4250, 1.3250, 1.3250, 1.2750]
        traj.compute_nearest_neighbors()
        # set atomic symbols to species
        traj.set_property('species', 'C', "species == 'A'")
        traj.set_property('species', 'H', "species == 'B'")
        # descriptor
        D = DscribeChemicalDescriptor(traj, backend=SOAP, sigma=0.1, rcut=3.0,
                                      lmax=7, nmax=7, rbf='gto', verbose=True)
        X = D.compute()

    def test_handle_nans(self):
        # trajectory and neighbors
        traj = Trajectory(os.path.join(self.data_dir, 'wahn_N1000.xyz'))
        traj.nearest_neighbors_method = 'fixed'
        traj.nearest_neighbors_cutoffs = [1.4250, 1.3250, 1.3250, 1.2750]
        traj.compute_nearest_neighbors()
        # arbitrary descriptor
        D = BondOrientationalDescriptor(traj)
        D.verbose = True
        D.accept_nans = False
        D.add_filter("species == 'B'", group=0)
        D.add_filter("species == 'B'", group=1)
        X = D.compute()


if __name__ == '__main__':
    unittest.main()
