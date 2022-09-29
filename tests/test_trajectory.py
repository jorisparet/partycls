#!/usr/bin/env python

import unittest
import os

from partycls import Trajectory
from partycls.descriptor import BondAngleDescriptor
from partycls import Workflow, ZScore, PCA, KMeans

from numpy import float32

try:
    import atooms
    HAS_ATOOMS = True
except ModuleNotFoundError:
    HAS_ATOOMS = False

class Test(unittest.TestCase):

    def setUp(self):
        data = os.path.join(os.path.dirname(__file__), '../data/')
        self.traj = Trajectory(os.path.join(data, 'dislocation.xyz'), fmt='xyz')

    def test_xyz(self):
        data = os.path.join(os.path.dirname(__file__), '../data/')
        traj = Trajectory(os.path.join(data, 'dislocation.xyz'), fmt='xyz')
        self.assertEqual(len(traj[0].particle), 27)

    def _test_read_write(self, fmt, suffix, backend):
        # Trajectory adds the suffix by itself, which requires passing it explicitly
        import os
        import shutil
        
        tmp = '/tmp/partycls_tests'
        try:
            os.makedirs(tmp)
        except OSError:
            pass
        output_path = os.path.join(tmp, 'traj.{}'.format(suffix))
        tw = Trajectory(self.traj.filename)
        tw.write(output_path=output_path, fmt=fmt, backend=backend)
        tr = Trajectory(output_path, fmt=fmt, backend=backend)
        for pi, pj in zip(tr[0].particle, tw[0].particle):
            self.assertAlmostEqual(pi.position[0], pj.position[0], places=5)
        try:
            shutil.rmtree(tmp)
        except:
            pass
        
    def test_read_write_native(self):
        self._test_read_write(fmt='xyz', suffix='xyz', backend=None)
        self._test_read_write(fmt='rumd', suffix='xyz.gz', backend=None)

    @unittest.skipIf(not HAS_ATOOMS, 'no atooms module')
    def test_read_write_atooms(self):
        self._test_read_write(fmt='xyz', suffix='xyz', backend='atooms')
        # This must fail
        try:
            self._test_read_write(fmt='inexistent', suffix='xyz', backend='atooms')
            self.fail('this should have failed')
        except:
            pass
    
    # TODO: implement writing trajectories with MDTraj
    @unittest.skip('Not implemented yet')
    def test_read_write_xyz_mdtraj(self):
        # This fails because it requires a topology file. Really?!
        self._test_read_write(fmt='xyz', suffix='xyz', backend='mdtraj')
        # This must fail
        try:
            self._test_read_write(fmt='inexistent', suffix='xyz', backend='mdtraj')
            self.fail('this should have failed')
        except:
            pass
        
    def test_additional_list_field(self):
        data = os.path.join(os.path.dirname(__file__), '../data/')
        traj = Trajectory(os.path.join(data, 'traj_with_neighbors.xyz'), 
                          additional_fields=['radius', 'neighbors'])       
        neighbors_all = traj.dump('neighbors')
        self.assertEqual(list(neighbors_all[0][0]), [1,2,3])
        self.assertEqual(list(neighbors_all[1][1]), [2,3])
        neighbors_A = traj.dump('neighbors', subset="species == 'A'")
        self.assertEqual(list(neighbors_A[0][0]), [1,2,3])
        self.assertEqual(list(neighbors_A[1][0]), [1,2])
        
    def test_additional_field_label(self):
        data = os.path.join(os.path.dirname(__file__), '../data/')
        traj = Trajectory(os.path.join(data, 'traj_with_labels.xyz'),
                                       additional_fields=['cluster'])
        # check if 'cluster' was correctly set to Particle.label
        expected_labels = [1,0,1,1,1,0]
        count = 0
        for system in traj:
            for particle in system.particle:
                self.assertEqual(particle.label, expected_labels[count],
                                 'incorrect particle label')
                count += 1
            

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

    def test_compute_fixed_cutoffs(self):
        # test 
        data = os.path.join(os.path.dirname(__file__), '../data/')
        traj = Trajectory(os.path.join(data, 'wahn_N1000.xyz'))
        # compute cutoffs
        traj.compute_nearest_neighbors_cutoffs(dr=0.1)
        self.assertEqual(set(map(float32, traj.nearest_neighbors_cutoffs)),
                         set(map(float32, [1.45, 1.35, 1.35, 1.25])),
                         'wrong computed cutoffs')

    def test_voronoi_signatures(self):
        data = os.path.join(os.path.dirname(__file__), '../data/')
        traj = Trajectory(os.path.join(data, 'wahn_N1000.xyz'))
        # compute voronoi signatures
        traj.set_property('radius', 0.525, subset="species == 'A'")
        traj.set_property('radius', 0.475, subset="species == 'B'")
        traj.compute_voronoi_signatures()
        # check a few signatures
        signatures = traj[0].dump('signature')
        self.assertEqual(signatures[5], '0_0_12', 'wrong Voronoi signature')
        self.assertEqual(signatures[20], '0_2_8_2', 'wrong Voronoi signature')

    def test_show(self):
        data = os.path.join(os.path.dirname(__file__), '../data/')
        traj = Trajectory(os.path.join(data, 'kalj_N150.xyz'), last=3)
        traj.show(backend='matplotlib', outfile='traj_show_matplotlib')

    def test_angular_zscore_pca_kmeans(self):
        D = BondAngleDescriptor(self.traj)
        D.add_filter("species == 'A'", group=0)
        X = D.compute()
        scaler = ZScore()
        X = scaler.scale(X)        
        reducer = PCA(n_components=3)
        Xred = reducer.reduce(X)
        clustering = KMeans(n_clusters=2, n_init=100)
        clustering.fit(Xred)
        # print('Fractions :', clustering.fractions, '(clustering alone)')
        
        # Same via workflow
        wf = Workflow(self.traj, descriptor='ba', scaling='zscore',
                      dim_reduction='pca', clustering='kmeans')
        wf.descriptor.add_filter('species == "A"', group=0)
        wf.dim_reduction.n_components = 3
        wf.clustering.n_init = 100
        wf.disable_output()
        wf.run()
        # print('Fractions :', clustering.fractions, '(via workflow)')
        self.assertEqual(set(clustering.fractions), set(wf.fractions))
        
if __name__ == '__main__':
    unittest.main()
