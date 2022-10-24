#!/usr/bin/env python

import unittest
import os

from partycls import Trajectory
from partycls.descriptors import BondAngleDescriptor
from partycls import Workflow, ZScore, PCA, KMeans

class Test(unittest.TestCase):

    def setUp(self):
        data = os.path.join(os.path.dirname(__file__), '../data/')
        self.traj = Trajectory(os.path.join(data, 'dislocation.xyz'), fmt='xyz')

    def test_show_matplotlib(self):
        from partycls.helpers import show_matplotlib
        s = self.traj[0]
        fig = show_matplotlib(s, 'species')
        fig = show_matplotlib(s, 'radius')
        show_matplotlib(s, 'radius', outfile='matplotlib_radius')
        os.remove("matplotlib_radius.png")
        
    def test_show_ovito(self):
        from partycls.helpers import show_ovito
        s = self.traj[0]
        # species
        out = 'ovito_species.png'
        fig = show_ovito(s, 'species', outfile=out)
        os.remove(out)
        # radius
        out = 'ovito_radius.png'
        fig = show_ovito(s, 'radius', outfile=out)
        os.remove(out)

    def test_show_3dmol(self):
        from partycls.helpers import show_3dmol
        from random import randint
        s = self.traj[0]
        for p in s.particle:
            p.label = randint(0,1)
        view = show_3dmol(s, 'label')

    def test_sort_merge_clusters(self):
        # Adapted from 4th tutorial
        import random
        import numpy as np
        from partycls.descriptors import BondAngleDescriptor
        from partycls import PCA, GaussianMixture
        from partycls.helpers import merge_clusters, sort_clusters

        np.random.seed(1)
        random.seed(1)
        
        # compute descriptor
        D = BondAngleDescriptor(self.traj)
        D.cutoffs = [1.25]
        X = D.compute()

        # dimensionality reduction
        redux = PCA(n_components=2)
        X_red = redux.reduce(X)

        # Perform a clustering on the reduced feature space
        # of the previous example
        C = GaussianMixture(n_clusters=4, n_init=50)
        C.fit(X_red)
        labels = C.labels

        # Sort clusters
        new_labels, new_centroids = sort_clusters(C.labels, C.centroids(X_red))
        self.assertEqual(list(new_labels), [0, 1, 0, 0, 3, 0, 3, 1, 0, 1, 0, 2, 2,
                                            1, 2, 0, 1, 2, 1, 0, 2, 1, 1, 2, 1, 1, 2])
        
        # Use weights from GMM to merge the clusters into `n_cluster_min`
        # this returns new weights and new labels
        weights = C.backend.predict_proba(X_red)
        new_weights, new_labels = merge_clusters(weights, n_clusters_min=2)

        # TODO: check they make sense
        self.assertEqual(new_labels, [0, 0, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        
if __name__ == '__main__':
    unittest.main()
