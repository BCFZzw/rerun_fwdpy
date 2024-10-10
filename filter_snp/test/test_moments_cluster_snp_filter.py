import unittest
import sys
sys.path.insert(1, '/home/alouette/projects/ctb-sgravel/alouette/Dz_sweep_final/simulation_review/filter_snp')
from cluster_snp_filter import *
import numpy as np
from moments import LD

class Test_moments_cluster_snp_filter(unittest.TestCase):

    def test_filtering_3_pos(self):
        L = 3
        n = 5
        np.random.seed(5553)
        G = np.random.randint(3, size=L * n).reshape(L, n)
        D2_pw, _, _, _ = LD.Parsing.compute_pairwise_stats(G, genotypes = True)
        pos_arr = [1, 401, 1400]
        ### only [1, 1400] remains
        filtered_pairs = list_pairs_within_threshold(pos_arr, 1000)
        filtered_index = bool_list_filtering(pos_arr, filtered_pairs)
        filtered_D2_pw = D2_pw[filtered_index]
        G_test = G[[0, 2], :]
        D2_check, _, _, _ = LD.Parsing.compute_pairwise_stats(G_test, genotypes = True)
        self.assertTrue(np.all(filtered_D2_pw == [D2_pw[1]]))
        self.assertTrue(np.all(filtered_D2_pw == D2_check))

if __name__ == '__main__':
    unittest.main()

