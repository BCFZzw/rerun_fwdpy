import numpy as np
import unittest
import moments.LD


class Test_simple_constrianed_genotype_count(unittest.TestCase):
    L = 5
    n = 5
    np.random.seed(5553)
    G = np.random.randint(3, size=L * n).reshape(L, n)
    # working as the previous unmodified version
    D2_pw, Dz_pw, pi2_pw, D_pw = moments.LD.Parsing.compute_pairwise_stats(G, genotypes = True)

    def test_catch_pos_array_dtype(self):
        with self.assertRaises(AssertionError):
            pos_array = np.array([1, 2, 3, 4, 5], dtype = np.int64)
            threshold = 2
            moments.LD.Parsing.compute_pairwise_stats(self.G, pos_array, genotypes = True, distance_constrained = threshold)

            
if __name__ == '__main__':
    unittest.main()