import count_genotypes_sparse_constrained as gsc_contrained
import numpy as np
import unittest


class Test_simple_constrianed_genotype_count(unittest.TestCase):
    pos_array = np.array([1, 2, 3, 4, 5], dtype = np.int32)
    G_dict = {0: {1: {1}, 2: {2, 3}, -1:{4}},
                1: {1: {0}, 2: {1, 2, 3}, -1: {4}},
                2: {1: {1, 2, 3}, 2: {0}, -1 : {4}},
                3: {1: {1}, 2: {0, 2}, -1:{4}},
                4: {1: {1, 2}, 2: {0}, -1:{4}}}
    missing = True 
    n = len(G_dict)


    def test_none_filtered(self):
        threshold = 0
        Count = gsc_contrained.count_genotypes_sparse(self.G_dict, self.n, self.missing)
        Bools = gsc_contrained.count_genotypes_distance_constrained(self.pos_array,threshold)
        Count_filtered = Count[Bools]
        self.assertTrue(np.all(Count_filtered == Count))
        self.assertTrue(sum(Bools) == len(Count))

    def test_inclusive_filtered(self):
        threshold = 1
        Count = gsc_contrained.count_genotypes_sparse(self.G_dict, self.n, self.missing)
        Bools = gsc_contrained.count_genotypes_distance_constrained(self.pos_array,threshold)
        Count_filtered = Count[Bools]
        self.assertTrue(np.all(Count_filtered == Count[[1, 2, 3, 5, 6, 8], :]))

    def test_all_filtered(self):
        threshold = 10
        Count = gsc_contrained.count_genotypes_sparse(self.G_dict, self.n, self.missing)
        Bools = gsc_contrained.count_genotypes_distance_constrained(self.pos_array,threshold)
        Count_filtered = Count[Bools]
        self.assertTrue(np.all(Count_filtered == None))

if __name__ == '__main__':
    unittest.main()