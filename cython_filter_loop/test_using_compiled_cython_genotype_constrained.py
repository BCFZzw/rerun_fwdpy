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

    def test_example(self):
        pos_array = np.array([1, 2, 3], dtype = np.int32)
        G_dict = {0: {1: {0, 1}, 2: {2}},
                1: {1: {0}, 2: {1, 2}},
                2: {1: {1, 2}, 2: {0}}}
        missing = False 
        n = len(G_dict)
        threshold = 0
        Count = gsc_contrained.count_genotypes_sparse(G_dict, n, missing)
        print(Count)
        Bools = gsc_contrained.count_genotypes_distance_constrained(pos_array, threshold)
        Count_filtered = Count[Bools]
        print(Bools)
        print(Count_filtered)
        #self.assertTrue(np.all(Count_filtered == None))

    def test_homogzygous_polymorphic_sites_mixed(self):
        ### 2 homozygous sites, 2 individual, 1 pair
        genotype = allel.GenotypeArray([[[0, 0], [0, 1]],
        [[0, 0], [1, 1]],
        [[1, 0], [0, 1]]])
        G = genotype.to_n_alt(fill = -1)
        G_dict, missing = moments.LD.Parsing._sparsify_genotype_matrix(G)
        n = len(G_dict)
        Counts = gsc_contrained.count_genotypes_sparse(G_dict, n, missing)
        print(G)
        D2_pw, Dz_pw, pi2_pw, D_pw = moments.LD.Parsing.compute_pairwise_stats(G, genotypes = True)
        print(D2_pw)
        ### return [np.nan], length is the number of pairs
        #self.assertTrue(np.all(np.isnan(D2_pw)))

if __name__ == '__main__':
    unittest.main()