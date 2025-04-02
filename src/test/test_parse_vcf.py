import unittest
import zarr
import sys
sys.path.insert(1, '/home/alouette/projects/ctb-sgravel/alouette/Dz_sweep_final/simulation_review/src')
from parse_vcf import *
import pandas as pd
import moments.LD


class Test_parse_vcf(unittest.TestCase):

    panel_file = "/home/alouette/projects/ctb-sgravel/alouette/1000Genome/population_panel/integrated_call_samples_v3.20130502.2504.ALL.panel"
    zarr_path = "/home/alouette/projects/ctb-sgravel/data/30x1000G_biallelic_strict_masked/zarrFormat/chr22"
    callset = zarr.open_group(zarr_path, mode='r')
    pos_array = allel.SortedIndex(callset['variants/POS'])
    callset_samples = callset["samples"] # to view, use [:], the code for individual is returned
    panel_df = pd.read_csv(panel_file, sep = "\t")

    def test_get_pops(self):
        populations = get_pops_from_superpop(self.panel_file, super_pop = "AFR")
        populations.sort()
        self.assertTrue(np.all(populations == ['ACB', 'ASW', 'ESN', 'GWD', 'LWK', 'MSL', 'YRI']))

    def test_catch_pop_superpop_both_specified(self):
        with self.assertRaises(ValueError):
            locate_panel_individuals(self.callset_samples, self.panel_file, pop = "GBR", super_pop = "EUR")

    def test_catch_pop_not_exist(self):
        with self.assertRaises(ValueError):
            locate_panel_individuals(self.callset_samples, self.panel_file, pop = "ABE")
            
    def test_subset_AFR(self):
        loc_samples = locate_panel_individuals(self.callset_samples, self.panel_file, super_pop = "AFR")
        self.assertTrue(len(loc_samples) == sum(self.panel_df.super_pop == "AFR"))
        self.assertTrue(len(loc_samples) == 661)
        samples = self.panel_df[self.panel_df["sample"].isin(self.callset["samples"][loc_samples])]
        self.assertTrue(samples["super_pop"].unique() == ['AFR'])
        for ind in self.callset["samples"][loc_samples]:
            self.assertTrue(self.panel_df[self.panel_df["sample"] == ind]["super_pop"].tolist()[0] == "AFR")

    def test_jackknife_YRI(self):
        loc_samples = locate_jackknife_individuals(self.callset_samples, self.panel_file, super_pop = "AFR", jackknife_pop = "YRI")
        self.assertTrue(len(loc_samples) == sum((self.panel_df.super_pop == "AFR") & (self.panel_df["pop"] != "YRI"))) # 553
        samples = self.panel_df[self.panel_df["sample"].isin(self.callset["samples"][loc_samples])]
        self.assertTrue(np.all(np.sort(samples["pop"].unique()) == ['ACB', 'ASW', 'ESN', 'GWD', 'LWK', 'MSL']))
        self.assertTrue(np.all(samples["super_pop"].unique() == ['AFR']))
        for ind in self.callset["samples"][loc_samples]:
            self.assertTrue(self.panel_df[self.panel_df["sample"] == ind]["pop"].tolist()[0] != "YRI")

    def test_catch_jackknife_not_exist(self):
        with self.assertRaises(ValueError):
            locate_jackknife_individuals(self.callset_samples, self.panel_file, super_pop = "AFR", jackknife_pop = "CHS")

    def test_all_panel(self):
        loc_samples = locate_panel_individuals(self.callset_samples, self.panel_file)
        self.assertTrue(len(loc_samples) == 2504)

    def test_no_snp_in_range(self):
        loc_region = locate_genotype_region(self.pos_array, 0, 1)
        self.assertTrue(loc_region == slice(None, 0, None))

    def test_no_snp_in_range_exclusive(self):
        """
        From the window side, it needs to be non-overlapping, exlusive on the right hand side
        """
        pos_arr = allel.SortedIndex([10])
        loc_region = locate_genotype_region(pos_arr, pos_end = 10)
        self.assertTrue(loc_region == slice(None, 0, None))
        pos_arr = allel.SortedIndex([10, 11])
        loc_region = locate_genotype_region(pos_arr, pos_end = 11)
        self.assertTrue(loc_region == slice(0, 2, None))
        self.assertTrue(np.all(pos_arr[loc_region] == pos_arr))
        pos_arr = allel.SortedIndex([1, 10])
        loc_region = locate_genotype_region(pos_arr, 1, 10)
        self.assertTrue(loc_region == slice(None, 0, None))
    

    def test_take_all_snps(self):
        loc_region = locate_genotype_region(self.pos_array)
        self.assertTrue(np.all(self.pos_array == self.pos_array[loc_region]))

    def test_no_snp_moment(self):
        genotype_012_dask, pos_array = get_genotype012(self.zarr_path, pos_start = 0, pos_end = 1, panel_file = self.panel_file)
        self.assertTrue(not genotype_012_dask.any())
        self.assertTrue(np.all(pos_array == []))
        ### this will return [] for all 
        D2_pw, Dz_pw, pi2_pw, D_pw = moments.LD.Parsing.compute_pairwise_stats(genotype_012_dask, genotypes = True)
        self.assertTrue(np.all(D2_pw == []))
        self.assertTrue(D2_pw.size == 0)

    def test_catch_reverse_order(self):
        with self.assertRaises(AssertionError):
            locate_genotype_region(self.pos_array, 100000000, 1)

    def test_negative_floating_value_works(self):
        loc_region = locate_genotype_region(self.pos_array, -1, 100000000.1)
        self.assertTrue(max(self.pos_array[loc_region]) <= 100000000.1)

    def test_maximum_value_all_region(self):
        loc_region = locate_genotype_region(self.pos_array, pos_start = None, pos_end = max(self.pos_array))
        self.assertTrue(np.all(self.pos_array[loc_region] == self.pos_array))

    def test_min_value_all_region(self):
        loc_region = locate_genotype_region(self.pos_array, pos_start = 0, pos_end = None)
        self.assertTrue(np.all(self.pos_array[loc_region] == self.pos_array))

    def test_none_all_region(self):
        loc_region = locate_genotype_region(self.pos_array, pos_start = None, pos_end = None)
        self.assertTrue(np.all(self.pos_array[loc_region] == self.pos_array))

    def test_inclusive_region(self):
        loc_region = locate_genotype_region(self.pos_array, pos_start = None, pos_end = 10666327)
        self.assertTrue(len(range(loc_region.stop)[loc_region]) == 2)

    def test_genotype_parsing(self):
        genotype_012_dask, pos_array = get_genotype012(self.zarr_path, pos_start = None, pos_end = 10666327, panel_file = self.panel_file)
        self.assertTrue(len(pos_array) == 2)
        loc_samples = locate_panel_individuals(self.callset_samples, self.panel_file)
        test_genotype = allel.GenotypeArray(self.callset['calldata/GT'][:2])
        test_genotype = test_genotype.take(loc_samples, axis = 1)
        test_genotype_012 = test_genotype.to_n_alt(fill=-1)
        self.assertTrue(np.all(test_genotype_012 == genotype_012_dask))

    def test_genotype_parsing_jackknife(self):
        genotype_012_dask, pos_array = get_genotype012(self.zarr_path, pos_start = None, pos_end = 10666327, panel_file = self.panel_file, super_pop = "AFR", jackknife_pop = "YRI")
        self.assertTrue(len(pos_array) == 2)
        loc_samples = locate_jackknife_individuals(self.callset_samples, self.panel_file, super_pop = "AFR", jackknife_pop = "YRI")
        test_genotype = allel.GenotypeArray(self.callset['calldata/GT'][:2])
        test_genotype = test_genotype.take(loc_samples, axis = 1)
        test_genotype_012 = test_genotype.to_n_alt(fill=-1)
        self.assertTrue(np.all(test_genotype_012 == genotype_012_dask))

    def test_connecting_with_moments(self):
        genotype_012, pos_array = get_genotype012(self.zarr_path, pos_start = None, pos_end = 10666327, panel_file = self.panel_file)
        D2_pw, Dz_pw, pi2_pw, D_pw = moments.LD.Parsing.compute_pairwise_stats(genotype_012, genotypes = True)
    
    def test_connecting_with_moments_constrained(self):
        threshold = 1
        genotype_012, pos_array = get_genotype012(self.zarr_path, pos_start = None, pos_end = 10666327, panel_file = self.panel_file)
        D2_pw_subset, Dz_pw_subset, pi2_pw_subset, D_pw_subset = moments.LD.Parsing.compute_pairwise_stats(genotype_012, pos_array, genotypes = True, distance_constrained = threshold)

if __name__ == '__main__':
    unittest.main()

            