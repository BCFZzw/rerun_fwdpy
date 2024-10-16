import unittest
import zarr
import sys
sys.path.insert(1, '/home/alouette/projects/ctb-sgravel/alouette/Dz_sweep_final/simulation_review/src')
from parse_vcf import *
import pandas as pd


class Test_parse_vcf(unittest.TestCase):

    panel_file = "/home/alouette/projects/ctb-sgravel/alouette/1000Genome/population_panel/integrated_call_samples_v3.20130502.2504.ALL.panel"
    zarr_path = "/home/alouette/projects/ctb-sgravel/data/30x1000G_biallelic_strict_masked/zarrFormat/chr22"
    callset = zarr.open_group(zarr_path, mode='r')
    pos_array = allel.SortedIndex(callset['variants/POS'])
    panel_df = pd.read_csv(panel_file, sep = "\t")

    def test_catch_pop_superpop_both_specified(self):
        with self.assertRaises(ValueError):
            locate_panel_individuals(self.callset, self.panel_file, pop = "GBR", super_pop = "EUR")


    def test_catch_pop_not_exist(self):
        with self.assertRaises(ValueError):
            locate_panel_individuals(self.callset, self.panel_file, pop = "ABE")
            

    def test_subset_AFR(self):
        loc_samples = locate_panel_individuals(self.callset, self.panel_file, super_pop = "AFR")
        self.assertTrue(len(loc_samples) == sum(self.panel_df.super_pop == "AFR"))
        for ind in self.callset["samples"][loc_samples]:
            self.assertTrue(self.panel_df[self.panel_df["sample"] == ind]["super_pop"].tolist()[0] == "AFR")

    def test_all_panel(self):
        loc_samples = locate_panel_individuals(self.callset, self.panel_file)
        self.assertTrue(len(loc_samples) == 2504)

    def test_catch_no_snp_in_range(self):
        with self.assertRaises(KeyError):
            locate_genotype_region(self.pos_array, 0, 1)

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


if __name__ == '__main__':
    unittest.main()

            