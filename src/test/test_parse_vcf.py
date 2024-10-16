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



if __name__ == '__main__':
    unittest.main()

            