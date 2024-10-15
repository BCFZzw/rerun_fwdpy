import unittest
import zarr
import sys
sys.path.insert(1, '/home/alouette/projects/ctb-sgravel/alouette/Dz_sweep_final/simulation_review/src')
from parse_vcf import *


class Test_parse_vcf(unittest.TestCase):

    panel_file = "/home/alouette/projects/ctb-sgravel/alouette/1000Genome/population_panel/integrated_call_samples_v3.20130502.2504.ALL.panel"
    zarr_path = "/home/alouette/projects/ctb-sgravel/data/30x1000G_biallelic_strict_masked/zarrFormat/chr22"
    callset = zarr.open_group(zarr_path, mode='r')

    def test_catch_pop_superpop_both_specified(self):
        with self.assertRaises(ValueError):
            locate_sample_individuals(self.callset, self.panel_file, pop = "GBR", super_pop = "EUR")



if __name__ == '__main__':
    unittest.main()

            