import unittest
import sys
sys.path.insert(1, '/home/alouette/projects/ctb-sgravel/alouette/Dz_sweep_final/simulation_review/src')
import windowing
import numpy as np 
import pandas as pd

class Test_windowing(unittest.TestCase):

    ### @TODO: add a test regarding rounding of position close to the end


    rate_map_chr1 = "/home/alouette/projects/ctb-sgravel/data/genetic_maps/HapMapII_GRCh38/genetic_map_Hg38_chr1.txt"
    chromosome = "chr1"
    prev_window_path = "/home/alouette/projects/ctb-sgravel/alouette/fwdpy_data/1000GfixedBin_raw/all_chr_ALLALL_LD.txt"
    window_df = pd.read_csv(prev_window_path, sep = "\t")
    window_df = window_df[window_df.chr == chromosome].sort_values(by = ["chr", "start"]).reset_index(drop = False)


    def test_0_windows_or_step_too_large(self):
        pos_start = 55550
        pos_end = 68965
        with self.assertRaises(ValueError):
            test = windowing.window_by_recombination(self.rate_map_chr1, pos_start = pos_start, pos_end = pos_end, rec_step = 0.8)
        
    def test_catch_reverse_range(self):
        pos_start = 55550
        pos_end = 68965
        with self.assertRaises(ValueError):
            test = windowing.window_by_recombination(self.rate_map_chr1, pos_start = pos_end, pos_end = pos_start)

    def test_catch_end_range_large(self):
        pos_end = 1e10
        with self.assertRaises(ValueError):
            test = windowing.window_by_recombination(self.rate_map_chr1, pos_end = pos_end)
        
    def test_only_1_window(self):
        pos_start = 55550
        pos_end = 68965
        ### range 0.04
        test = windowing.window_by_recombination(self.rate_map_chr1, pos_start = pos_start, pos_end = pos_end)
        self.assertTrue(np.shape(test)[0] == 1)

        
    def test_region_window(self):
        pos_start = 68965
        pos_end = 308604
        ### the last window will not be included 
        test = windowing.window_by_recombination(self.rate_map_chr1, pos_start = None, pos_end = None)



if __name__ == '__main__':
    unittest.main()

            