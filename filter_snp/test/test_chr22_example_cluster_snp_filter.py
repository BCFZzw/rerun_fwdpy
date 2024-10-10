import unittest
import sys
sys.path.insert(1, '/home/alouette/projects/ctb-sgravel/alouette/Dz_sweep_final/simulation_review/filter_snp')
from cluster_snp_filter import *
import numpy as np
from moments import LD
import zarr 
import allel
import pandas as pd

class Test_chr22_example_cluster_snp_filter(unittest.TestCase):

    zarr_path = "/home/alouette/projects/ctb-sgravel/data/30x1000G_biallelic_strict_masked/zarrFormat/chr22"
    callset = zarr.open_group(zarr_path, mode='r')
    pos_array = allel.SortedIndex(callset['variants/POS'])
    snp_array = callset['calldata/GT']
    pos_first_10 = pos_array[:10]
    ###[10661612, 10666327, 10666344, 10666378, 10669347, 10669364, 10680028, 10680035, 10680045, 10680047]###
    ### Difference:
    ###[ 4715    17    34  2969    17 10664     7    10     2]
    snp_first_10 = snp_array[:10]

    panel_2504 = "/home/alouette/projects/ctb-sgravel/alouette/1000Genome/population_panel/integrated_call_samples_v3.20130502.2504.ALL.panel"
    panel_df = pd.read_csv(panel_2504, sep = "\t", usecols = [0, 1, 2])
    samples = callset['samples'][:]
    samples_list = list(samples)
    samples_callset_index = np.where(np.in1d(samples_list, panel_df['sample']))[0]
    panel_df['callset_index'] = samples_callset_index
    loc_samples_unrelated= panel_df.callset_index.values

    gt_region = allel.GenotypeArray(snp_first_10)
    gt_unrelated = gt_region.take(loc_samples_unrelated, axis=1)
    gt_unrelated_012 = gt_unrelated.to_n_alt(fill=-1)

    D2_pw, _, _, _ = LD.Parsing.compute_pairwise_stats(gt_unrelated_012, 
    genotypes = True)

    def test_filtering_5bp_last_removed(self):
        filtered_pairs = list_pairs_within_threshold(self.pos_first_10, 5)
        filtered_index = bool_list_filtering(self.pos_first_10, filtered_pairs)
        filtered_D2_pw = self.D2_pw[filtered_index]
        self.assertTrue(len(filtered_D2_pw) == (len(self.D2_pw) -1))
        self.assertTrue(np.all(filtered_D2_pw == self.D2_pw[:-1]))

    def test_filtering_2bp_none_removed(self):
        filtered_pairs = list_pairs_within_threshold(self.pos_first_10, 2)
        filtered_index = bool_list_filtering(self.pos_first_10, filtered_pairs)
        filtered_D2_pw = self.D2_pw[filtered_index]
        self.assertTrue(np.all(filtered_D2_pw == self.D2_pw))
    
    def test_filtering_1000000bp_all_removed(self):
        filtered_pairs = list_pairs_within_threshold(self.pos_first_10, 1000000)
        filtered_index = bool_list_filtering(self.pos_first_10, filtered_pairs)
        filtered_D2_pw = self.D2_pw[filtered_index]
        self.assertTrue(np.all(filtered_D2_pw == []))

if __name__ == '__main__':
    unittest.main()

