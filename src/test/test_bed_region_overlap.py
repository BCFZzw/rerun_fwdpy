import unittest
import sys
sys.path.insert(1, '/home/alouette/projects/ctb-sgravel/alouette/Dz_sweep_final/simulation_review/src')
from bed_region_overlap import *
import pybedtools
import numpy as np 
import pandas as pd

class Test_parse_vcf(unittest.TestCase):

    bed1 = "test1.bed"
    bed1_non_overlap = "test1_non_overlap.bed"
    bed2 = "test2.bed"
    bed2_freq = "test2_freq.bed"
    pybedtools.helpers.set_bedtools_path(path='/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/bedtools/2.31.0/bin/')

    def test_intersect_2_bed_u(self):
        df = intersect_2_bed(self.bed1, self.bed2, u = True)
        assert np.all(df == pd.DataFrame({"chrom" : ["chr1"],
            "start" : [10], "end": [100]
        }))

    def test_intersect_2_bed_wo(self):
        df = intersect_2_bed(self.bed1, self.bed2, wo = True)
        assert np.all(df == pd.DataFrame({"chrom" : ["chr1", "chr1"],
            "start" : [10, 10], "end": [100, 100], "name":  ["chr1", "chr1"], "score" : [10, 20], "strand": [20, 30], "thickStart" : [10, 10]
        }))

    def test_intersect_2_bed_v(self):
        bed1 = pybedtools.BedTool(self.bed1)
        bed2 = pybedtools.BedTool(self.bed2)
        df = intersect_2_bed(bed1, bed2, v = True)
        assert np.all(df == pd.DataFrame({"chrom" : ["chr2", "chr3", "chr3"],
            "start" : [10, 10, 10], "end": [200, 20, 100]
        }))

    def test_intersect_windows(self):
        df = read_bedfile(self.bed1)
        intersect_df = intersect_windows(df, self.bed2, v = True)
        assert np.all(intersect_df == pd.DataFrame({"chrom" : ["chr2", "chr3", "chr3"],
            "start" : [10, 10, 10], "end": [200, 20, 100]
        }))

    def test_rename_df(self):
        df = intersect_2_bed(self.bed1, self.bed2, wo = True)
        cols = ["Chr", "Start", "End", "Chr2", "Start2", "End2", "Score"]
        rename_df(df, cols)
        assert np.all(df == pd.DataFrame({"Chr" : ["chr1", "chr1"],
            "Start" : [10, 10], "End": [100, 100], "Chr2":  ["chr1", "chr1"], "Start2" : [10, 20], "End2": [20, 30], "Score" : [10, 10]
        }))

    def test_group_intersect_results(self):
        df = intersect_2_bed(self.bed1, self.bed2, wo = True)
        cols = ["Chr", "Start", "End", "Chr2", "Start2", "End2", "Score"]
        rename_df(df, cols)
        group_df = group_intersect_results(df, {"Score": "first"})
        assert np.all(group_df == pd.DataFrame({"Chr" : ["chr1"],
            "Start" : [10], "End": [100], "Score" : [10]
        }))
        

    def test_check_overlapping_features(self):
        assert check_overlapping_features(self.bed1) == True
        assert check_overlapping_features(self.bed2) == False



if __name__ == '__main__':
    unittest.main()

            