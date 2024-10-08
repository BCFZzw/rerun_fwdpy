import unittest
import sys
sys.path.insert(1, '/home/alouette/projects/ctb-sgravel/alouette/Dz_sweep_final/filter_snp')
from cluster_snp_filter import *
import numpy as np

class Test_cluster_snp_filter(unittest.TestCase):

    def test_sort_array(self):
        arr1 = [1, 2, 3, 4, 5]
        self.assertTrue((sort_array(arr1) == arr1).all())
        arr2 = [2, 3, 1, 5, 4]
        self.assertTrue((sort_array(arr2) == arr1).all())
        arr3 = [1, 1, 2, 2, 3]
        self.assertTrue((sort_array(arr3) == arr3).all())
        arr4 = [2, 1, 3, 2, 1]
        self.assertTrue((sort_array(arr4) == arr3).all())
        arr3 = [[1, 2], 3, 4, 5]
        ### catch assertion errors

    def test_index_in_pairwise_calculation(self):
        arr1 = [0, 1, 2, 3, 4]
        n1 = len(arr1)
        self.assertTrue(index_in_pairwise_calculation(n1, 0, 1) = )


if __name__ == '__main__':
    unittest.main()