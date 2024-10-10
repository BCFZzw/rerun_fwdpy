import unittest
import sys
import itertools
sys.path.insert(1, '/home/alouette/projects/ctb-sgravel/alouette/Dz_sweep_final/simulation_review/filter_snp')
from cluster_snp_filter import *
import numpy as np
from moments import LD

class Test_cluster_snp_filter(unittest.TestCase):

    def test_sort_already_sorted_array(self):
        arr = [1, 2, 3, 4, 5]
        arr_sorted = arr
        self.assertTrue((sort_array(arr) == arr_sorted).all())

    def test_sort_unsorted_array(self):
        arr = [2, 3, 1, 5, 4]
        arr_sorted = [1, 2, 3, 4, 5]
        self.assertTrue((sort_array(arr) == arr_sorted).all())
    
    def test_sort_multi_element_array(self):
        arr = [1, 1, 2, 2, 3]
        arr_sorted = [1, 1, 2, 2, 3]
        self.assertTrue((sort_array(arr) == arr_sorted).all())

    def test_sort_catch_multi_dim_array(self):
        arr = [[1, 2], [3, 4]]
        ### test example [[1, 2], 3, 4, 5] inhomogeneous array error difficult to catch
        ### skipping for now
        try:
            sort_array(arr)
        except AssertionError:
            self.assertTrue(True)
            return
        self.assertTrue(False)
    
    def test_find_index(self):
        arr =  [2, 3, 1, 5, 4]
        self.assertTrue(find_index(arr, 1) == 2)
    
    def test_find_index_catch_missing(self):
        arr = [0, 1, 2, 3, 4, 5]
        try:
            find_index(arr, 7)
        except AssertionError:
            self.assertTrue(True)
            return
        self.assertTrue(False)
        
    def test_find_index_edge_inhomogenous(self):
        ### will not throw an error
        arr = [0, [1, 1], 2, 3, 4, 5]
        self.assertTrue(find_index(arr, 2) == 2)
    
    def test_find_index_np_array(self):
        arr = np.array([0, 1, 2, 3, 4, 5])
        self.assertTrue(find_index(arr, 3) == 3)

    def test_find_index_multiple_instance(self):
        arr = np.array([0, 1, 1, 3, 4, 5])
        self.assertTrue(find_index(arr, 1) == 1)

    def test_index_in_sorted_list1(self):
        arr = [2, 3, 4, 9, 10]
        pairwise_tuple_list = list(itertools.combinations(arr, 2))
        self.assertTrue(pairwise_tuple_list[index_in_pairwise_list(arr, (2, 3))] == (2, 3))
        
    def test_index_in_sorted_list_bigger_first(self):    
        arr = [2, 3, 4, 9, 10]
        pairwise_tuple_list = list(itertools.combinations(arr, 2))
        self.assertTrue(pairwise_tuple_list[index_in_pairwise_list(arr, (9, 2))] == (2, 9))

    def test_list_pairs_too_close_all(self):
        arr = [1, 100, 150, 180, 250, 400, 600]
        list_pairs = list_pairs_within_threshold(arr, 1000)
        ### should be the same as itertools, order is different
        check_pairs_tuple = list(itertools.combinations(arr, 2))
        check_pairs_list = [[i[0], i[1]] for i in check_pairs_tuple]
        for test in list_pairs:
            self.assertTrue(np.any(check_pairs_list == test))
        self.assertTrue(len(list_pairs) == len(check_pairs_list))

    def test_list_pairs_too_close_exclusive(self):
        arr = [1, 100, 150, 399, 400 , 401, 1400]
        list_pairs = list_pairs_within_threshold(arr, 1000)
        ### should be the same as itertools, order is different
        check_pairs_tuple = list(itertools.combinations(arr, 2))
        check_pairs_list = [[i[0], i[1]] for i in check_pairs_tuple]
        check_pairs_list.remove([1, 1400])
        check_pairs_list.remove([100, 1400])
        check_pairs_list.remove([150, 1400])
        check_pairs_list.remove([399, 1400])
        check_pairs_list.remove([400, 1400])
        for test in list_pairs:
            self.assertTrue(np.any(check_pairs_list == test))
        self.assertTrue(len(list_pairs) == len(check_pairs_list))

    def test_list_pairs_too_close_none(self):
        arr = [1, 1001, 2001, 3001, 4001 , 5001, 6001]
        list_pairs = list_pairs_within_threshold(arr, 1000)
        self.assertTrue(list_pairs == [])

if __name__ == '__main__':
    unittest.main()