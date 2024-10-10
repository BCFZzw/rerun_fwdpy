import itertools
import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
#from scipy.sparse import csr_array

def sort_array(array_1d: list) -> np.array:
    array_np = np.array(array_1d)
    assert array_np.ndim == 1
    return np.sort(array_np)

def find_index(arr: list, a: int) -> int:
    ### avoiding np.where because it is not ideal for 1d array
    assert a in arr
    ### using python list built-in function
    arr_list = list(arr)
    ### the first instance is returned, assuming only 1 instance is returned
    idx = arr_list.index(a)
    return idx

def index_in_pairwise_list(arr: list, tuple_pair: list):
    ###@TODO array needs to be non repeating and ascending
    #assert type(tuple_pair) == tuple
    assert len(tuple_pair) == 2
    assert np.array(arr).ndim == 1
    a = tuple_pair[0]
    b = tuple_pair[1]
    index_a = find_index(arr, a)
    index_b = find_index(arr, b)
    n = len(arr)
    assert n >= 2
    if index_a > index_b:
        index_a, index_b = index_b, index_a
    assert index_b <= n
    assert index_a < index_b
    assert index_a >= 0

    # index starts from 0
    ### simplify mathmatically, cast to int
    index_pairwise = int(((n - 1) + (n - index_a))*index_a/2)

    index_pairwise += index_b - index_a - 1
    return index_pairwise

#def bool_index_list(pairs: int) -> list:
#    assert n > 0
#    index_list = [True]*pairs
#    return index_list

#def filter_clustered_SNPs

def list_pairs_within_threshold(pos_array: list, threshold: int) -> list:
    ### numpy sliding window view to generate the difference faster than for loop
    ### stop when all values larger than threshold
    assert np.array(pos_array).ndim == 1 
    filter_pair_list = []
    j = 2
    while j <= len(pos_array):
        ### retain the first and last position in the sliding window
        pair_list_j = sliding_window_view(pos_array, j)[:, [0, -1]]
        pair_list_j_diff = np.diff(pair_list_j)[:, 0]
        filter_pair_j_bool = pair_list_j_diff <= threshold
        ### if no more distance withint threshold
        if sum(filter_pair_j_bool) == 0:
            break
        filter_pair_list.extend(pair_list_j[filter_pair_j_bool])
        j = j + 1
    return filter_pair_list
    