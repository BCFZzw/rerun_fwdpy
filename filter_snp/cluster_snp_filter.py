import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

def sort_array(array_1d: list) -> np.array:
    """
    Sort a 1d array.
    """
    array_np = np.array(array_1d)
    assert array_np.ndim == 1
    return np.sort(array_np)

def find_index(arr: list, a: int) -> int:
    """
    Find the index of element in a array.
    Assertion error if the element is not found.
    """
    ### avoiding np.where because it is not ideal for 1d array
    assert a in arr
    ### using python list built-in function
    arr_list = list(arr)
    ### the first instance is returned, assuming only 1 instance is returned
    idx = arr_list.index(a)
    return idx

def index_in_pairwise_list(arr: list, pair_ab: list) -> int:
    """
    Find the index of a given pair [a, b] in the list of a pairwise comparison 
    of a sorted array. 
    The order of the list is the same as itertools.combinations(, 2).
    Assertion error if the query pair is not found in the given array.
    """
    assert len(pair_ab) == 2
    assert np.all(sort_array(arr) == arr)
    a = pair_ab[0]
    b = pair_ab[1]
    ### Assertion error raised if element not found.
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

def list_pairs_within_threshold(pos_array: list, threshold: int) -> list:
    """
    Filter any sorted position pairs < a distance defined by the user.
    The filtered pairs are listed in a list: [[pos_a, pos_b], [pos_c, pos_d] ...].
    Return an empty list if no pairs are filtered.
    """
    ### using numpy sliding window view to get all the pairs
    assert np.all(sort_array(pos_array) == pos_array)
    assert min(pos_array) > 0
    assert threshold > 0
    filter_pair_list = []
    j = 2
    while j <= len(pos_array):
        ### retain the first and last position in the sliding window
        pair_list_j = sliding_window_view(pos_array, j)[:, [0, -1]]
        pair_list_j_diff = np.diff(pair_list_j)[:, 0]
        filter_pair_j_bool = pair_list_j_diff < threshold
        ### if no more pairs need to be filtered
        if sum(filter_pair_j_bool) == 0:
            break
        filter_pair_list.extend(pair_list_j[filter_pair_j_bool])
        j = j + 1
    return filter_pair_list

def bool_list_filtering(pos_array: list, filtered_pairs: list) -> list:
    """
    Generate a 1d boolean array to easily apply filtering to other pairwise arrays.
    The boolean array has the size of total number of pairs from the position array.
    The order is the same as itertools.combinations(, 2) and moments.LD calculations.
    False at the indices for filtered pairs, True for retained pairs.
    Assertion error if any pairs is not found in the given array.
    """
    assert len(filtered_pairs) > 0
    assert len(pos_array) >= 2
    size_pos_array = len(pos_array)
    total_pairs = int(size_pos_array * (size_pos_array -1)/2)
    index_list = [True]*total_pairs
    for pairs in filtered_pairs:
        index_filtered_pairs = index_in_pairwise_list(pos_array, pairs)
        index_list[index_filtered_pairs] = False
    return index_list
