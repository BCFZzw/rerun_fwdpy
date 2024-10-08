import itertools
import numpy as np

def sort_array(array_1d: list) -> np.array:
    array_np = np.array(array_1d)
    assert array_np.ndim == 1
    return np.sort(array_np)

def index_in_pairwise_calculation(n, index_a, index_b):
    if index_a > index_b:
        index_a, index_b = index_b, index_a
    assert index_b <= n
    assert index_a < index_b
    assert index_a >= 0
    # assume index starts from 0
    # iterate the number of elements before a, each - 1
    index_pairwise = 0
    for i in range(0, index_a):
        index_pairwise += n - i - 1

    index_pairwise += index_b - index_a - 1
    return index_pairwise

def pairwise_array(array_1d_sorted: np.array) -> np.ndarray:
    assert array_1d_sorted.ndim == 1
    ### check if strictly increasing and sorted
    assert sum(np.diff(array_1d_sorted) > 0)
    pairwise_tuple_list = list(itertools.combinations(array_1d_sorted, 2))
    pairwise_array_list = np.array([list(pairwise_tuple) for pairwise_tuple in pairwise_tuple_list])
    num_element = len(array_1d_sorted)
    assert np.shape(pairwise_array_list) == (num_element*(num_element -1)/2, 2)
    return pairwise_array_list


def boolean_pairwise_filter_list(pairwise_2darray: np.ndarray, threshold: int) -> np.array:
    assert pairwise_2darray.ndim == 2
    pairwise_diff_2darray = np.diff(pairwise_2darray, axis = 1)
    pairwise_diff_array = np.array([array[0] for array in pairwise_diff_2darray])
    pariwise_boolean_array = pairwise_diff_array >= threshold
    assert pariwise_boolean_array.ndim == 1
    return pariwise_boolean_array


def filter_pairwise_array(LD_array: np.array, boolean_array: np.array) -> np.array:
    assert LD_array.ndim == 1
    assert boolean_array.ndim == 1
    assert len(LD_array) == len(boolean_array)
    return LD_array[boolean_array]