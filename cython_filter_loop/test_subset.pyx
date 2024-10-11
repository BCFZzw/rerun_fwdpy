#define NPY_NO_DEPRECATED_API NPY_1_25_API_VERSION
import numpy as np 
cimport numpy as np

cpdef test_subsetting():
    """
    Similar to count_genotypes, but using the sparse genotype representation instead
    """
    cdef np.ndarray[np.uint8_t, ndim = 1] bool_array
    bool_array = np.array([True, False, True, False, True], dtype = np.bool_)

    print(bool_array)

    cdef np.ndarray[np.int8_t, ndim = 1] filter_array
    filter_array = np.array([1, 2, 3, 4, 5], dtype = np.int8)
    
    return filter_array[bool_array]