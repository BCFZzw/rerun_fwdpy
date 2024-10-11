import numpy as np 
cimport numpy as np

cpdef tally_sparse(dict G1, dict G2, int n, missing=False):
    """
    G1 and G2 are dictionaries with sample indices of genotypes 1 and 2
    and -1 if missing is True
    n is the diploid sample size
    """
    cdef int n22, n21, n20, n2m, n12, n11, n10, n1m, n02, n01, n00, nm1, nm2, nm
    
    if missing == True:
        # account for missing genotypes
        n22 = (G1[2] & G2[2]).__len__()
        n21 = (G1[2] & G2[1]).__len__()
        n2m = (G1[2] & G2[-1]).__len__()
        n20 = (G1[2]).__len__()-n22-n21-n2m
        n12 = (G1[1] & G2[2]).__len__()
        n11 = (G1[1] & G2[1]).__len__()
        n1m = (G1[1] & G2[-1]).__len__()
        n10 = (G1[1]).__len__()-n12-n11-n1m
        nm2 = (G1[-1] & G2[2]).__len__()
        nm1 = (G1[-1] & G2[1]).__len__()
        n02 = (G2[2]).__len__()-n22-n12-nm2
        n01 = (G2[1]).__len__()-n21-n11-nm1
        # total possible is n-len(set of either missing)
        nm = len(G1[-1].union(G2[-1]))
        n00 = (n-nm)-n22-n21-n20-n12-n11-n10-n02-n01
    else:
        n22 = (G1[2] & G2[2]).__len__()
        n21 = (G1[2] & G2[1]).__len__()
        n20 = (G1[2]).__len__()-n22-n21
        n12 = (G1[1] & G2[2]).__len__()
        n11 = (G1[1] & G2[1]).__len__()
        n10 = (G1[1]).__len__()-n12-n11
        n02 = (G2[2]).__len__()-n22-n12
        n01 = (G2[1]).__len__()-n21-n11
        n00 = n-n22-n21-n20-n12-n11-n10-n02-n01
    return (n22, n21, n20, n12, n11, n10, n02, n01, n00)



cpdef count_genotypes_sparse_distance_constrained(dict G_dict,  int n, np.ndarray[np.int32_t, ndim=1] pos_array, int threshold, missing=False):
    """
    
    """
    cdef int L = len(G_dict)
    
    cdef np.ndarray[np.int32_t, ndim=2] Counts = np.empty((L*(L-1)//2, 9), dtype=np.int32)
    cdef int c = 0
    cdef int i,j
    cdef int pair_distance
    
    for i in range(L-1):
        for j in range(i+1,L):
            pair_distance = pos_array[j] - pos_array[i]
            if (pair_distance < threshold):
                continue
            Counts[c] = tally_sparse(G_dict[i], G_dict[j], n, missing=missing)
            c += 1
    return Counts