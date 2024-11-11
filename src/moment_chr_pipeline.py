import parse_vcf
import numpy as np
import allel
import moments.LD


def initialize_list(n, dtype) -> np.array:
    arr = np.empty(n)
    arr.fill( np.nan)
    return arr


def record_moments_LD(zarr_path, window_list, threshold = 0, panel_file = None, super_pop = None):
    n_windows = len(window_list)
    assert n_windows >= 1 
    D2_list = initialize_list(n_windows, float)
    Dz_list = initialize_list(n_windows, float)
    D_list = initialize_list(n_windows, float)
    pi2_list = initialize_list(n_windows, float)

    D2_filter_list = initialize_list(n_windows, float)
    Dz_filter_list = initialize_list(n_windows, float)
    D_filter_list = initialize_list(n_windows, float)
    pi2_filter_list = initialize_list(n_windows, float)

    pair_count_list = initialize_list(n_windows, int)
    pair_count_filter_list = initialize_list(n_windows, int)

    for i in range(n_windows):
        pos_start = window_list[i][0]
        pos_end = window_list[i][1]
        genotype_012, pos_array = parse_vcf.get_genotype012(zarr_path, pos_start = pos_start, pos_end = pos_end, panel_file = panel_file, super_pop = super_pop)
        if (len(pos_array) == 0):
            continue
        D2_pw, Dz_pw, pi2_pw, D_pw = moments.LD.Parsing.compute_pairwise_stats(genotype_012, genotypes = True)
        for arr, stats in zip([D2_list, Dz_list, D_list, pi2_list], [D2_pw, Dz_pw, D_pw, pi2_pw]):
            arr[i] = np.sum(stats)
        pair_count_list[i] = len(D2_pw)
        D2_pw, Dz_pw, pi2_pw, D_pw = moments.LD.Parsing.compute_pairwise_stats(genotype_012, genotypes = True, pos_array = pos_array, distance_constrained = threshold)
        for arr, stats in zip([D2_filter_list, Dz_filter_list, D_filter_list, pi2_filter_list], [D2_pw, Dz_pw, D_pw, pi2_pw]):
            arr[i] = np.sum(stats)
        pair_count_filter_list[i] = len(D2_pw)
    
    LD_unfiltered_dict = {"D2_pw": D2_list, "Dz_pw": Dz_list, "D_pw": D_list, 
    "pi2_pw": pi2_list, "pair_count": pair_count_list, "windows": window_list}
    LD_filtered_dict = {"D2_pw_subset": D2_filter_list, "Dz_pw_subset": Dz_filter_list, 
    "D_pw_subset": D_filter_list, "pi2_pw_subset": pi2_filter_list, "pair_count": pair_count_filter_list, "windows": window_list}

    return LD_unfiltered_dict, LD_filtered_dict