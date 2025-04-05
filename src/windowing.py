import pandas as pd 
import msprime as ms
import numpy as np
import allel


def _msprime_read_HapMap(rec_map_path: str) -> pd.DataFrame:
    """
    Require msprime >= 1.0.0
    Read the recombination map using msprime functionality.
    Can uses the functionality: get_cumulative_mass(), find_index()
    Column 1-2, positions; left: inclusive, right: exclusive
    """
    ms_RateMap = ms.RateMap.read_hapmap(rec_map_path, position_col=1, rate_col=2)
    return ms_RateMap


def _get_bp_from_cum_cM(cM_list: list, ms_RateMap: ms.RateMap) -> np.ndarray:
    """
    Following msprime RateMap get_cumulative_mass functionality for the bp to cumulative cM convertion.
    """
    M_array = np.array(cM_list)/100 ### convert to Morgan
    if np.any(M_array < 0) or np.any(M_array > ms_RateMap.total_mass):
        raise ValueError(f"Cannot have cumulative cM < 0 or > {ms_RateMap.total_mass * 100}")
    pos_array = np.interp(M_array, np.nancumsum(rec_map_df.mass), ms_RateMap.left)
    pos_array = pos_array.astype(int)
    return pos_array


def _get_cum_cM_from_bp(pos_list: list, ms_RateMap: ms.RateMap) -> np.ndarray:
    """
    Using msprime RateMap functionality
    """
    return ms_RateMap.get_cumulative_mass(pos_list)*100


def _scikit_allel_windowing(snp_array: list, snp_pos: list, window_pos_list:list) -> list:
    """
    Master window function, based on the input bp, give back the subset windows in a list
    This function is the core function and used for other detailed window_by_* functions.
    :param snp_array: The SNP array read by scikit-allel format.
    :param window_index_array: The list of tuple (start, end) index of each windows of the SNP array.
    Return a list of slices
    """
    index_snp_pos = allel.SortedIndex(snp_pos)
    if np.any(window_pos_list) < 0 or np.max(window_pos_list) > np.max(pos_array) or np.min(window_pos_list) > np.min(pos_array):
        raise KeyError(f"Window positions out of bounds")
    if not np.all(np.diff(window_pos_list) > 0):
        raise KeyError(f"Window positions are not in asending order")
    window_list = [index_snp_pos.locate_range(start, end) for start, end in zip(window_pos_list[:-1], window_pos_list[1:]) ]
    ### next step requires loading the data, do it when in need
    return window_list


def window_by_reombination(snp_array: list, pos_array: list, rec_map: ms.RateMap, rec_start: float, rec_end: float, rec_step = 0.04) -> list:
    """
    Require msprime >= 1.0.0 to use RateMap
    Call this function to get the window list separated by rec_distance
    How to handle the last position if rec subsetting is not possible?
    """
    ### end has to be larger than start
    if np.any(rec_start < rec_end):
        raise ValueError(f"Input recombination range is invalid")
    ### if given rec distance is smaller or equal to the the step 
    ### return the original array without subsetting
    if rec_end - rec_start <= rec_step:
        rec_list = [rec_start, rec_end]                    
    else:
        rec_list = np.arange(rec_start, rec_end, rec_step)
        ### check the last position
        ### floating point issue, the last position can be larger than end
        if rec_list[-1] > rec_end:
            rec_list[-1] = rec_end
    ### get bp positions
    bp_array = _get_bp_from_cum_cM(rec_list)
    ### call _windowing function
    window_list = _scikit_allel_windowing(snp_array, pos_array, bp_array)
    return window_list
    
    