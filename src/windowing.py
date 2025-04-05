import pandas as pd 
import msprime as ms
import numpy as np


def msprime_read_HapMap(path: str, pos_col = 1, rate_col = 2) -> ms.RateMap:
    """
    Require msprime >= 1.0.0
    Read the recombination map using msprime functionality.
    Can uses the functionality: get_cumulative_mass(), find_index()
    Column 1-2, positions; left: inclusive, right: exclusive
    The HapMap format input must be a file with a header and at least 2 columns, delimited by white space.
    The 2 columns must be Position(bp), and Rate(cM/Mb). Default for them to take the second and third columns. 
    The return will be a msprime RateMap object, with position in bp, and rate in M/bp.
    """
    ms_RateMap = ms.RateMap.read_hapmap(path, position_col= pos_col , rate_col = rate_col)
    return ms_RateMap


def get_bp_from_cum_cM(cM_list: list, ms_RateMap: ms.RateMap) -> np.ndarray:
    """
    Following msprime RateMap get_cumulative_mass functionality for the bp to cumulative cM convertion.
    """
    M_array = np.array(cM_list)/100 ### convert to Morgan
    if np.any(M_array < 0) or np.any(M_array > ms_RateMap.total_mass):
        raise ValueError(f"Cannot have cumulative cM < 0 or > {ms_RateMap.total_mass * 100}")
    ### np.interp(query, given_x, given_f(x))
    pos_array = np.interp(M_array, ms_RateMap.get_cumulative_mass(ms_RateMap.right), ms_RateMap.right)
    pos_array = pos_array.astype(int)
    return pos_array


def get_cum_cM_from_bp(pos_list: list, ms_RateMap: ms.RateMap) -> np.ndarray:
    """
    Using msprime RateMap functionality
    """
    return ms_RateMap.get_cumulative_mass(pos_list)*100




def window_by_reombination(snp_array: list, pos_array: list, rec_map: ms.RateMap, rec_start: float, rec_end: float, rec_step = 0.04) -> pd.DataFrame:
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
    

def window_to_bed(savePath: str, window_df: pd.DataFrame, chromosome = "."):
    """
    Save the windows to a bed file with header and 4 columns: chr, window_start, window_end, cum_rec_rate_at_start(cM)
    """
    df = window_df.copy()
    df.insert(loc=0, column='chr', value= chromosome)
    df.to_csv(savePath, sep = "\t", index = False, header = ["chr", "window_start", "window_end", "cum_rec_rate_at_start(cM)"])