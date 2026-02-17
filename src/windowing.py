import pandas as pd 
import msprime as ms
import numpy as np


def msprime_read_Plink(path: str) -> ms.RateMap:
    """
    Require msprime >= 1.0.0
    Read Plink format, which is a cumulative sum of recombiantion rates.
    """
    df = pd.read_csv(path, names = ["chr", "dot", "cum_mass", "pos"], sep = " ")
    pos = df.pos
    rate = df.cum_mass.diff() / df.pos.diff() / 100 
    ### Adding 0 for pos and 0 for rate as padding directly read RateMap
    pos = [0] + pos.to_list()
    rate = [0] + rate[1:].to_list()
    ms_RateMap = ms.RateMap(position= pos , rate = rate)
    return ms_RateMap



def msprime_read_HapMap(path: str) -> ms.RateMap:
    """
    Require msprime >= 1.0.0
    Read the recombination map using msprime functionality.
    Can uses the functionality: get_cumulative_mass(), find_index()
    Column 1-2, positions; left: inclusive, right: exclusive
    The HapMap format input must be a file with a header and at least 2 columns, delimited by white space.
    The 2 columns must be Position(bp), and Rate(cM/Mb). Default for them to take the second and third columns. 
    The return will be a msprime RateMap object, with position in bp, and rate in M/bp.
    """
    ms_RateMap = ms.RateMap.read_hapmap(path, position_col= 1 , rate_col = 2)
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


def window_by_recombination(rec_map_path: str, rate_map_type: str, rec_step = 0.04, pos_start = None, pos_end = None) -> pd.DataFrame:
    """
    Return a dataframe of non-overlapping windows of fixed recombination rate calculated from the provided rate map. Units are in centi-Morgan.
    The user can specify a position start and end to calculate the windows.
    The first window will be equal to the rate map start or position start, whichever is larger.
    The last window will be last full-length window of rec_step size, equal or smaller than rate map end or position end.
    """
    if rate_map_type == "HapMap":
        ms_RateMap = msprime_read_HapMap(rec_map_path)
    elif rate_map_type == "Plink":
        ms_RateMap = msprime_read_Plink(rec_map_path)
    else:
        raise ValueError("The type of recombination map is not 'HapMap' or 'Plink'.")
    if (pos_start is None):
        pos_start = 0
    if (pos_end is None):
        pos_end = ms_RateMap.sequence_length
    if pos_end < pos_start:
        raise ValueError("Input positions are invalid.")
    if pos_end > ms_RateMap.sequence_length:
        raise ValueError("Position beyond rate map end.")
    ### in RateMap.slice(), right side is exclusive, plus 1 to include
    ### except when it is already the entire length
    trim_RateMap = ms_RateMap.slice(pos_start, min(pos_end + 1, ms_RateMap.sequence_length), trim = False)
    ### generate fixed rec rate windows
    cM_list = np.arange(0, trim_RateMap.total_mass * 100, rec_step)
    ### assert if step size is too large
    if len(cM_list) <= 1:
        raise ValueError("No windows will be generated. Step size can be too large")
    ### get bp positions
    bin_list = get_bp_from_cum_cM(cM_list, trim_RateMap)
    ### If only one window is generated (with the initial 0)
    if len(cM_list) == 2:
        window_df = pd.DataFrame({"window_start": bin_list[0], "window_end": bin_list[1], "cum_cM": cM_list[1]}, index = [0])
    else:
        window_df = pd.DataFrame({"window_start": bin_list[:-1], "window_end": bin_list[1:], "cum_cM": cM_list[1:]})
    return window_df
    

def window_to_bed(savePath: str, window_df: pd.DataFrame, chromosome = "."):
    """
    Save the windows to a bed file with header and 4 columns: chr, window_start, window_end, cum_rec_rate_at_start(cM)
    """
    df = window_df.copy()
    df.insert(loc=0, column='chr', value= chromosome)
    df.to_csv(savePath, sep = "\t", index = False, header = ["chr", "window_start", "window_end", "cum_rec_rate_at_start(cM)"], float_format='%.2f')