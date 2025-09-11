import pandas as pd
import numpy as np
import pybedtools

#### In script please specify bedtools location
### On compute canada this is 
### pybedtools.helpers.set_bedtools_path(path='/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/bedtools/2.31.0/bin/')


def read_bedfile(path: str) -> pd.DataFrame :
    ### Requires bed file with no header
    bed = pybedtools.BedTool(path)
    bed_df = bed.to_dataframe()
    return bed_df 

def convert_dataframe_to_bed(df: pd.DataFrame):
    ### currenlty taking 3-4 columns data
    bed = pybedtools.BedTool.from_dataframe(df)
    return bed

def intersect_2_bed(bed1, bed2, **kwargs) -> pd.DataFrame:
    """
    General function to scall pybedtools intersect 
    **kwargs to account all non specific parameters for bedtools
    Require pybedtools.BedTool or str format for bedfile inputs.
    """
    bed1 = pybedtools.BedTool(bed1)
    bed2 = pybedtools.BedTool(bed2)
    intersect_bed = bed1.intersect(bed2, **kwargs)
    return intersect_bed

def bed_to_dataframe(bed, **kwargs):
    return bed.to_dataframe(**kwargs)

def replace_na(df, col, col_type):
    """
    Sometimes bedtools intersect will return "." to indicate no data.
    When converting to the dataframe, the "." is seen as a string, and often makes the column object type.
    This function converts "." as None in the column and re-convert the type for proper analysis
    """
    df_copy = df.copy()
    df_copy[col] = df_copy[col].replace(".", None).astype(col_type)
    return df_copy


def intersect_windows(df, bed, rename_cols = None, **kwargs) -> pd.DataFrame:
    """
    Assume non-overlapped features in the bed file.
    Using kwargs to specify which bed operations to use.
    """
    df_bed = convert_dataframe_to_bed(df)
    if check_overlapping_features(bed):
        raise ValueError("Bed file contains overlapping features.")
    intersect_df = bed_to_dataframe(intersect_2_bed(df_bed, bed, **kwargs), names = rename_cols)
    return intersect_df

def group_intersect_results(intersect_df, group_ops: dict):
    """
    Using group_ops to specify the dict in pd.DataFrame.agg() functions. Grouping will be done using the first three columns, which is the ["chrom", "start", "end"] in the default pybedtools intersect output . Grouping Dict should be of form {column_index : function operation}
    """
    cols = intersect_df.columns
    assert len(cols) > 3
    position_cols = list(cols[:3])
    sorted_keys = list(group_ops.keys())
    ### sort the column index in order to make sure output is in order
    agg_dict = {}
    for key in sorted_keys:
        if not (key in cols[3:]):
            raise ValueError("Aggregate keys not in ungroupped columns")
        agg_dict[key] = group_ops[key]
    ### the first 3 columns will always be this
    group_df = intersect_df.groupby(position_cols).agg(agg_dict).reset_index()
    return group_df


def rename_df(df, rename_cols: list):
    """
    bedtools by default names the merged file with bed columns 1-12 (what is after 1-12?)
    """
    assert np.shape(df.columns) == np.shape(rename_cols)
    df.columns = rename_cols

    

def intersect_windows_get_overlap(df, bed, group_ops: dict) -> pd.DataFrame:
    """
    Deprecated
    -wo is used in bedtools intersect to get the ratio of overlap.
    Assumes non-overlapped bed track features.
    """
    intersect_df = intersect_windows(df, bed, wo = True)
    ### column 6 is the number of bp overlapped
    ratio_df = group_intersect_results(intersect_df, group_ops = group_ops)
    ratio_df.columns = ["chr", "window_start", "window_end", "overlap"]
    ratio_df["window_size"] = ratio_df["window_end"] - ratio_df["window_start"]
    ratio_df["overlap_ratio"] = ratio_df["overlap"]/ratio_df["window_size"]
    return ratio_df

def intersect_windows_get_overlap_freq(df, bedfile, group_ops:dict) -> pd.DataFrame:
    """
    Deprecated
    -wo is used in bedtools intersect to get the ratio of overlap and maximum frequency of the overalpped tracks.
    Assumes non-overlapped bed track features.
    """
    intersect_df = intersect_windows(df, bed, wo = True)
    ### column 6 is frequency
    ### column 7 is the number of bp overlapped
    freq_df = group_intersect_results(intersect_df, group_ops = group_ops)
    freq_df.columns = ["chr", "window_start", "window_end", "freq", "overlap"]
    freq_df["window_size"] = freq_df["window_end"] - freq_df["window_start"]
    freq_df["overlap_ratio"] = freq_df["overlap"]/freq_df["window_size"]
    return freq_df


def check_overlapping_features(bed) -> bool:
    """
    Check if all bed file features are non-overlapping.
    Bed regions are left-open and right-closed.
    Features are checked per chromosomes.
    Taking the first 3 columns as chrom, start, end for checking overlaps.
    """
    df = read_bedfile(bed)
    cols = df.columns
    col1 = cols[0]
    col2 = cols[1]
    col3 = cols[2]
    for chrom in df[col1].unique():
        chr_region_start = df[df[col1] == chrom][col2].to_list()
        chr_region_end = df[df[col1] == chrom][col3].to_list()
        intervals = [pd.Interval(*regions, closed = "right") for regions in zip(chr_region_start, chr_region_end)]
        interval_arr =  pd.arrays.IntervalArray(intervals)
        ### sort index
        interval_index = pd.IntervalIndex(interval_arr)
        sorted_index = interval_index.sort_values()
        sorted_interval_arr = pd.arrays.IntervalArray(sorted_index)
        if not sorted_interval_arr.is_non_overlapping_monotonic:
            return True
    return False


def instersect_Nea_tracks(df, bedfile):
    LD_window = convert_dataframe_to_bed(df)
    intersect_Nea = intersect_2_bed(LD_window, bedfile, wo = True)
    ### sum all overlapped bases
    intersect_Nea.columns = ["chr", "window_start", "window_end", "bed_chr", "bed_start", "bed_end", "freq", "overlap"]
    intersect_window = intersect_Nea.groupby(["chr", "window_start", "window_end"]).agg({ "overlap": "sum", "freq": "max"}).reset_index()
    intersect_window["window_size"] = intersect_window.window_end - intersect_window.window_start
    intersect_window["overlap_ratio"] = intersect_window["overlap"]/intersect_window.window_size
    return intersect_window




#def window_overlap_correlation()


# def differential_plot