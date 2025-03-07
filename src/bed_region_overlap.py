import pandas as pd
import numpy as np
import pybedtools

#### In script please specify bedtools location
### On compute canada this is 
### pybedtools.helpers.set_bedtools_path(path='/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/bedtools/2.31.0/bin/')


def read_bedfile(path: str) -> pd.DataFrame :
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
    intersect = bed1.intersect(bed2, **kwargs)
    intersect_df = intersect.to_dataframe()
    return intersect_df

def intersect_windows(LD_bed, bed) -> pd.DataFrame:
    intersect_df = intersect_2_bed(LD_bed, bed, wo = True)
    


def instersect_Nea_tracks(df, bedfile):
    LD_window = convert_dataframe_to_bed(df)
    intersect_Nea = intersect_2_bed(LD_window, bedfile, wo = True).to_dataframe()
    ### sum all overlapped bases
    intersect_Nea.columns = ["chr", "window_start", "window_end", "bed_chr", "bed_start", "bed_end", "freq", "overlap"]
    intersect_window = intersect_Nea.groupby(["chr", "window_start", "window_end"]).agg({ "overlap": "sum", "freq": "max"}).reset_index()
    intersect_window["window_size"] = intersect_window.window_end - intersect_window.window_start
    intersect_window["overlap_ratio"] = intersect_window["overlap"]/intersect_window.window_size
    return intersect_window

#def window_overlap_corelation()