import pandas as pd 
import allel
import zarr
import numpy as np
#import dask

def _read_panel(panel_file):
    """
    Read panel list, expect tab-delimited file and the headers "sample", "pop", "super_pop" present.
    """
    panel_df = pd.read_csv(panel_file, sep = "\t")
    assert "sample" in panel_df.columns
    assert "pop" in panel_df.columns
    assert "super_pop" in panel_df.columns
    return panel_df

def locate_panel_individuals(callset_samples, panel_file, pop = None, super_pop = None):
    """
    Locate individuals specified in the panel files in the VCF callset data. Individuals belonging to a specific population or superpopulation can be sampled.
    """
    panel_df = _read_panel(panel_file)
    samples_list = list(callset_samples[:]) 
    if (pop is not None) and (super_pop is not None):
        raise ValueError("One of population or superpopulation can be specified. Not Both.")
    if (pop is not None):
        panel_df = panel_df[panel_df["pop"] == pop]
    if (super_pop is not None):
        panel_df = panel_df[panel_df["super_pop"] == super_pop]
    if len(panel_df) == 0:
        raise ValueError("Sampled population/superpopulation not found in panel.")
    samples_callset_index = np.where(np.in1d(samples_list, panel_df["sample"]))[0]
    ### no one missing
    if -1 in samples_callset_index:
        raise ValueError("Some panel individuals are missing in the vcf data.")
    panel_df['callset_index'] = samples_callset_index
    loc_samples = panel_df.callset_index.values
    return loc_samples

def locate_genotype_region(pos_array, pos_start: int, pos_end: int):
    """
    Locate the indices of a region in the genotype position array.
    """
    if pos_start is None:
        pos_start = 1
    if pos_end is None:
        pos_end = max(pos_array) + 1
    assert pos_end > pos_start
    try:
        loc_region = pos_array.locate_range(pos_start, pos_end)
    except KeyError:
        raise KeyError("No variants in the given region!")
    return loc_region


def _read_zarr_callset(zarr_path):
    """
    Read zarr callset, expect the genotype to be under "calldata/GT", 
    positions under "variants/POS", and samples under "samples"
    """
    callset = zarr.open_group(zarr_path, mode='r')
    genotype_zarr = callset['calldata/GT']
    callset_pos = callset['variants/POS']
    callset_samples = callset["samples"]
    return callset, genotype_zarr, callset_pos, callset_samples


def scikit_allele_parse_genotypes(zarr_path, pos_start = None, pos_end = None, panel_file = None, pop = None, super_pop = None):
    callset, genotype_zarr, callset_pos, callset_samples = _read_zarr_callset(zarr_path)
    pos_array = allel.SortedIndex(callset_pos)
    genotype_dask = allel.GenotypeDaskArray(genotype_zarr)

    if (pos_start is not None) or (pos_end is not None):
        loc_region = locate_genotype_region(pos_array, pos_start, pos_end)
        pos_array = pos_array[loc_region]
        ### when using indices: take, when using boolean: compress
        genotype_dask = genotype_dask.take(loc_region, axis = 0)

    if panel_file is not None:
        loc_samples = locate_panel_individuals(callset_samples, panel_file, pop, super_pop)
        genotype_dask = genotype_dask.take(loc_samples, axis = 1)
        
    genotype_012_dask = genotype_dask.count_alleles()
    
    return genotype_012_dask, pos_array
