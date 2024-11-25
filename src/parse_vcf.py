import pandas as pd 
import allel
import zarr
import numpy as np

def _read_panel(panel_file) -> pd.DataFrame :
    """
    Read panel list, expect tab-delimited file and the headers "sample", "pop", "super_pop" present.
    """
    panel_df = pd.read_csv(panel_file, sep = "\t")
    assert "sample" in panel_df.columns
    assert "pop" in panel_df.columns
    assert "super_pop" in panel_df.columns
    return panel_df

def get_pops_from_superpop(panel_file, super_pop) -> np.array :
    """
    List all populations within a superpopulation.
    """
    panel_df = _read_panel(panel_file)
    assert super_pop in panel_df["super_pop"].unique()
    return panel_df[panel_df.super_pop == super_pop]["pop"].unique()

def locate_jackknife_individuals(callset_samples, panel_file, super_pop, jackknife_pop):
    """
    Locate the indices of leave-one-population-out individuals. 
    Used in jackknife estimation.
    """
    loc_super_pop_samples = locate_panel_individuals(callset_samples, panel_file, super_pop = super_pop)
    loc_pop_samples = locate_panel_individuals(callset_samples, panel_file, pop = jackknife_pop)
    if len(np.intersect1d(loc_super_pop_samples, loc_pop_samples)) == 0:
        raise ValueError("Jackknife population does not exist in the specified superpopulation.")
    if len(np.intersect1d(loc_super_pop_samples, loc_pop_samples)) != len(loc_pop_samples):
        raise ValueError("Jackknife population does not fully overalp the specified superpopulation.")
    loc_jackknife = np.setxor1d(loc_super_pop_samples, loc_pop_samples, assume_unique = False)
    assert len(loc_super_pop_samples) - len(loc_pop_samples) == len(loc_jackknife)
    return loc_jackknife


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
    ### If missing
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
    ### no SNPs in region
    if np.all((pos_array > pos_start) & (pos_array < pos_end) == False):
        return slice(None, -0, None)
    else:
        loc_region = pos_array.locate_range(pos_start, pos_end)
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

def select_variants(genotype_dask, pos_array, pos_start, pos_end):
    """
    Subset genotype based on given positions.
    """
    loc_region = locate_genotype_region(pos_array, pos_start, pos_end)
    pos_array_var = pos_array[loc_region]
    ### when using indices: take, when using boolean: compress
    genotype_var = genotype_dask.take(loc_region, axis = 0)
    return genotype_var, pos_array_var

def select_sample(genotype_dask, callset_samples, panel_file, pop = None, super_pop = None):
    """
    Subset genotype based on a given population or superpopulation.
    """
    loc_samples = locate_panel_individuals(callset_samples, panel_file, pop, super_pop)
    genotype_sample = genotype_dask.take(loc_samples, axis = 1)
    return genotype_sample

def select_sample_jackknife(genotype_dask, callset_samples, panel_file, super_pop, jackknife_pop):
    """
    For Jackknife, subsect genotype based on leave-one-out populations from a superpopulation.
    """
    loc_jackknife = locate_jackknife_individuals(callset_samples, panel_file, super_pop = super_pop, jackknife_pop = jackknife_pop)
    genotype_jackknife = genotype_dask.take(loc_samples, axis = 1)
    return genotype_jackknife


def select_sample_variants(zarr_path, pos_start = None, pos_end = None, panel_file = None, pop = None, super_pop = None, jackknife_pop = None):
    """
    Return the genotype array of selected samples and individuals.
    Either population or superpopulation can be specified.
    """
    callset, genotype_zarr, callset_pos, callset_samples = _read_zarr_callset(zarr_path)
    pos_array = allel.SortedIndex(callset_pos)
    genotype_dask = allel.GenotypeDaskArray(genotype_zarr)

    if (pos_start is not None) or (pos_end is not None):
        genotype_dask, pos_array = select_variants(genotype_dask, pos_array, pos_start, pos_end)

    if panel_file is not None:
        if jackknife_pop is not None:
            assert super_pop is not None
            genotype_dask = select_sample_jackknife(genotype_dask, callset_samples, panel_file, super_pop, jackknife_pop)
        else:
            genotype_dask = select_sample(genotype_dask, callset_samples, panel_file, pop = pop, super_pop = super_pop)
    elif (jackknife_pop is not None) or (super_pop is not None) or (pop is not None):
        raise ValueError("Panel file is not given when selecting individuals.")
    
    return genotype_dask, pos_array



def get_genotype012(zarr_path, pos_start = None, pos_end = None, panel_file = None, pop = None, super_pop = None, jackknife_pop = None):
    genotype_dask, pos_array = select_sample_variants(zarr_path, pos_start = pos_start, pos_end = pos_end, panel_file = panel_file, pop = pop, super_pop = super_pop, jackknife_pop = jackknife_pop)

    ### using .compute() to load from dask
    # Count the number of alternative alleles for biallelic sites
    # homogzygous reference is 0
    genotype_012 = genotype_dask.compute().to_n_alt(fill=-1)
    pos_array_np = np.array(pos_array) # required numpy by constrain function
    return genotype_012, pos_array_np


def get_genotype(zarr_path, pos_start = None, pos_end = None, panel_file = None, pop = None, super_pop = None, jackknife_pop = None):
    genotype_dask, pos_array = select_sample_variants(zarr_path, pos_start = pos_start, pos_end = pos_end, panel_file = panel_file, pop = pop, super_pop = super_pop, jackknife_pop = jackknife_pop)

    genotype = genotype_dask.compute()
    pos_array_np = np.array(pos_array)
    return genotype, pos_array_np