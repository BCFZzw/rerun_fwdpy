import pandas as pd 
import allel
import zarr

def _read_panel(panel_file):
    """
    Read panel list, expect tab-delimited file and the headers "sample", "pop", "super_pop" present.
    """
    panel_df = pd.read_csv(panel_file, sep = "\t")
    assert "sample" in panel_df.columns
    assert "pop" in panel_df.columns
    assert "super_pop" in panel_df.columns
    return panel_df

def locate_sample_individuals(callset, panel_file, pop = None, super_pop = None):
    """
    Locate individuals specified in the panel files in the VCF callset data. Individuals belonging to a specific population or superpopulation can be sampled.
    """
    panel_df = _read_panel(panel_file)
    samples_list = list(callset["samples"][:]) 
    if (pop is not None) and (super_pop is not None):
        raise ValueError("One of population or superpopulation can be specified. Not Both.")
    if (pop is not None):
        panel_df = panel_df[panel_df.pop == pop]
    if (super_pop is not None):
        panel_df = panel_df[panel_df.super_pop == super_pop]
    assert len(panel_df) > 0
    samples_callset_index = np.where(np.in1d(samples_list, panel_df.sample))[0]
    ### no one missing
    assert -1 not in samples_callset_index
    panel_df['callset_index'] = samples_callset_index
    loc_samples = panel_df.callset_index.values
    return loc_samples

