import pandas as pd 
import allel
import zarr

def _read_panel(panel_file):
    """
    Read panel list, expect space-delimited file and the headers "SampleID", "Population", "Superpopulation" present.
    """
    panel_df = pd.read_csv(panel_file, sep = " ")
    assert ["SampleID", "Population", "Superpopulation"] in panel_df.columns
    return panel_df

def locate_sample_individuals(callset, panel_file, population = None, superpopulation = None):
    """
    Locate individuals specified in the panel files in the VCF callset data. Individuals belonging to a specific population or superpopulation can be sampled.
    """
    panel_df = _read_panel(panel_file)
    samples_list = list(samples) 
    if (population is not None) and (superpopulation is not None):
        raise ValueError("One of population or superpopulation can be specified. Not Both.")
    if (population is not None):
        panel_df = panel_df[panel_df.Population == population]
    if (superpopulation is not None):
        panel_df = panel_df[panel_df.Superpopulation == superpopulation]
    assert len(panel_df) > 0
    samples_callset_index = np.where(np.in1d(samples_list, panel_df.SampleID))[0]
    ### no one missing
    assert -1 not in samples_callset_index
    panel_df['callset_index'] = samples_callset_index
    loc_samples = panel_df.callset_index.values
    return loc_samples

