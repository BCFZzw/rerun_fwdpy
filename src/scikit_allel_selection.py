import allel
import os
import numpy as np
import pandas as pd
import parse_vcf
np.seterr(divide='ignore', invalid='ignore') ### ignore allel fst.py numpy undefined division 

### Probably double-check another paper with how scanning and windows are handled

def genotype_dask_to_allele_count(genotype_dask, pos_array, pos_start, pos_end):
    genotype_var_dask, pos_var_array = parse_vcf.select_variants(genotype_dask, pos_array, pos_start, pos_end)
    if len(pos_var_array) == 0:
        return np.array([]), np.array([])
    genotype_var = genotype_var_dask.compute()
    pos_array_np = np.array(pos_var_array)
    genotype_ac = genotype_var.count_alleles()
    return genotype_ac, pos_array_np

def allel_PBS(zarr_path, panel_file, pop1_ID, pop2_ID, pop3_ID, step = 100000):
    ### Typically done on the single allele level, and normalized with allele windowing methods 
    ### Iteratively go through the allele lists in steps, avoid reading all at once
    ### Also record allele frequency metrics
    genotype_pop1, pos_array = parse_vcf.select_sample_variants(zarr_path, panel_file = panel_file, super_pop = pop1_ID)
    genotype_pop2, _ = parse_vcf.select_sample_variants(zarr_path, panel_file = panel_file, super_pop = pop2_ID)
    genotype_pop3, _ = parse_vcf.select_sample_variants(zarr_path, panel_file = panel_file, super_pop = pop3_ID)

    n_step = (pos_array[-1] - pos_array[0])//step + 1

    PBS_list = []
    pos_list = []

    ### Looping via all windows
    for i in range(n_step):
        pos_start = pos_array[0] + i*step
        if i == n_step :
            pos_end = None
        else:
            pos_end = pos_array[0] + (i+1)*step - 1

        genotype_ac_pop1, pos_ac = genotype_dask_to_allele_count(genotype_pop1, pos_array = pos_array, pos_start = pos_start, pos_end = pos_end)
        if len(pos_ac) == 0 :
            continue
        genotype_ac_pop2, _ = genotype_dask_to_allele_count(genotype_pop2, pos_array = pos_array, pos_start = pos_start, pos_end = pos_end)
        genotype_ac_pop3, _ = genotype_dask_to_allele_count(genotype_pop3, pos_array = pos_array, pos_start = pos_start, pos_end = pos_end)

        PBS_array = allel.pbs(genotype_ac_pop1, genotype_ac_pop2, genotype_ac_pop3, window_size = 1, window_step = 1)
        assert len(PBS_array) == len(pos_ac)
        PBS_list.append(PBS_array) 
        pos_list.append(pos_ac)

    PBS_dict = {"PBS_list": PBS_list, "position": pos_list, "populations": [pop1_ID, pop2_ID, pop3_ID]}

    return PBS_dict