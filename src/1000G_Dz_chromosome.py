import argparse
import sys 
sys.path.insert(1, '/home/alouette/projects/ctb-sgravel/alouette/Dz_sweep_final/simulation_review/src')
from moment_chr_pipeline import *
import os
import pandas as pd 
import numpy as np

parser = argparse.ArgumentParser(
                    prog='Dz_moment_chromosome_pipeline',
                    description='place_holder'
                   )
parser.add_argument('-c', '--chromosome', dest = "chr", help = "The chromosome to parse LD.", required = True)
parser.add_argument('-s', '--super_pop', dest = "super_pop", help = "Super population to sample.", required = False, default = None)
args = parser.parse_args()

chromosome = args.chr
super_pop = args.super_pop
threshold = 1000
save_path = "/home/alouette/projects/ctb-sgravel/alouette/Dz_sweep_final/simulation_review/src/output"
zarr_master_path = "/home/alouette/projects/ctb-sgravel/data/30x1000G_biallelic_strict_masked/zarrFormat/"
zarr_path = os.path.join(zarr_master_path, chromosome)
panel_file = "/home/alouette/projects/ctb-sgravel/alouette/1000Genome/population_panel/integrated_call_samples_v3.20130502.2504.ALL.panel"
window_path = "/home/alouette/projects/ctb-sgravel/alouette/fwdpy_data/1000GfixedBin_raw/all_chr_ALLALL_LD.txt"
window_df = pd.read_csv(window_path, sep = "\t")
window_df = window_df[window_df.chr == chromosome].sort_values(by = ["chr", "start"]).reset_index(drop = False)
window_list = list(zip(window_df.start.tolist(), window_df.end.tolist()))

super_pop_prefix = super_pop
if super_pop is None:
    super_pop_prefix = "ALL"


LD_unfiltered_dict, LD_filtered_dict = record_moments_LD(zarr_path, window_list, panel_file = panel_file, threshold = threshold, super_pop = super_pop_prefix)

np.save(os.path.join(save_path, ".".join((chromosome, super_pop_prefix, "unfiltered.dict.npy"))), LD_unfiltered_dict, allow_pickle = True)
np.save(os.path.join(save_path, ".".join((chromosome, super_pop_prefix, "filtered.dict.npy"))), LD_filtered_dict, allow_pickle = True)
