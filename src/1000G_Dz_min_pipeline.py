from parse_vcf import *
import moments.LD
import pandas as pd 
import numpy as np
import os


zarr_path = "/home/alouette/projects/ctb-sgravel/data/30x1000G_biallelic_strict_masked/zarrFormat/chr22"
panel_file = "/home/alouette/projects/ctb-sgravel/alouette/1000Genome/population_panel/integrated_call_samples_v3.20130502.2504.ALL.panel"
window_path = "/home/alouette/projects/ctb-sgravel/alouette/fwdpy_data/1000GfixedBin_raw/all_chr_ALLALL_LD.txt"
save_path = "/home/alouette/projects/ctb-sgravel/alouette/Dz_sweep_final/simulation_review/src"
threshold = 1000

D2_list = []
Dz_list = []
D_list = []
pi2_list = []
D2_subset_list = []
Dz_subset_list = []
D_subset_list = []
pi2_subset_list = []

pair_count_list = []
pair_filtered_list = []

window_df = pd.read_csv(window_path, sep = "\t")
window_df = window_df[window_df.chr == "chr22"].sort_values(by = ["chr", "start"]).reset_index(drop = False)

i = 0
for pos_start, pos_end in zip(window_df.start.tolist(), window_df.end.tolist()):
    #if i == 50:
    #    break
    try:
        genotype_012, pos_array = scikit_allele_parse_genotypes(zarr_path, pos_start = pos_start, pos_end = pos_end, panel_file = panel_file)
    except KeyError:
        for arr in [D2_list, Dz_list, D_list, pi2_list, D2_subset_list, Dz_subset_list, D_subset_list, pi2_subset_list]:
            arr.append(np.nan)
        pair_count_list.append(0)
        pair_filtered_list.append(0)
        continue
    D2_pw, Dz_pw, pi2_pw, D_pw = moments.LD.Parsing.compute_pairwise_stats(genotype_012, genotypes = True)
    pair_count_list.append(len(D2_pw))
    for arr, stats in zip([D2_list, Dz_list, D_list, pi2_list], [D2_pw, Dz_pw, D_pw, pi2_pw]):
            arr.append(np.mean(stats))
    D2_pw_subset, Dz_pw_subset, pi2_pw_subset, D_pw_subset = moments.LD.Parsing.compute_pairwise_stats(genotype_012, pos_array, genotypes = True, distance_constrained = threshold)
    for arr, stats in zip([D2_subset_list, Dz_subset_list, D_subset_list, pi2_subset_list], [D2_pw_subset, Dz_pw_subset, D_pw_subset, pi2_pw_subset]):
        ### if all pairwise are filtered 
        if (len(stats) == 0):
            arr.append(np.nan)
        else:
            arr.append(np.mean(stats))
    pair_filtered_list.append(len(D2_pw_subset))
    i = i + 1


LD_unfiltered_dict = {"D2_pw": D2_list, "Dz_pw": Dz_list, 
"D_pw": D_list, "pi2_pw": pi2_list, "pair_count": pair_count_list, "windows": list(zip(window_df.start.tolist(), window_df.end.tolist()))}
LD_filtered_dict = {"D2_pw_subset": D2_subset_list, "Dz_pw_subset": Dz_subset_list, 
"D_pw_subset": D_subset_list, "pi2_pw_subset": pi2_subset_list, "pair_filtered_count": pair_filtered_list, "windows": list(zip(window_df.start.tolist(), window_df.end.tolist()))}

np.save(os.path.join(save_path, "output", "chr22_all_LD_unfiltered.dict.npy"), LD_unfiltered_dict, allow_pickle = True)
np.save(os.path.join(save_path, "output", "chr22_all_LD_filtered.dict.npy"), LD_filtered_dict, allow_pickle = True)