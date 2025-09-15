import numpy as np
import os
import glob
import pandas as pd 

def read_filtered(np_save_path, chromosome):
    LD_dict = np.load(np_save_path, allow_pickle = True)
    Dz_pw = LD_dict.item().get("Dz_pw_subset")
    D2_pw = LD_dict.item().get("D2_pw_subset")
    pi2_pw = LD_dict.item().get("pi2_pw_subset")
    D_pw = LD_dict.item().get("D_pw_subset")
    pair_count = LD_dict.item().get("pair_count")
    windows = LD_dict.item().get("windows")
    df = pd.DataFrame({"window_start": [window[0] for window in windows],
                      "window_end": [window[1] for window in windows],
                       "D2": D2_pw,
                       "Dz": Dz_pw,
                      "D": D_pw,
                      "pi2": pi2_pw,
                      "snp_pairs": pair_count,
                      })
    df["sigD2"] = df.D2/df.pi2
    df["sigDz"] = df.Dz/df.pi2
    df["window_size"] = df.window_end - df.window_start
    df["chr"] = chromosome
    return df

def read_jackknife(np_save_path, chromosome):
    LD_dict = np.load(np_save_path, allow_pickle = True)
    Dz_pw = LD_dict.item().get("Dz_pw")
    D2_pw = LD_dict.item().get("D2_pw")
    pi2_pw = LD_dict.item().get("pi2_pw")
    D_pw = LD_dict.item().get("D_pw")
    pair_count = LD_dict.item().get("pair_count")
    jackknife_pop = LD_dict.item().get("params")["jackknife_pop"]
    windows = LD_dict.item().get("windows")
    df = pd.DataFrame({"window_start": [window[0] for window in windows],
                      "window_end": [window[1] for window in windows],
                       "D2": D2_pw,
                       "Dz": Dz_pw,
                      "D": D_pw,
                      "pi2": pi2_pw,
                      })
    
    df["_".join(["sigD2", "jk", jackknife_pop])] = df.D2/df.pi2
    df["_".join(["sigDz", "jk", jackknife_pop])] = df.Dz/df.pi2
    df.drop(["D2", "Dz", "D", "pi2"], inplace = True, axis = 1)
    df["chr"] = chromosome
    return df

def concatenate_filtered(df_list):
    df_combined = pd.concat(df_list, axis = 0)
    df_combined = df_combined.sort_values(by = ["chr", "window_start"]).reset_index(drop = True)
    return df_combined

def combine_jackknife(jackkife_df_list):
    merge_df = jackkife_df_list[0].copy()
    for df_jackknife in jackkife_df_list[1:]:
        merge_df = merge_df.merge(df_jackknife, on = ["chr", "window_start", "window_end"])
    return merge_df

def jackknife_std(jackknife_df, col_prefix):
    """
    Is this working as expected?
    Require the callable column to filter out noise
    """
    df = jackknife_df.copy()
    jackknife_cols = [cols for cols in df.columns if cols.startswith(col_prefix)]
    n_jackknife = len(jackknife_cols)
    ### Row-wise jackknife and estimate the variance by *n/(n-1)
    mean_jackknife = df[jackknife_cols].mean(axis = 1)
    
    diff_jackknife = df[jackknife_cols].sub(mean_jackknife, axis = 0)
    sq_diff_jackknife = np.square(diff_jackknife)
    
    sum_sq_diff_jackknife = np.sum(sq_diff_jackknife, axis = 1)
    var_jackknife = sum_sq_diff_jackknife/(n_jackknife)*(n_jackknife-1)
    return np.sqrt(var_jackknife)

def sort_jackknifle_cols(jackknife_df, default_cols: list):
    """
    Hard coded for the dataset
    """
    df = jackknife_df.copy()
    cols = list(df.columns)
    ### assert all data_cols are present
    assert set(default_cols).difference(cols) == set()
    sigDz_cols = [col for col in cols if col.startswith("sigDz_jk")]
    sigD2_cols = [col for col in cols if col.startswith("sigD2_jk")]
    sort_cols = default_cols + sigDz_cols + sigD2_cols
    other_cols = [col for col in cols if col not in sort_cols]
    sort_cols = sort_cols + other_cols
    df = df[sort_cols]
    return df

#def jackknife_quantile(jackknife_df)

population_path="/home/alouette/projects/ctb-sgravel/alouette/Dz_sweep_final/simulation_review/src/output/"
jackknife_path = "/home/alouette/projects/ctb-sgravel/alouette/Dz_sweep_final/simulation_review/src/jackknife_output"
callable_path = "/home/alouette/projects/ctb-sgravel/alouette/Dz_sweep_final/blacklist/LD_windows.callable_overlap.merged.bed"
output_path = "/home/alouette/projects/ctb-sgravel/alouette/Dz_sweep_final/simulation_review/src/jackknife_script/jackknife_combined"


populations = ["AMR", "SAS"]
all_pops = "ALL"
ALL_df = []
pop_df_list = []
default_cols = ['chr', 'window_start', 'window_end', 'window_size', 'snp_pairs', 'D2', 'Dz', 'D', 'pi2', 'sigD2', 'sigDz']
callable_df = pd.read_csv(callable_path, sep = "\t")
merge_df_list =[]


### read ALL
for chr in range(1, 23):
    chromosome = "chr" + str(chr)
    all_pops_path = os.path.join(population_path, all_pops + "_chr", ".".join([chromosome, all_pops, "filtered.dict.npy"]))
    ALL_df.append(read_filtered(all_pops_path, chromosome))

ALL_merged = concatenate_filtered(ALL_df)
    
### read populations, jacknife need to be added per each chromosome (avoid listing all the subpops)
for superpop in populations:
    chr_df_list = []
    for chr in range(1, 23):
        chromosome = "chr" + str(chr)
        pop_path = os.path.join(population_path, superpop + "_chr", ".".join([chromosome, superpop, "filtered.dict.npy"]))
        df = read_filtered(pop_path, chromosome)
        ### read jackknife
        pop_jk_path = [ f for f in os.listdir(os.path.join(jackknife_path, superpop)) if f.startswith(".".join([chromosome, superpop, "jackknife"]))]
        pop_jk_path.sort()
        chr_jk_list = []
        for infile in pop_jk_path:
            chr_jk_list.append(read_jackknife(os.path.join(jackknife_path, superpop, infile), chromosome))
        jackknife_df = combine_jackknife(chr_jk_list)
        
        merge_df = df.merge(jackknife_df, on = ["chr", "window_start", "window_end"])
        chr_df_list.append(merge_df)
    ### Add all chromosomes together
    pop_df = pd.concat(chr_df_list, axis = 0)
    ### Sort columns
    pop_df = sort_jackknifle_cols(pop_df, default_cols)
    ### Calculate std from jackknife across all chromosomes
    pop_df["std_jk_sigDz"] = jackknife_std(pop_df, "sigDz_jk")
    pop_df["std_jk_sigD2"] = jackknife_std(pop_df, "sigD2_jk")
    ### Sort by chromosomes
    pop_df["sort_chr"] = pop_df.chr.str.split("chr").str[1].astype(int)
    pop_df = pop_df.sort_values(["sort_chr", "window_start"]).reset_index(drop = True)
    #jackknife_df = jackknife_df[[cols for cols in jackknife_df.columns if not cols.startswith("sigD2_jk")]]
    pop_df_list.append(pop_df)


### Save
for superpop, superpop_df in zip(populations, pop_df_list):
    superpop_df.to_csv(os.path.join(output_path, superpop + "_LD.jackknife_combined.tsv"), sep = "\t", index = False)

#ALL_merged.to_csv(os.path.join(output_path,  "ALL_LD.jackknife_combined.tsv"), sep = "\t", index = False)

### merge all files and add callable regions
### snp pairs might be the same throughout 
#callable_df = pd.read_csv(callable_path, sep = "\t")
#merge_df = ALL_merged.merge(AFR_merged, on = ["chr", "window_start", "window_end", "window_size"], suffixes = ["", "_AFR"])
#merge_df = merge_df.merge(EAS_merged, on = ["chr", "window_start", "window_end", "window_size"], suffixes = ["", "_EAS"])
#merge_df = merge_df.merge(EUR_merged, on = ["chr", "window_start", "window_end", "window_size"], suffixes = ["", "_EUR"])
#merge_df = merge_df.merge(callable_df, on = ["chr", "window_start", "window_end"])

#merge_df.to_csv(os.path.join(output_path, "ALL_superpop_merged_LD.jackknife.tsv"), sep = "\t", index = False)