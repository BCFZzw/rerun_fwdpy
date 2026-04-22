import allel
import numpy as np
import moments.LD

def getFileList(inputPath, suffix):
    """
    List and sort all files inside a directory, note that it can include undesired files
    sorting is done naively
    @param: File directory {string}
    @return: {list}
    """
    fileList  = os.listdir(inputPath)
    fileList.sort()
    return [f for f in fileList if f.endswith(suffix)]



def tskit_to_allel(tree, genotype = True):
    data = allel.HaplotypeArray(tree.genotype_matrix())
    if genotype:
        data = data.to_genotypes(ploidy = 2)
    pos_array = allel.SortedIndex(tree.sites_position)
    return data, pos_array

def locate_window_region(pos_array, win_start = None, win_end = None):
    """
    Locate the indices of a region in the genotype position array.
    Improvement can be added when there are SNPs on the edge of the windows.
    If no win_start or win_end is given, select the first or last position of the pos array.
    If there is no position within the win_start or win_end, return an empty slice(None, 0, None).
    """
    if win_start is None:
        win_start = 1
    if win_end is None:
        win_end = max(pos_array) + 1
    assert win_end > win_start
    ### no SNPs in region
    ### exclusive at both sides at the moment, better to make 1 inclusive
    ### make sure there is no double counting of the SNPs per windows
    if np.all((pos_array > win_start) & (pos_array < win_end) == False):
        return slice(None, 0, None)
    else:
        loc_region = pos_array.locate_range(win_start, win_end)
        return loc_region


def select_variants(genotype, pos_array, pos_start, pos_end):
    """
    Subset genotype based on given positions.
    """
    loc_region = locate_window_region(pos_array, pos_start, pos_end)
    pos_array_var = pos_array[loc_region]
    ### when using indices: take, when using boolean: compress
    genotype_var = genotype[loc_region]
    return genotype_var, pos_array_var


def windowed_diversity(genotype: allel.GenotypeArray, pos_array: np.array, window_list: list):
    ac = genotype.count_alleles()
    pi, windows, n_bases, counts = allel.windowed_diversity(pos_array, ac, windows = window_list)
    return pi


def norm_iHS(haplotype: allel.HaplotypeArray, pos_array: np.array):
    ### iHS is calculated for each variant, and then normalized by allele count
    score = allel.ihs(haplotype, pos_array)
    ac = haplotype.count_alleles()
    aac = ac[:, 1]
    #norm_score, bins = allel.standardize_by_allele_count(score, aac)
    return score


def windowed_tajima_D(genotype: allel.GenotypeArray, pos_array: np.array, window_list: list):
    ac = genotype.count_alleles()
    tajima_D_list, windows, counts = allel.windowed_tajima_d(pos_array, ac, windows = window_list)
    return tajima_D_list


def window_roh(genotype: allel.GenotypeArray, pos_array: np.array, contig_size):
    """
    Runs of homozygosity in fixed sized window. Calculated by poisson HMM.
    Fast but reduced resolution.
    """
    genotype_vector = allel.GenotypeVector(genotype)
    df_roh, fraction_roh = allel.roh_poissonhmm(genotype, pos_array, contig_size = contig_size)
    return df_roh, fraction_roh 


def window_LD(genotype: allel.GenotypeArray, pos_array: np.array, window_list: list):
    ### scikit allel window utilities cannot handle this
    stats = {"D2": [], "Dz": [], "pi2": []}
    n_wins = len(window_list)
    for i in range(n_wins):
        genotype_win, pos_win = select_variants(genotype, pos_array, window_list[i][0], window_list[i][1])
        genotype_012 = genotype_win.to_n_alt(fill = -1)
        ### threshold and filtering like the real data?
        D2_pw, Dz_pw, pi2_pw, D_pw = moments.LD.Parsing.compute_pairwise_stats(genotype_012, genotypes = True)
        stats["D2"].append(np.sum(D2_pw))
        stats["Dz"].append(np.sum(Dz_pw))
        stats["pi2"].append(np.sum(pi2_pw))
    return stats["D2"], stats["Dz"], stats["pi2"]



#def pop_branching_stats(genotype1, genotype2, genotype3, )

#def Fst