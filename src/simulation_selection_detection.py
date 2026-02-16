import allel
import numpy as np

def tskit_to_allel(tree):
    haplotype = allel.HaplotypeArray(tree.genotype_matrix())
    genotype = haplotype.to_genotypes(ploidy = 2)
    pos_array = allel.SortedIndex(tree.sites_position)
    return genotype, pos_array


def locate_window_region(pos_array, win_start = None, win_end = None):
    """
    Locate the indices of a region in the genotype position array.
    Improvement can be added when there are SNPs on the edge of the windows.
    If no win_start or win_end is given, select the first or last position of the pos array.
    If there is no position within the win_start or win_end, return an empty slice(None, 0, None).
    """
    if pos_start is None:
        pos_start = 1
    if pos_end is None:
        pos_end = max(pos_array) + 1
    assert pos_end > pos_start
    ### no SNPs in region
    ### exclusive at both sides at the moment, better to make 1 inclusive
    ### make sure there is no double counting of the SNPs per windows
    if np.all((pos_array > pos_start) & (pos_array < pos_end) == False):
        return slice(None, 0, None)
    else:
        loc_region = pos_array.locate_range(pos_start, pos_end)
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


def genotype_allele_counts(genotype: allel.GenotypeArray):
    return genotype.count_alleles()


def windowed_tajima_D(genotype: allel.GenotypeArray, pos_array: np.array, window_list: list):
    genotype_ac = genotype_allele_counts(genotype)
    tajima_D_list, windows, counts = allel.windowed_tajima_d(pos_array, genotype_ac, windows = window_list)
    return tajima_D_list, windows, counts

def window_roh(genotype: allel.GenotypeArray, pos_array: np.array, contig_size):
    """
    Runs of homozygosity in fixed sized window. Calculated by poisson HMM.
    Fast but reduced resolution.
    """
    genotype_vector = allel.GenotypeVector(genotype)
    df_roh, fraction_roh = allel.roh_poissonhmm(genotype, pos_array, contig_size = contig_size)
    return df_roh, fraction_roh 

#def pop_branching_stats(genotype1, genotype2, genotype3, )

#def Fst