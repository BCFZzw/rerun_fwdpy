import allel
import numpy as np


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