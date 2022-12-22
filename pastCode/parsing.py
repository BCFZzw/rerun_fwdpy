from turtle import shape
import allel
import numpy as np
import itertools
import pandas as pd
import os 
import fwdpy11
import tskit as ts
    

def segmentation(allele, distanceMatrix, maxLen, boundaries):
    """
    Separate a given allele matrix into equal genomic distance segments
    @param: allele {np.narray}. "calldata/GT" matrix extracted from scikit allele of a vcf
    @param: distanceMatrix. {np.narray} "variants/POS" matrix extracted from scikit allele of a vcf 
    @param: maxLen {int}. the maximum length from a vcf file
    @param: boudnaries {int}. the intended separation segments of the genotypes given allele's genomic region
    @returns segment {np.narray} [shape=(boundaries, _, indNum, 2)], segmentName {np.array} [shape=boundaries]. return the moment accepted segment format and associated names in two arrays. 
    """
    segment = np.empty(dtype = "O", shape = boundaries)
    #segmentName = np.empty(dtype = ("str", len(str(boundaries))+1), shape = boundaries) #the str storage size determined by the length of boundaries, boundaries = 100 => str lenght of 3, +1 for R at the beginning
    for i in range(1, boundaries + 1):
        if i == 1:
            index = distanceMatrix < i*maxLen/boundaries
        elif i == boundaries:
            index = distanceMatrix > (i-1)*maxLen/boundaries
        else:
            index = (distanceMatrix > (i-1)*maxLen/boundaries) & (distanceMatrix < i*maxLen/boundaries)
        
        segment[i-1] = (allele[index, :])
        #segmentName[i-1] = "R"+str(i)
    return segment#, segmentName


def splitData(npyFile, straddle = False, name = False):
    pdList = [pd.DataFrame.from_records(npyFile[i]) for i in range(0 , len(npyFile))]
    D2 = pd.DataFrame([entry["D2"] for entry in pdList]).reset_index(drop = True)
    Dz = pd.DataFrame([entry["Dz"] for entry in pdList]).reset_index(drop = True)
    Pi2 = []
    if not straddle:
        Pi2 = pd.DataFrame([entry["Pi2"] for entry in pdList]).reset_index(drop = True)
    D = pd.DataFrame([entry["D"] for entry in pdList]).reset_index(drop = True)
    if name:
        name = pd.DataFrame([entry["name"] for entry in pdList]).reset_index(drop = True)
        name = name.iloc[0] # get the first row of region names, all the rest of the data are the same
    else:
        name = ""
    return D2, Dz, Pi2, D, name

    
def shapeTransform(allel):
    """
    Return a genotype array suitable for moment parsing
    Add within each bracket, haplotype 0|0 = 0, 0|1 = 1, 1|1 = 2
    @param allel {list/numpy.narray}: haplotype matrix from allel 
    @return {numpy.narray}: 0, 1, 2 genotype matrix for moment parsing
    """
    return np.sum(allel, axis =2)


def pairedSegments(segment, segmentName, numPaired = 2):
    """
    Used to compute between region stastics for omega statistics. 
    @update: Include same region comparison.
    @oaram segment {np.narray}: list of genotype array of a given vcf
    @param segmentName {np.array}: list of name associateed with regions
    @param numPaired {int}: default 2 for moments' input
    @returns {list} return list of paired segments and associated segment name in the same order
    """
    pairedSeg = list(itertools.combinations_with_replacement(segment, numPaired))
    pairedSegName = list(itertools.combinations_with_replacement(segmentName, numPaired))
    return pairedSeg, pairedSegName


def readVcfs(vcfFileList, inputPath):
    """
    read vcf from list provided (with input directory or not, if not inputPath = "")
    @param vcfFileList {list of strings}
    @param inputPath {string}
    @return {np.narray}
    """
    legnth = len(vcfFileList)
    genotype = np.empty(dtype = "O", shape = legnth)
    distance = np.empty(dtype = "O", shape = legnth)
    for i in range(legnth):
        distance[i], genotype[i],  = readVcf(vcfFileList[i], inputPath)
    
    return distance, genotype

def readVcf(fileName, inputPath): #checked
    """
    read a single vcf (with input directory or not, if not inputPath = "")
    @param fileName {string .trees}
    @param inputPath {string}
    @return {np.narray two fields}
    @TODO: extract contig length, slightly hard and might require hardcoding
    """
    read = allel.read_vcf(inputPath + fileName, fields = ["variants/POS", "calldata/GT"])
    distance = read["variants/POS"]
    momentGenotype = shapeTransform(read["calldata/GT"])
    return distance, momentGenotype


    
def treeData(simSavePath, nsim, nbins, nsamples, mode = "site"):
    #@TODO: the data structure needs to be re-considered
    treeList = getFileList(simSavePath, ".trees")
    tree = ts.load(simSavePath + treeList[0])
    pop = fwdpy11.DiploidPopulation.create_from_tskit(tree)
    windows = np.linspace(0, pop.tables.genome_length, nbins) #equal spacing between distances

    treeData = np.empty(dtype = 'O', shape = (nsim,nbins-1))
    for i in range(0, nsim):
        tree = ts.load(simSavePath + treeList[i])
        pop = fwdpy11.DiploidPopulation.create_from_tskit(tree)
        
        treeData[i, :] = tree.diversity(
            sample_sets=[pop.alive_nodes[:nsamples]], windows=windows, mode=mode
        ).flatten()
    
    return windows, treeData


def getFileList(simSavePath, suffix):
    """
    List and sort all files inside a directory, note that it can include undesired files
    sorting is done naively
    @param: File directory {string}
    @return: {list}
    """
    fileList  = os.listdir(simSavePath)
    fileList.sort()
    return [f for f in fileList if f.endswith(suffix)]

def getSimLength(simSavePath):
    treeList = getFileList(simSavePath, ".trees")
    assert len(treeList) > 0
    return ts.load(simSavePath + treeList[0]).sequence_length

def combineFiles(statsSavePath):
    """
    Concatenate all npy files within a folder into one npy stats
    @return: np.ndarray{dict}
    """
    stats = np.array([])
    fileList = getFileList(statsSavePath, ".npy")
    for files in fileList:
        stats = np.concatenate((stats, (np.load(statsSavePath + files, allow_pickle = True))))
    return stats
    

