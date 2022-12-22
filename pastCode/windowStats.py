import numpy as np
import moments.LD
import allel
import os
import pandas as pd

def pandasParsing(pos, maxLen, binSize):
    """
    Updated to the correct one, still using pandas to subset quickly
    """
    binEdge = np.array([i * maxLen/binSize for i in range(0, binSize + 1)])

    df = pd.DataFrame(pos, columns =  ["positions"])
    #@TODO: kept non polymorphic but skip with moments
    labels = pd.cut(pos, binEdge, right = False)
    df["labels"] = labels
    indicesList = list(df.groupby("labels").indices.values()) #need to add list to avoid non-accessible issue

    return df, indicesList


def correctSegment(genotype, indiceList, index):
    indice = indiceList[index]
    if len(indice) == 0:
        return None

    return np.take(genotype, indice, axis = 0)


def segmentation(allele, distanceMatrix, maxLen, boundaries):
    """
    Wrong, deprecated
    -------------------------
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

def shapeTransform(allel):
    """
    Return a genotype array suitable for moment parsing
    -------previous function uses np.sum(axis = 2), updated to scikit-allel function

    """
    return allel.to_n_alt(fill=-1)

def readVcf(fileName, inputPath): #checked
    """
    read a single vcf (with input directory or not, if not inputPath = "")
    @param fileName {string .trees}
    @param inputPath {string}
    @return {np.narray two fields}
    @TODO: extract contig length, slightly hard and might require hardcoding
    """
    read = allel.read_vcf(os.path.join(inputPath ,fileName), fields = ["variants/POS", "calldata/GT"])
    distance = read["variants/POS"]
    momentGenotype = shapeTransform(read["calldata/GT"])
    return distance, momentGenotype
            
def binStats(segment):
    """
    For one vcf, use given boundaries to separate alleles based on equal distance into bins, achieved via parsing.segmentation(). Further used in allBinStats.
    @TODO: make a class for statistics? that way segments can be inherited 
    """
    size = len(segment)
    D2 = np.empty(shape = size)
    Dz = np.empty(shape = size)
    Pi2 = np.empty(shape = size)
    D = np.empty(shape = size)
    #name = np.empty(shape = size, dtype = ("str", 10))
    for i in range(0, size):
        #TODO change to parsing
        region = segment #originally multiple segments together, remove [i]
        stats = moments.LD.Parsing.compute_pairwise_stats(region, genotypes = True)
        
        D2[i] = np.sum(stats[0])
        Dz[i] = np.sum(stats[1])
        Pi2[i] = np.sum(stats[2])
        D[i] = np.sum(stats[3])
    #    name[i] = segmentName[i]

    resultDict = {"D2":D2, 
                  "Dz":Dz, 
                  "Pi2":Pi2, 
                  "D":D,}
                 #"name":name}
    return resultDict
        
def allBinStats(vcf_list, inputPath, outputPath, maxLen = 1e7, outName = "windowStats.npy",  boundaries = 25, save = True):
    """
    Process binStats for all vcfs. Save format may be improved?
    """
    length = len(vcf_list)
    resultList = np.empty(dtype = "O", shape = length)
    
    for i in range(0, length):
        distance, genotype = readVcf(vcf_list[i], inputPath)
        segment = segmentation(genotype, distance, maxLen, boundaries)
        result =  binStats(segment)
        resultList[i] = result
    

    if save:
        np.save(outputPath + outName, resultList)

    return resultList


def pipeline(inputVcf, savePath, maxLen = 1e7, binSize = 25, save = True):
    callset = allel.read_vcf(inputVcf, fields = ["variants/POS", "calldata/GT"])
    pos = callset["variants/POS"]
    genotype = callset["calldata/GT"]
    
    df, indiceList = pandasParsing(pos, maxLen, binSize)
    assert len(indiceList) == binSize
    for i in range(binSize):
        segment = correctSegment(genotype, indiceList, i)
        result =  binStats(segment)
