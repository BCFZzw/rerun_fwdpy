import numpy as np
import moments.LD
import parsing as ps


            
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
        #region = ps.shapeTransform(segment[i])
        region = segment[i]
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
    #distanceList, genotypeList = ps.readVcfs(vcf_list, inputPath)
    resultList = np.empty(dtype = "O", shape = length)
    
    for i in range(0, length):
        distance, genotype = ps.readVcf(vcf_list[i], inputPath)
        segment, segmentName = ps.segmentation(genotype, distance, maxLen, boundaries)
        result =  binStats(segment, segmentName)
        resultList[i] = result
    

    if save:
        np.save(outputPath + outName, resultList)

    return resultList


    
