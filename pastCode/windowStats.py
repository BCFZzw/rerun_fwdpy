import numpy as np
import moments.LD
import parsing as ps
import tracemalloc


def dictLdScore(LdInput):
    """
    LD results from moments.LD, [D2, Dz, Pi2, D], inidiviudal variant pair level
    sum of all 50 individuals' LD scores within one simulations, assign to a dictionary
    """
    #assert input should be length 4
    assert len(LdInput) == 4
    LdDict = {
        "D2": np.sum(LdInput[0]),
        "Dz": np.sum(LdInput[1]),
        "Pi2": np.sum(LdInput[2]),
        "D": np.sum(LdInput[3])}
    return LdDict

            
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
        
def allBinStats(vcf_list, inputPath, outputPath, maxLen = 1e7, outName = "windowStats.npy",  boundaries = 25, straddling = False, save = True):
    """
    Process binStats for all vcfs. Save format may be improved?
    """
    length = len(vcf_list)
    #distanceList, genotypeList = ps.readVcfs(vcf_list, inputPath)
    resultList = np.empty(dtype = "O", shape = length)
    straddleList = np.empty(dtype = "O", shape = length)
    
    for i in range(0, length):
        distance, genotype = ps.readVcf(vcf_list[i], inputPath)
        #genotype = genotypeList[i]
        #distance = distanceList[i]
        segment, segmentName = ps.segmentation(genotype, distance, maxLen, boundaries)
        if straddling:
            straddle = binStraddling(segment, segmentName)
            straddleList[i] = straddle
        else:
            result =  binStats(segment, segmentName)
            resultList[i] = result
    

    
    if straddling and save:
        np.save(outputPath + outName, straddleList)
    elif save:
        np.save(outputPath + outName, resultList)

    return resultList, straddleList


    
#no longer important anymore
def binStraddling(segment, segmentName):
    """"
    Deprecated
    For one vcf, test the straddling region
    @TODO: the regions will need to adjust to more custome range, and not just equal splitting, one way to do so is to feed in segment and segmentName without passing all parameters
    """
    
    pairedSeg, pairedSegName = ps.pairedSegments(segment, segmentName, numPaired=2)
    size = len(pairedSeg)
    D2 = np.empty(shape = size)
    Dz = np.empty(shape = size)
    #Pi2 = np.empty(shape = size)
    D = np.empty(shape = size)
    name = np.empty(shape = size, dtype = ("str", 10)) #parsing funciton has a soft method for adjusting dimensions
    tracemalloc.start()
    for i in range(0, size):
        #TODO: change to parsing
        #region1 = ps.shapeTransform(pairedSeg[i][0])
        #region2 = ps.shapeTransform(pairedSeg[i][1])
        region1 = pairedSeg[i][0]
        region2 = pairedSeg[i][1]

        stats = moments.LD.Parsing.compute_pairwise_stats_between(region1, region2, genotypes = True)
        print(tracemalloc.get_traced_memory())
        #for omega statistics, average of ratios instead of ratios of averages
        #EXCLUDE zero Pi2, and respective statistics, this steps takes about 1 min
        nonzeroIndex = np.argwhere(stats[2] != 0)
        nonzeroD2 = stats[0][nonzeroIndex]
        nonzeroDz = stats[1][nonzeroIndex]
        nonzeroPi2 = stats[2][nonzeroIndex]
        nonzeroD = stats[3][nonzeroIndex] 

        D2[i] = np.sum(nonzeroD2/nonzeroPi2)
        Dz[i] = np.sum(nonzeroDz/nonzeroPi2)
        D[i] = np.sum(nonzeroD/nonzeroPi2)
        name[i] = "".join(pairedSegName[i])

        assert (not np.isnan(D2[i])) and (not np.isinf(D2[i]))
        assert (not np.isnan(Dz[i])) and (not np.isinf(Dz[i]))
        assert (not np.isnan(D[i])) and (not np.isinf(D[i]))
        
    
    tracemalloc.stop() #takes around 50G

    #for one vcf file, return a dictionry of this format
    resultDict = {"D2":D2, 
                  "Dz":Dz, 
                  "D":D,
                  "name": name}
    return resultDict

    