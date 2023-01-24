import allel
import sys
import os
import numpy as np
import pandas as pd
import moments

job = int(sys.argv[1])
vcfPath = sys.argv[2]
savePath = sys.argv[3]


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

def pandasParsing(subsetPos, binEdgeBp):
    """
    Achieve segmentation via pandas's groupby and cut functions
    Bin range, position data is saved to csv meta file
    Also utilizes a few function from scikit allel for subsetting
    """

    df = pd.DataFrame(subsetPos, columns =  ["positions"])
    #@TODO: kept non polymorphic but skip with moments
    labels = pd.cut(subsetPos, binEdgeBp, right = False)
    df["labels"] = labels
    groupDf = df.groupby("labels")
    df["counts"] = groupDf.transform("count")

    return df

def parseAndMoment(df, subsetGeno, filePrefix = "", savePath = ""):
    """
    Achieve segmentation via pandas's groupby and cut functions
    Bin range, position data is saved to csv meta file
    @TODO: improve performance by unsorting groupby
    """

    size = len(df)
    D2 = np.empty(shape = size)
    Dz = np.empty(shape = size)
    Pi2 = np.empty(shape = size)
    D = np.empty(shape = size)
    medianPos = np.empty(shape = size)
    snpCounts = np.empty(shape = size)
    
    indicesList = list(df.groupby("labels").indices.values() )
    counter = 0

    for indices in indicesList: #@TODO possibly not optimal
        if len(indices) < 2:
            #@TODO: save as 0 or sth else, same with polymorphisms
            continue
        
        data = np.take(subsetGeno, indices, axis = 0) #optimize?
        assert np.shape(data)[0] == len(indices)
        stats = moments.LD.Parsing.compute_pairwise_stats(data, genotypes = True)
        # when genotypes = True, use 0, 1, 2
        #I remember this is faster than df[idx] = ... each time 
        D2[counter] = np.sum(stats[0]) #sum up every possible pair values
        Dz[counter] = np.sum(stats[1])
        Pi2[counter] = np.sum(stats[2])
        D[counter] = np.sum(stats[3])
        medianPos[counter] = np.median(df["positions"][indices]).astype(int)
        snpCounts[counter] = df["counts"][indices[0]]
        counter = counter + 1
    
    D2 = D2[:counter]
    Dz = Dz[:counter]
    Pi2 = Pi2[:counter]
    D = D[:counter]
    medianPos = medianPos[:counter]
    snpCounts = snpCounts[:counter]

    statsDf = pd.DataFrame({"medPos": medianPos,
                            "snpCounts": snpCounts,
                            "D2":D2,
                            "Dz": Dz,
                            "Pi2": Pi2,
                            "D": D
                            })


    return statsDf


def simPipeline(callset, maxLen, binSize):
    genotype = allel.GenotypeChunkedArray(callset['calldata/GT'])
    #daskGenotype = allel.GenotypeDaskArray(genotype)
    position = allel.SortedIndex(callset['variants/POS'])
    assert np.shape(genotype)[0] == np.shape(position)[0]

    binEdgeBp = [i* maxLen/binSize for i in range(0, binSize + 1)]
    
    df = pandasParsing(position, binEdgeBp)
    print("Bins created: "+ str(df.labels.nunique()))
    print("Loading complete")

    geno012 = genotype.to_n_alt(fill=-1)

    print("conversion complete")

    statsDf = parseAndMoment(df, geno012, "", "")

    return statsDf





genomeLength = 1e7
binSize = 25

sweepList  = getFileList(os.path.join(vcfPath, "sweep"), ".vcf")
neutralList  = getFileList(os.path.join(vcfPath, "neutral"), ".vcf")


for i in range((job-1) * 5, job*5):
    sweepFile = sweepList[i]
    neutralFile = neutralList[i]

    sweepVcf = os.path.join(vcfPath, "sweep", sweepFile)
    neutralVcf = os.path.join(vcfPath, "neutral", neutralFile)


    callset = allel.read_vcf(sweepVcf)

    statsDf = simPipeline(callset, genomeLength, binSize)
    break
    statsDf.to_csv(os.path.join(savePath, "sweep", "sweep_" + vcfID  + ".csv"), index = False)

    
    callset = allel.read_vcf(neutralVcf)

    statsDf = simPipeline(callset, genomeLength, binSize)
    
    statsDf.to_csv(os.path.join(savePath, "neutral", "neutral_" + vcfID  + ".csv"), index = False)


