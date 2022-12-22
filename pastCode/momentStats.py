import windowStats as ws
import sys
import numpy as np
import os

start = int(sys.argv[1])
end = int(sys.argv[2])
genomeLength = int(1e7)

#read in from a file listing all file names

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

filePath = "/home/alouette/projects/ctb-sgravel/alouette/rerun_fwdpy_data/vcf"
vcfPath = os.path.join(filePath, "sweep")
NSvcfPath = os.path.join(filePath, "neutral")



SWvcfList = getFileList(vcfPath, ".vcf")[start:end]
NSvcfList = getFileList(NSvcfPath, ".vcf")[start:end]

statsPath = "/home/alouette/projects/ctb-sgravel/alouette/rerun_fwdpy_data/new_vcf_old_moment_LD"

statsSWPath = os.path.join(statsPath, "sweep")
statsNSPath = os.path.join(statsPath, "neutral")

assert os.path.isdir(statsSWPath)
assert os.path.isdir(statsNSPath)


nbins = 25
fileName = str(nbins) + " Window Stats.npy"


ws.allBinStats(SWvcfList, vcfPath, statsSWPath, maxLen = genomeLength, outName = "[" + str(start) + "-" + str(end) + "] sweep " + fileName, boundaries = nbins, save = True)
ws.allBinStats(NSvcfList, NSvcfPath, statsNSPath, maxLen = genomeLength, outName = "[" + str(start) + "-" + str(end) + "] non sweep " + fileName, boundaries = nbins, save = True)


