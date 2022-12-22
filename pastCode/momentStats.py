import windowStats as ws
import sys
import numpy as np
import os
import parsing as ps

start = int(sys.argv[1])
end = int(sys.argv[2])
condition= sys.argv[3]
genomeLength = int(sys.argv[4])*1e7

#read in from a file listing all file names

fileNamePath = "/home/alouette/projects/ctb-sgravel/alouette/fwdpy_data/prunedVcfs/"
vcfPath = "/home/alouette/projects/ctb-sgravel/alouette/fwdpy_data/prunedVcfs/" + condition + "/sweep/"
NSvcfPath = "/home/alouette/projects/ctb-sgravel/alouette/fwdpy_data/prunedVcfs/" + condition + "/no_sweep/"
#TODO: better ways to auto adjust this



SWvcfList = ps.getFileList(vcfPath, ".vcf")[start:end]
NSvcfList = ps.getFileList(NSvcfPath, ".vcf")[start:end]

statsSWPath = "/home/alouette/projects/ctb-sgravel/alouette/fwdpy_data/statistics/" + condition + "/sweep/"
statsNSPath = "/home/alouette/projects/ctb-sgravel/alouette/fwdpy_data/statistics/" + condition + "/no_sweep/"

assert os.path.isdir(statsSWPath)
assert os.path.isdir(statsNSPath)


nbins = 25 * int(sys.argv[4])
fileName = str(nbins) + " Window Stats.npy"



ws.allBinStats(SWvcfList, vcfPath, statsSWPath, maxLen = genomeLength, outName = "[" + str(start) + "-" + str(end) + "] sweep " + fileName, straddling = False, boundaries = nbins, save = True)
ws.allBinStats(NSvcfList, NSvcfPath, statsNSPath, maxLen = genomeLength, outName = "[" + str(start) + "-" + str(end) + "] non sweep " + fileName, straddling = False, boundaries = nbins, save = True)


