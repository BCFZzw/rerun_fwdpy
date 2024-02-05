import allel
import sys
import os
import numpy as np
import pandas as pd
import moments
import argparse

parser = argparse.ArgumentParser(
                    prog='Sweep-neutral-simulation',
                    description='This program simulates for a hard sweep and generate respective neutral simualtions'
                   )

parser.add_argument('-p', '--parallel', dest = "thread", help = "Add parallel option", default = 0)
parser.add_argument('-j', '--job-id', dest = "jobID", help = "which job to run", default = 0)
args = parser.parse_args()
thread = int(args.thread)
job_id = args.jobID

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


#sweep_vcf = "/home/alouette/projects/ctb-sgravel/alouette/simulation.0.21.0/VCF/sweep/sweep_5563.vcf"
#neutral_vcf = "/home/alouette/projects/ctb-sgravel/alouette/simulation.0.21.0/VCF/neutral/neutral_5563.vcf"
inputPath= "/home/alouette/projects/ctb-sgravel/alouette/Dz_sweep_simulation/VCF/sweep"
savePath = "/home/alouette/projects/ctb-sgravel/alouette/Dz_sweep_simulation/moments/sweep"
input_file_list = getFileList(inputPath, ".vcf")
input_file_job = np.array_split(input_file_list, thread)[int(job_id)]

#base = sweep_vcf.split("/")[-1].strip(".vcf")
binSize = 25
genome_length = 1e7
binEdge = [i* genome_length/binSize for i in range(0, binSize + 1)]
bin_index = list(range(len(binEdge)-1))

D2_dict = {}
Dz_dict = {}
pi2_dict = {}
D_dict = {}
for i in bin_index:
    D2_dict[str(i)] = []
    Dz_dict[str(i)] = []
    pi2_dict[str(i)] = []
    D_dict[str(i)] = []

for vcf in input_file_job:
    callset = allel.read_vcf(os.path.join(inputPath, vcf))
    genotype = allel.GenotypeChunkedArray(callset['calldata/GT'])
    position = allel.SortedIndex(callset['variants/POS'])
    
    for i in bin_index:
        bin_start = binEdge[i]
        bin_end = binEdge[i+1]
        loc_region = position.locate_range(bin_start, bin_end)
        gt_region = allel.GenotypeArray(genotype[loc_region])
        gt_region_012 = gt_region.to_n_alt(fill=-1)
        ### genotypes encoded as 0, 1, 2
        D2_pw, Dz_pw, pi2_pw, D_pw = moments.LD.Parsing.compute_pairwise_stats(gt_region_012, genotypes = True)
        D2_dict[str(i)].append(np.sum(D2_pw))
        Dz_dict[str(i)].append(np.sum(Dz_pw))
        pi2_dict[str(i)].append(np.sum(pi2_pw))
        D_dict[str(i)].append(np.sum(D_pw))
    print("finished " + vcf)

pd.DataFrame(D2_dict).to_csv(os.path.join(savePath, "D2.simulation.bin." + job_id + ".csv.gz"), index = False, sep ="\t", compression = "gzip")
pd.DataFrame(Dz_dict).to_csv(os.path.join(savePath, "Dz.simulation.bin." + job_id + ".csv.gz"), index = False, sep ="\t", compression = "gzip")
pd.DataFrame(pi2_dict).to_csv(os.path.join(savePath, "pi2.simulation.bin." + job_id + ".csv.gz"), index = False, sep ="\t", compression = "gzip")
pd.DataFrame(D_dict).to_csv(os.path.join(savePath, "D.simulation.bin." + job_id + ".csv.gz"), index = False, sep ="\t", compression = "gzip")
