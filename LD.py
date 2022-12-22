import allel
import zarrMoment1000G as zm
import sys
import os
import numpy as np

def simPipeline(callset):
    genotype = allel.GenotypeChunkedArray(callset['calldata/GT'])
    #daskGenotype = allel.GenotypeDaskArray(genotype)
    position = allel.SortedIndex(callset['variants/POS'])
    assert np.shape(genotype)[0] == np.shape(position)[0]

    binEdgeBp = [i*4e5 for i in range(0, 26)]
    
    df = zm.pandasParsing(position, binEdgeBp)
    print("Bins created: "+ str(df.labels.nunique()))
    print("Loading complete")

    geno012 = genotype.to_n_alt(fill=-1)

    statsDf = zm.parseAndMoment(df, geno012, "", "", False)

    return statsDf


vcfPath = "/home/alouette/projects/ctb-sgravel/alouette/fwdpy_data/vcfs/Nielsen_1.0"
savePath = "/home/alouette/projects/ctb-sgravel/alouette/rerun_fwdpy_data/previous_run_LD"
job = int(sys.argv[1])

for i in range((job-1) * 5, job*5):
    vcfID = str(i)

    sweepVcf = os.path.join(vcfPath, "sweep", "[" + vcfID + "]sweep.vcf")
    neutralVcf = os.path.join(vcfPath, "no_sweep", "[" + vcfID + "]no_sweep.vcf")


    callset = allel.read_vcf(sweepVcf)

    statsDf = simPipeline(callset)

    statsDf.to_csv(os.path.join(savePath, "sweep", "sweep_" + vcfID  + ".csv"), index = False)

    callset = allel.read_vcf(neutralVcf)

    statsDf = simPipeline(callset)

    statsDf.to_csv(os.path.join(savePath, "no_sweep", "neutral_" + vcfID  + ".csv"), index = False)


