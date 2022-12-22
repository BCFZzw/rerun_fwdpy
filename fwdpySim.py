import fwdpy11
import numpy as np
import msprime
import fwdpy11.conditional_models
import fwdpy11.tskit_tools
import sys
import os

nPop = 500 #ancestral size
nSample = 500
NielsenRecRate = 500 
NielsenMutRate =0.002 #4Nmu, mu is the rate of mutation per bp per generation (bp)
NielsenAlpha = 500
simGenomeLength = int(1e7)
addMutLeft = simGenomeLength/2 
addMutRight = simGenomeLength/2+1 
fwdpySimLen = 100
savePath = "/home/alouette/projects/ctb-sgravel/alouette/rerun_fwdpy_data"
job = int(sys.argv[1])


fwdpyAlpha = NielsenAlpha #2Ns
fwdpyMuRate = NielsenMutRate/4/nPop*simGenomeLength # fwdpy takes haploid genome mutation rate per generation
fwdpyRecRate = NielsenRecRate/4/nPop #genome wide recombination rate
msprimeRecRate = NielsenRecRate/4/nPop/simGenomeLength




def setup(seed, simLen = fwdpySimLen):
    # Dropping mutations requires existing
    # ancestry, which we can get either
    # from a burn-in or from msprime.
    initial_ts = msprime.sim_ancestry(
        samples=nSample,
        population_size=nPop,
        recombination_rate=msprimeRecRate,
        random_seed=seed+1,
        sequence_length=simGenomeLength
    )

    # Build the pop from msprime output
    pop = fwdpy11.DiploidPopulation.create_from_tskit(initial_ts)

    # Set up basic model parameters
    pdict = {
        "recregions": [fwdpy11.PoissonInterval(0, simGenomeLength, fwdpyRecRate, discrete=True)],
        # Here, 2 means that fitness is multiplicative
        # over 1, 1+hs, 1+2s.
        "gvalue": fwdpy11.Multiplicative(2.0),
        "rates": (0, 0, None),
        "prune_selected": False,
        "simlen": simLen,
    }
    params = fwdpy11.ModelParams(**pdict)

    return pop, params

for i in range((job-1) * 5, job*5):
    seed = 100 + i
    rng = fwdpy11.GSLrng(seed)
    pop, params = setup(seed)

    mutation_data = fwdpy11.conditional_models.NewMutationParameters(
        frequency=fwdpy11.conditional_models.AlleleCount(1),
        data=fwdpy11.NewMutationData(effect_size=fwdpyAlpha / 2 / pop.N, dominance=1),
        position=fwdpy11.conditional_models.PositionRange(left=addMutLeft, right=addMutRight),
    )

    output = fwdpy11.conditional_models.selective_sweep(
        rng, 
        pop,
        params,
        mutation_data,
        fwdpy11.conditional_models.GlobalFixation()
    )

    assert output.pop.generation == params.simlen # miss leading
    assert pop.generation == 0

    FIXATION_TIME = output.pop.fixation_times[0]

    nmuts = fwdpy11.infinite_sites(rng, output.pop, fwdpyMuRate)

    print(f"{nmuts} mutations added to sweep")
    
    ts = output.pop.dump_tables_to_tskit()

    ts.dump(os.path.join( savePath, "sim", "sweep", "sweep_"  + str(i) +  ".trees"))
    with open(os.path.join( savePath, "vcf", "sweep", "sweep_"  + str(i) +  ".vcf"), "w") as vcfFile:
        ts.write_vcf(vcfFile, individuals =  [i for i in range(50)])

    neuPop, neuParams = setup(seed, FIXATION_TIME)

    fwdpy11.evolvets(rng, neuPop, neuParams, simplification_interval = 100)

    nmuts = fwdpy11.infinite_sites(rng, neuPop, fwdpyMuRate)
    print(f"{nmuts} mutations added to neutral")

    ts = neuPop.dump_tables_to_tskit()
    with open(os.path.join( savePath, "vcf", "neutral", "neutral_"  + str(i) +  ".vcf"), "w") as vcfFile:
        ts.write_vcf(vcfFile, individuals =  [i for i in range(50)])

    #ts.dump(os.path.join( savePath, "neutral", "neutral_"  + str(i) +  ".trees"))