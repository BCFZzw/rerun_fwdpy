import yaml
import fwdpy11 
import numpy as np 
import fwdpy11.tskit_tools
import fwdpy11.conditional_models
import sys
import os
import msprime 
import argparse
from pathlib import Path

### adapted from fwdpy11 manual on selective sweep, https://molpopgen.github.io/fwdpy11/short_vignettes/incomplete_sweep.html

parser = argparse.ArgumentParser(
                    prog='Sweep-neutral-simulation',
                    description='This program simulates for a hard sweep and generate respective neutral simualtions'
                   )

parser.add_argument('-p', '--parallel', dest = "parallel_seed", help = "Add parallel option", default = 0)
parser.add_argument('-y', '--yaml_param', dest = "yaml_param", help = "Parameter set in the YAML format for sweep simulations. If not given, feed parameters using options.", default = 0)
parser.add_argument('-s', '--save_path', dest = "save_path", help = "Savepath for all the outputs.")
args = parser.parse_args()
parallel = args.parallel_seed
yaml_input = args.yaml_param
savePath = args.save_path

### making directories, python >=3.5
Path(savePath).mkdir(parents=True, exist_ok=True)
Path(os.path.join(savePath, "neutral", "Trees")).mkdir(parents=True, exist_ok=True)
Path(os.path.join(savePath, "neutral", "VCF")).mkdir(parents=True, exist_ok=True)
Path(os.path.join(savePath, "sweep", "Trees")).mkdir(parents=True, exist_ok=True)
Path(os.path.join(savePath, "sweep", "VCF")).mkdir(parents=True, exist_ok=True)

param_yaml = yaml.safe_load(open(yaml_input))
#/home/alouette/projects/ctb-sgravel/alouette/Dz_sweep_final/simulation_review/parameters.yaml"
#savePath = "/home/alouette/projects/ctb-sgravel/alouette/Dz_sweep_simulation"

rec_rate = param_yaml["simulationParam"]["recombinationRate"]
mut_rate = param_yaml["simulationParam"]["mutationRate"]
alpha = param_yaml["simulationParam"]["alpha"]
sim_length = param_yaml["simulationParam"]["genomeLength"] # will get loaded as int
sim_sample = param_yaml["simulationParam"]["nSample"] # will get loaded as int

sim_time = param_yaml["userParam"]["simulationTime"]
time_post_fixation = param_yaml["userParam"]["timeAfterFixation"]

SEED = param_yaml["userParam"]["Seed"] + int(parallel) # will get loaded as int
iteration = param_yaml["userParam"]["iteration"]

def convert_fwdpy11params(rec_rate, mut_rate, alpha, sim_sample, sim_length):
    """
    Convert Nielsen param to fwdpy11 standard params
    2Ns => 2Ns
    4Nmu => mu*L
    4NLp => Lp
    """
    fwdpy_alpha = alpha #2Ns
    fwdpy_mut_rate = mut_rate/4/sim_sample*sim_length # fwdpy takes haploid genome mutation rate per generation
    fwdpy_rec_rate = rec_rate/4/sim_sample #genome wide recombination rate
    return fwdpy_alpha, fwdpy_mut_rate, fwdpy_rec_rate

def burnin_ancestry(sim_sample, sim_length, fwdpy_rec_rate, SEED):
    """
    For msprime, seed needs to be bigger than 1
    In msprime. recombination rates are per "unit" (bp) Lp => p
    The burnin population needs to be re-simulated for each iteration
    """
    initial_ts = msprime.sim_ancestry(
        samples = sim_sample, 
        population_size = sim_sample, 
        sequence_length=sim_length, 
        recombination_rate = fwdpy_rec_rate/sim_length, # p
        random_seed = SEED, 
        discrete_genome = True)
    pop = fwdpy11.DiploidPopulation.create_from_tskit(initial_ts)
    assert pop.N == sim_sample
    assert pop.tables.genome_length == sim_length
    return pop

def fwdpy_pdict(sim_length, fwdpy_rec_rate, sim_time):
    pdict = {
        # Here, the rec rate is number of events per region, Lp, not per bp!
        "recregions": [fwdpy11.PoissonInterval(0, int(sim_length), fwdpy_rec_rate, discrete=True)],
        "gvalue": fwdpy11.Multiplicative(2.0),
        "rates": (0, 0, None),
        "prune_selected": False,
        "simlen": sim_time,
    }
    params = fwdpy11.ModelParams(**pdict)
    return params

def write_trees_vcf(ts, nIndividual, prefix):
    ### output to tskit trees
    ts.dump(prefix  + ".trees")
    ### output to vcf
    with open(prefix +  ".vcf", "w") as vcf_file:
        ts.write_vcf(vcf_file, individuals =  [i for i in range(nIndividual)])

fwdpy_alpha, fwdpy_mut_rate, fwdpy_rec_rate = convert_fwdpy11params(rec_rate, mut_rate, alpha, sim_sample, sim_length)
params = fwdpy_pdict(sim_length, fwdpy_rec_rate, sim_time)
mut_pos_left = sim_length/2
mut_pos_right = sim_length/2 + 1


for i in range(0, iteration):
    seed = SEED + i 
    rng = fwdpy11.GSLrng(seed)
    pop = burnin_ancestry(sim_sample, sim_length, fwdpy_rec_rate, seed)
    ### drop the mutation exactly in the middle
    mutation_data = fwdpy11.conditional_models.NewMutationParameters(
        frequency=fwdpy11.conditional_models.AlleleCount(1),
        data=fwdpy11.NewMutationData(effect_size=fwdpy_alpha / 2 / pop.N, dominance=1),
        position=fwdpy11.conditional_models.PositionRange(left=mut_pos_left, right=mut_pos_right),
    )

    output = fwdpy11.conditional_models.selective_sweep(
        rng, 
        pop,
        params,
        mutation_data,
        fwdpy11.conditional_models.GlobalFixation(),
        return_when_stopping_condition_met = True #very important parameters, it will stop when sweep fixes
    )
    assert output.pop.generation <= params.simlen # miss leading
    assert pop.generation == 0

    ### not necessary for neutral after normalizing by polymorphism 
    #FIXATION_TIME = output.pop.fixation_times[0]
    nmuts = fwdpy11.infinite_sites(rng, output.pop, fwdpy_mut_rate)
    print(f"{nmuts} mutations added to sweep")
    sweep_prefix = os.path.join(savePath, "sweep", "sweep_" + str(seed))
    sweep_ts = output.pop.dump_tables_to_tskit()
    write_trees_vcf(sweep_ts, 50, sweep_prefix)
    
    ### using the smae fixation time to generate the neutral popualtion
    neutral_pop = burnin_ancestry(sim_sample, sim_length, fwdpy_rec_rate, seed)

    fwdpy11.evolvets(rng, neutral_pop, params, simplification_interval = 100)
    nmuts = fwdpy11.infinite_sites(rng, neutral_pop, fwdpy_mut_rate)
    print(f"{nmuts} mutations added to neutral")
    neutral_prefix = os.path.join(savePath, "neutral", "neutral_" + str(seed))
    neutral_ts = neutral_pop.dump_tables_to_tskit()
    write_trees_vcf(neutral_ts, 50, neutral_prefix)






