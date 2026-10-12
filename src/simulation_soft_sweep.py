import fwdpy11 
import numpy as np 
import fwdpy11.conditional_models
import os
import msprime 
import argparse
import yaml
from simulation_fwdpy11 import *



parser = argparse.ArgumentParser(
                    prog='fwdpy11_soft_sweep_simulation',
                   )
parser.add_argument('-s', '--savedir', dest = "savedir", required = True)
parser.add_argument('-y', '--yaml', dest = "yaml", required = True)
args = parser.parse_args()



config_path = args.yaml
savePath = args.savedir
os.makedirs(savePath, exist_ok=True)


### set-up
with open(config_path) as f:
    config = yaml.safe_load(f)

Ne = int(config["simulationParam"]["N"])
sim_region = int(config["simulationParam"]["genomeLength"])
rec_rate = float(config["simulationParam"]["recombinationRate"])
mut_rate = float(config["simulationParam"]["mutationRate"])
sel_coeff = float(config["simulationParam"]["selectionCoeff"])
master_seed = int(config["userParam"]["Seed"])
simlen = int(config["userParam"]["simulationTime"])
n_replicates = int(config["userParam"]["iteration"])
# soft_sweep specific
start_freq = config["simulationParam"]["startFrequency"]


mode = config["stop_condition"]["mode"]
stop_freq = config["stop_condition"]["targetFrequency"]
post_fix_gen = config["stop_condition"]["generationPostFix"]


pdict = {
    "recregions": [fwdpy11.PoissonInterval(0, sim_region, sim_region*rec_rate, discrete=True)],
    "gvalue": fwdpy11.Multiplicative(2.0),
    "rates": (0, 0, None),
    "prune_selected": False,
    "simlen": simlen,
    "demography": fwdpy11.ForwardDemesGraph.tubes([Ne], burnin=10)
}
params = fwdpy11.ModelParams(**pdict)




for freqs in start_freq:
    start_freq_min = float(freqs[0])
    start_freq_max = float(freqs[1])

    sweep_site = fwdpy11.conditional_models.NewMutationParameters(
        frequency=fwdpy11.conditional_models.FrequencyRange(start_freq_min, start_freq_max),
        data=fwdpy11.NewMutationData(effect_size = sel_coeff, dominance=1),
        position=fwdpy11.conditional_models.PositionRange(left=(sim_region/2-0.5), right=(sim_region/2 + 0.5)),
    )

    rng = np.random.RandomState(master_seed)
    seeds = rng.randint(1, 2**31, size=(n_replicates*2, 4))

    iteration = 0
    sim_i = 0
    while iteration < n_replicates:
        ancestry_seed, simulation_seed, mutation_seed, _ = seeds[sim_i]
        initial_ts = msprime.sim_ancestry(
            samples=Ne,
            population_size=Ne,
            recombination_rate=rec_rate,
            random_seed=ancestry_seed,
            sequence_length=sim_region,
        )
        pop = fwdpy11.DiploidPopulation.create_from_tskit(initial_ts)

        ### Check for initial_ts that have a suitable tree to add mutation at this location
        ### If an appropriate tree is in the center, given that the interval is narrow
        tree = initial_ts.at(sim_region/2, tracked_samples=[i for i in initial_ts.samples()])
        freqs = np.array([tree.num_tracked_samples(u) for u in tree.nodes()])/len(initial_ts.samples())
        if not any(start_freq_min <= f <= start_freq_max for f in freqs):
            ### not possible to add standing variation, move to next iteration
            sim_i = sim_i+1
            continue

        output = sweep_fixation(pop, params, sweep_site, simulation_seed)
        outfile = os.path.join(savePath, f"sweep_soft_freq_{start_freq_min}_{master_seed}_rep_{sim_i}_fixation.trees")
        save_tree_with_mutations(output, outfile, mut_rate * sim_region, mutation_seed)
        sim_i = sim_i+1
        iteration = iteration +1









