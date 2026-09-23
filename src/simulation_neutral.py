import fwdpy11 
import numpy as np 
import os
import msprime 
import copy
import argparse
import yaml


### Update: directly save the tree sequences, without converting to VCFs.
### Read tree sequences using tskit for moments.

parser = argparse.ArgumentParser(
                    prog='Neutral_simulations',
                   )
parser.add_argument('-y', '--yaml', dest = "yaml", required = True)
parser.add_argument('-s', '--savedir', dest = "savedir", required = True)
args = parser.parse_args()

### Adapted from fwdpy11 manual on selective sweep, https://molpopgen.github.io/fwdpy11/short_vignettes/selective_sweep.html




def neutral_simulation(pop, fwdpy11_params, mut_rate, sim_seed, mut_seed):
    pop_evolved = copy.deepcopy(pop)
    fwdpy11.evolvets(fwdpy11.GSLrng(sim_seed), pop_evolved, fwdpy11_params, simplification_interval = 10)
    pop_with_mut = copy.deepcopy(pop_evolved)
    nmuts = fwdpy11.infinite_sites(fwdpy11.GSLrng(mut_seed), pop_with_mut, mut_rate)
    print(f"{nmuts} mutations added to neutral simulations")
    return pop_with_mut



config_path = args.yaml
savePath = args.savedir
os.makedirs(savePath, exist_ok=True)


### set-up
with open(config_path) as f:
    config = yaml.safe_load(f)

scaling_factor = int(config["scaling"])
ne = int(config["model"]["ne"])
ne_scaled= int(ne / scaling_factor)
sim_region = int(config["model"]["genome"])
rec_rate = float(config["model"]["rec_rate"] * scaling_factor)
mut_rate = float(config["model"]["mut_rate"] * scaling_factor)
n_replicates = int(config["simulation"]["replicates"])
burnin_period = int(config["simulation"]["burnin"])
start_seed = int(config["simulation"]["seed"])
simlen = int(config["simulation"]["simlen"])

rng = np.random.RandomState(start_seed)
seeds = rng.randint(1, 2**31, size=(n_replicates, 3))


pdict = {
        "recregions": [fwdpy11.PoissonInterval(0, sim_region, sim_region*rec_rate, discrete=True)],
        "gvalue": fwdpy11.Multiplicative(2.0),
        "rates": (0, 0, None),
        "prune_selected": False,
        "simlen": simlen,
        "demography": fwdpy11.ForwardDemesGraph.tubes([ne_scaled], burnin=burnin_period)
}
params = fwdpy11.ModelParams(**pdict)




for sim_i in range(n_replicates):
    ancestry_seed, simulation_seed, mutation_seed = seeds[sim_i]
    initial_ts = msprime.sim_ancestry(
        samples=ne_scaled,
        population_size=ne_scaled,
        recombination_rate=rec_rate,
        random_seed=ancestry_seed,
        sequence_length=sim_region,
    )
    pop = fwdpy11.DiploidPopulation.create_from_tskit(initial_ts)
    
    neutral_pop = neutral_simulation(pop, params, mut_rate * sim_region, simulation_seed, mutation_seed)
    outfile = os.path.join(savePath, "_".join(["neutral", str(start_seed), "rep" + str(sim_i)]) + ".trees")
    neutral_ts = neutral_pop.dump_tables_to_tskit()
    neutral_ts.dump(outfile)
    