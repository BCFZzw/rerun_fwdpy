import fwdpy11 
import numpy as np 
import fwdpy11.conditional_models
import os
import msprime 
import copy
import argparse
import yaml



parser = argparse.ArgumentParser(
                    prog='Soft_sweep_simulations',
                   )
parser.add_argument('-s', '--savedir', dest = "savedir", required = True)
parser.add_argument('-y', '--yaml', dest = "yaml", required = True)
args = parser.parse_args()




def add_neu_mutations(pop, mut_rate, mut_seed):
    pop_with_mut = copy.deepcopy(pop)
    nmuts = fwdpy11.infinite_sites(fwdpy11.GSLrng(mut_seed), pop_with_mut, mut_rate)
    return nmuts, pop_with_mut



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
sel_coeff = float(config["selection"] * scaling_factor)
starting_freq1 = float(config["standing_var"]["left"])
starting_freq2 = float(config["standing_var"]["right"])
n_replicates = int(config["simulation"]["replicates"])
burnin_period = int(config["simulation"]["burnin"])
start_seed = int(config["simulation"]["seed"])
simlen = int(config["simulation"]["simlen"])

rng = np.random.RandomState(start_seed)
### For the soft-sweep case, some times no tree with the desired configuration is present in initial_ts
### Move on to the next set of seeds
seeds = rng.randint(1, 2**31, size=(n_replicates*2, 3))



pdict = {
        "recregions": [fwdpy11.PoissonInterval(0, sim_region, sim_region*rec_rate, discrete=True)],
        "gvalue": fwdpy11.Multiplicative(2.0),
        "rates": (0, 0, None),
        "prune_selected": False,
        "simlen": simlen,
        "demography": fwdpy11.ForwardDemesGraph.tubes([ne_scaled], burnin=burnin_period)
}
params = fwdpy11.ModelParams(**pdict)

sweep_site = fwdpy11.conditional_models.NewMutationParameters(
    frequency=fwdpy11.conditional_models.FrequencyRange(starting_freq1, starting_freq2),
    data=fwdpy11.NewMutationData(effect_size = sel_coeff, dominance=1),
    position=fwdpy11.conditional_models.PositionRange(left=(sim_region/2), right=(sim_region/2 + 1)),
)

time_at_fixation = []

skip_iteration = 0
for sim_i in range(n_replicates):
    ancestry_seed, simulation_seed, mutation_seed = seeds[sim_i+skip_iteration]
    initial_ts = msprime.sim_ancestry(
        samples=ne_scaled,
        population_size=ne_scaled,
        recombination_rate=rec_rate,
        random_seed=ancestry_seed,
        sequence_length=sim_region,
    )
    ### Check if an appropriate tree is on site
    tree = initial_ts.at(sim_region/2, tracked_samples=[i for i in initial_ts.samples()])
    freqs = np.array([tree.num_tracked_samples(u) for u in tree.nodes()])/len(initial_ts.samples())
    if not any(starting_freq1 <= f <= starting_freq2 for f in freqs):
        ### not possible to add standing variation, move to next iteration
        skip_iteration = skip_iteration+1
        continue

    pop = fwdpy11.DiploidPopulation.create_from_tskit(initial_ts)


    output = fwdpy11.conditional_models.selective_sweep(
            fwdpy11.GSLrng(simulation_seed),
            pop,
            params,
            sweep_site,
            fwdpy11.conditional_models.GlobalFixation(),
            return_when_stopping_condition_met = True ### stopping when sweep fixed
    )
    assert output.pop.generation == output.pop.fixation_times[0]
    fixation_time = output.pop.fixation_times[0]
    time_at_fixation.append(fixation_time)
    print(f"Evolved {fixation_time} generation till fixation")

    outfile = os.path.join(savePath, "_".join(["soft_sweep", str(start_seed), "rep" + str(sim_i), "fixation"]) + ".trees")
    nmuts, pop_with_mut = add_neu_mutations(output.pop, mut_rate * sim_region, mutation_seed)
    print(f"{nmuts} mutations added to soft sweep at fixation")
    
    sweep_ts = pop_with_mut.dump_tables_to_tskit()
    sweep_ts.dump(outfile)

    
with open(os.path.join(savePath, "soft_sweep_fixation_time.txt"), "w") as f:
    f.write("\n".join(map(str, time_at_fixation)) + "\n")






