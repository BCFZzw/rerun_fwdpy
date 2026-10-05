import fwdpy11 
import numpy as np 
import fwdpy11.conditional_models
import os
import msprime 
import copy
import argparse
import yaml



parser = argparse.ArgumentParser(
                    prog='fwdpy11_sweep_simulation',
                   )
parser.add_argument('-o', '--output', dest = "output", required = True)
parser.add_argument('-y', '--yaml', dest = "yaml", required = True)
args = parser.parse_args()


class IncompleteSweep(object):
    def __init__(self, stop_frequency: float):
        if not 0 < stop_frequency < 1:
            raise ValueError(f"stop_frequency must be in (0, 1), got {stop_frequency}")
        self.stop_frequency = stop_frequency

    def __call__(
        self, pop: fwdpy11.DiploidPopulation, index: int, key: tuple
    ) -> fwdpy11.conditional_models.SimulationStatus:
        if pop.mutations[index].key != key:
            # it is fixed or lost, neither of 
            # which we want
            return fwdpy11.conditional_models.SimulationStatus.Restart
        if pop.mcounts[index] == 0:
            return fwdpy11.conditional_models.SimulationStatus.Restart
        
        # Terminate the first time we see the 
        # variant get about a freq above the frequency
        if pop.mcounts[index] / 2 / pop.N >= self.stop_frequency:
            return fwdpy11.conditional_models.SimulationStatus.Success
        # make sure there's a valid return value
        return fwdpy11.conditional_models.SimulationStatus.Continue


def sweep_fixation(pop, seed):
    output = fwdpy11.conditional_models.selective_sweep(
            fwdpy11.GSLrng(seed),
            pop,
            params,
            sweep_site,
            fwdpy11.conditional_models.GlobalFixation(),
            return_when_stopping_condition_met = True ### stopping when sweep fixed
    )
    assert output.pop.generation == output.pop.fixation_times[0]
    print(f"Evolved {fixation_time} generation till fixation")
    return output

def sweep_incomplete(pop, seed, stop_frequency):
    output = fwdpy11.conditional_models.selective_sweep(
            fwdpy11.GSLrng(seed),
            pop,
            params,
            sweep_site,
            IncompleteSweep(stop_frequency),
            return_when_stopping_condition_met = True ### stopping when sweep fixed
    )
    print(f"Evolved {output.pop.generation} generation till reaching frequency {stop_frequency}")
    return output

def sweep_post_fix(pop, seed1, seed2, params_post_fix):
    output = sweep_fixation(pop, seed1)
    fixation_time = output.pop.fixation_times[0]
    pop_fixed = copy.deepcopy(output.pop)
    post_fix_output = fwdpy11.evolvets(fwdpy11.GSLrng(seed2), pop_fixed, params_post_fix, simplification_interval = 100)
    print(f"Continue evolved post fixation for {post_fix_output.pop.generation - fixation_time} generation")
    return post_fix_output


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

Ne = int(config["simulationParam"]["N"])
sim_region = int(config["simulationParam"]["genomeLength"])
rec_rate = float(config["simulationParam"]["recombinationRate"])
mut_rate = float(config["simulationParam"]["mutationRate"])
sel_coeff = float(config["simulationParam"]["selectionCoeff"])
seed = int(config["userParam"]["Seed"])
simlen = int(config["userParam"]["simulationTime"])
n_replicates = int(config["userParam"]["iteration"])

sampling = int(config["sampling"]["nSample"])
output_sampling = False

if sampling != Ne:
    output_sampling = True

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


sweep_site = fwdpy11.conditional_models.NewMutationParameters(
    frequency=fwdpy11.conditional_models.AlleleCount(1),
    data=fwdpy11.NewMutationData(effect_size = sel_coeff, dominance=1),
    position=fwdpy11.conditional_models.PositionRange(left=(sim_region/2), right=(sim_region/2 + 1)),
)

if mode == "post_fixation":
    pdict_post_fix = copy.deepcopy(pdict)
    pdict_post_fix["simlen"] = post_fix_gen
    params_post_fix = fwdpy11.ModelParams(**pdict_post_fix)
else:
    params_post_fix = None


rng = np.random.RandomState(seed)
seeds = rng.randint(1, 2**31, size=(n_replicates, 3))


time_at_fixation = []

for sim_i in range(n_replicates):
    ancestry_seed, simulation_seed, mutation_seed = seeds[sim_i]
    initial_ts = msprime.sim_ancestry(
        samples=Ne,
        population_size=Ne,
        recombination_rate=rec_rate,
        random_seed=ancestry_seed,
        sequence_length=sim_region,
    )
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

    outfile = os.path.join(savePath, "_".join(["sweep", str(seed), "rep" + str(sim_i), "fixation"]) + ".trees")
    nmuts, pop_with_mut = add_neu_mutations(output.pop, mut_rate * sim_region, mutation_seed)
    print(f"{nmuts} mutations added to sweep at fixation")
    
    sweep_ts = pop_with_mut.dump_tables_to_tskit()
    sweep_ts.dump(outfile)
    
with open(os.path.join(savePath, "hard_sweep_fixation_time.txt"), "w") as f:
    f.write("\n".join(map(str, time_at_fixation)) + "\n")






